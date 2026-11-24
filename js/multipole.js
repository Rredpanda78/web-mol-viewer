/* Web Mol Viewer — atomic multipoles (OpenMM AmoebaMultipoleForce conventions) and their electrostatic potential.
 *   mp = { d:    Float64Array(3N)  dipoles in the local frame (e Å),
 *          Q:    Float64Array(9N)  quadrupoles in the local frame (e Å², traceless, as stored by OpenMM / Tinker: Θ/3),
 *          ax:   Int32Array(3N)    z, x, y axis atoms (structure index, −1 = none),
 *          type: Uint8Array(N)     axis type (OpenMM enum) }
 *   MolMultipole.axisType(kz, kx, ky)      signs of the axis types in the force field -> axis type
 *   MolMultipole.lab(mp, X)                -> { d: Float64Array(3N), Q: Float64Array(6N: xx xy xz yy yz zz) } in the lab frame for coordinates X (3N, Å)
 *   MolMultipole.potential(x, y, z, q, d, Q, rmin) -> V(px, py, pz) in e/Å; d and Q are lab-frame arrays for those atoms (or null)
 *   MolMultipole.dipole(x, y, z, q, d)     -> [mx, my, mz] in e Å about the centre of the atoms (d may be null)
 */
'use strict';

const MolMultipole = (() => {
  const ZTHENX = 0, BISECTOR = 1, ZBISECT = 2, THREEFOLD = 3, ZONLY = 4, NOAXIS = 5;

  /* same rules as OpenMM's force-field reader: a negative axis type marks the bisector-style frames */
  function axisType(kz, kx, ky) {
    let t = ZTHENX;
    if (kz === 0) t = NOAXIS;
    if (kz !== 0 && kx === 0) t = ZONLY;
    if (kz < 0 || kx < 0) t = BISECTOR;
    if (kx < 0 && ky < 0) t = ZBISECT;
    if (kz < 0 && kx < 0 && ky < 0) t = THREEFOLD;
    return t;
  }

  const sub = (X, a, b) => [X[3 * a] - X[3 * b], X[3 * a + 1] - X[3 * b + 1], X[3 * a + 2] - X[3 * b + 2]];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const unit = u => { const n = Math.sqrt(dot(u, u)) || 1; return [u[0] / n, u[1] / n, u[2] / n]; };
  const add = (u, v) => [u[0] + v[0], u[1] + v[1], u[2] + v[2]];

  /* rotate every atom's dipole and quadrupole from its local frame to the lab frame */
  function lab(mp, X) {
    const N = mp.type.length; const d = new Float64Array(3 * N), Q = new Float64Array(6 * N);
    for (let i = 0; i < N; i++) {
      const o = 3 * i, p = 9 * i; const t = mp.type[i];
      let l0 = mp.d[o], l1 = mp.d[o + 1], l2 = mp.d[o + 2];
      const q = [mp.Q[p], mp.Q[p + 1], mp.Q[p + 2], mp.Q[p + 3], mp.Q[p + 4], mp.Q[p + 5], mp.Q[p + 6], mp.Q[p + 7], mp.Q[p + 8]];
      if (!l0 && !l1 && !l2 && !q.some(v => v)) continue;
      const iz = mp.ax[o], ix = mp.ax[o + 1], iy = mp.ax[o + 2];
      let R;
      if (t === NOAXIS || iz < 0) R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      else {
        let vz = unit(sub(X, iz, i)), vx, vy;
        if (t === ZONLY || ix < 0) vx = Math.abs(vz[0]) < 0.866 ? [1, 0, 0] : [0, 1, 0];
        else {
          vx = sub(X, ix, i);
          if (t === BISECTOR) { vx = unit(vx); vz = unit(add(vz, vx)); }
          else if (t === ZBISECT && iy >= 0) { vx = unit(vx); vy = unit(sub(X, iy, i)); vx = unit(add(vx, vy)); }
          else if (t === THREEFOLD && iy >= 0) { vx = unit(vx); vy = unit(sub(X, iy, i)); vz = unit(add(vz, add(vx, vy))); }
        }
        const c = dot(vz, vx); vx = unit([vx[0] - c * vz[0], vx[1] - c * vz[1], vx[2] - c * vz[2]]); vy = cross(vz, vx);
        R = [vx, vy, vz];
        if (t === ZTHENX && iy >= 0) {                              // chiral centre with the other handedness: mirror the y axis
          const vol = dot(cross(sub(X, iz, iy), sub(X, ix, iy)), sub(X, i, iy));
          if (vol < 0) { l1 = -l1; q[1] = -q[1]; q[3] = -q[3]; q[5] = -q[5]; q[7] = -q[7]; }
        }
      }
      const L = [l0, l1, l2];
      for (let a = 0; a < 3; a++) d[o + a] = L[0] * R[0][a] + L[1] * R[1][a] + L[2] * R[2][a];
      const G = (a, b) => { let s = 0; for (let k = 0; k < 3; k++) for (let m = 0; m < 3; m++) s += R[k][a] * R[m][b] * q[3 * k + m]; return s; };
      const u = 6 * i; Q[u] = G(0, 0); Q[u + 1] = G(0, 1); Q[u + 2] = G(0, 2); Q[u + 3] = G(1, 1); Q[u + 4] = G(1, 2); Q[u + 5] = G(2, 2);
    }
    return { d, Q };
  }

  /* potential of charges + dipoles + quadrupoles: q/r + μ·R/r³ + Σ Q_ab (3 R_a R_b − r² δ_ab)/r⁵, R = point − atom */
  function potential(x, y, z, q, d, Q, rmin) {
    const n = x.length; const lo = rmin || 0;
    if (!d && !Q) return (px, py, pz) => { let s = 0; for (let k = 0; k < n; k++) { const dx = px - x[k], dy = py - y[k], dz = pz - z[k]; const r = Math.sqrt(dx * dx + dy * dy + dz * dz); s += q[k] / (r < lo ? lo : r); } return s; };
    return (px, py, pz) => {
      let s = 0;
      for (let k = 0; k < n; k++) {
        const dx = px - x[k], dy = py - y[k], dz = pz - z[k]; let r = Math.sqrt(dx * dx + dy * dy + dz * dz); if (r < lo) r = lo;
        const i1 = 1 / r, i3 = i1 * i1 * i1; s += q[k] * i1;
        if (d) { const o = 3 * k; s += (d[o] * dx + d[o + 1] * dy + d[o + 2] * dz) * i3; }
        if (Q) {
          const u = 6 * k; const xx = Q[u], xy = Q[u + 1], xz = Q[u + 2], yy = Q[u + 3], yz = Q[u + 4], zz = Q[u + 5];
          const rqr = xx * dx * dx + yy * dy * dy + zz * dz * dz + 2 * (xy * dx * dy + xz * dx * dz + yz * dy * dz);
          s += (3 * rqr * i1 * i1 - (xx + yy + zz)) * i3;
        }
      }
      return s;
    };
  }

  /* dipole moment about the centre of the atoms (e Å); meaningful for a neutral molecule */
  function dipole(x, y, z, q, d) {
    const n = x.length; let cx = 0, cy = 0, cz = 0; for (let k = 0; k < n; k++) { cx += x[k]; cy += y[k]; cz += z[k]; }
    cx /= n || 1; cy /= n || 1; cz /= n || 1; const m = [0, 0, 0];
    for (let k = 0; k < n; k++) {
      m[0] += q[k] * (x[k] - cx); m[1] += q[k] * (y[k] - cy); m[2] += q[k] * (z[k] - cz);
      if (d) { m[0] += d[3 * k]; m[1] += d[3 * k + 1]; m[2] += d[3 * k + 2]; }
    }
    return m;
  }

  return { axisType, lab, potential, dipole, NOAXIS };
})();
