/* Web Mol Viewer — UI and viewer logic (3Dmol.js does the drawing). */
'use strict';

(() => {
  const $ = id => document.getElementById(id);
  const viewer = $3Dmol.createViewer($('viewer'), { backgroundColor: 'white' });
  const files = [];            // entries: { name, kind, frames, i, measures, sets, traj, xyz, energy, scan, topo }
  let cur = null;              // active entry
  let atoms = [];              // atoms of the frame on screen (all of them)
  let vis = [];                // the atoms that are drawn (solvent filter, no virtual sites)
  let picked = [];             // atom indices picked so far for the measurement being made
  let mode = 'distance';       // measuring tool: distance (2 atoms), angle (3), dihedral (4)
  let timer = null;
  let plotMap = null;          // how the plot maps frames to x (set by drawPlot, used by clicks)
  const NEED = { distance: 2, angle: 3, dihedral: 4 };
  const MCOL = { distance: '#1f4e79', angle: '#2e7d32', dihedral: '#7030a0' };
  const INTRA = '#d62728', INTER = '#1f77b4';
  const HARTREE = 2625.4996;   // kJ/mol

  const STYLE = {
    ballstick: { stick: { radius: 0.14, colorscheme: 'Jmol' }, sphere: { scale: 0.27, colorscheme: 'Jmol' } },
    stick: { stick: { radius: 0.18, colorscheme: 'Jmol' } },
    line: { line: { linewidth: 2, colorscheme: 'Jmol' } },
    sphere: { sphere: { colorscheme: 'Jmol' } },
  };
  const WATER = { ballstick: { stick: { radius: 0.07, colorscheme: 'Jmol' }, sphere: { scale: 0.16, colorscheme: 'Jmol' } },
    stick: { stick: { radius: 0.08, colorscheme: 'Jmol' } }, line: { line: { linewidth: 1, colorscheme: 'Jmol' } }, sphere: { sphere: { scale: 0.6, colorscheme: 'Jmol' } } };

  /* ---------- geometry helpers ---------- */
  const sub = (a, b) => [a.x - b.x, a.y - b.y, a.z - b.z];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const norm = u => Math.sqrt(dot(u, u));
  const dist = (a, b) => norm(sub(a, b));
  const angle = (a, b, c) => { const u = sub(a, b), v = sub(c, b); return Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (norm(u) * norm(v))))) * 180 / Math.PI; };
  function dihedral(a, b, c, d) {
    const b0 = sub(a, b), b1 = sub(c, b), b2 = sub(d, c); const n1 = norm(b1); const u = b1.map(x => x / n1);
    const v = b0.map((x, i) => x - dot(b0, u) * u[i]), w = b2.map((x, i) => x - dot(b2, u) * u[i]);
    return Math.atan2(dot(cross(u, v), w), dot(v, w)) * 180 / Math.PI;
  }
  const mid = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 });
  const tag = a => a.elem + (a.index + 1);
  const dark = () => $('dark').checked;
  const lsize = () => Math.max(6, Math.min(40, parseFloat($('label-size').value) || 14));      // atom / charge labels; distances are 2 pt larger
  const nF = e => (e ? (e.traj ? e.traj.nframes : e.frames.length) : 0);
  const signed = (v, d) => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(d);

  /* ---------- topology of an entry (built once, from the first frame) ---------- */
  const WATER_RE = /^(HOH|WAT|SOL|H2O|DOD|TIP\d?P?|T[345]P|SPCE?|OPC3?)$/i;
  function buildTopo(e, A) {
    const n = A.length; const virt = new Uint8Array(n), water = new Uint8Array(n), res = new Int32Array(n);
    const elem = A.map(a => a.elem || 'X');
    let r = -1, key = null;
    A.forEach((a, i) => {
      const nm = (a.atom || '').trim();
      virt[i] = /^(EP|LP|M|X|DU|Q)$/i.test(elem[i]) || /^(MW|EPW|EP\d*|LP\d*)$/i.test(nm) ? 1 : 0;
      water[i] = WATER_RE.test((a.resn || '').trim()) ? 1 : 0;
      const k = (a.chain || '') + '|' + a.resi + '|' + (a.resn || ''); if (k !== key) { key = k; r++; } res[i] = r;
    });
    const nWater = water.reduce((s, v) => s + v, 0); const hasSolvent = nWater > 0 && nWater < n;
    const bonds = A.map((a, i) => (a.bonds || []).filter(j => !virt[j] && !virt[i]));
    const mol = new Int32Array(n).fill(-1); let m = 0;            // molecules = connected components
    for (let i = 0; i < n; i++) {
      if (mol[i] >= 0) continue; const q = [i]; mol[i] = m;
      while (q.length) { const u = q.pop(); bonds[u].forEach(v => { if (mol[v] < 0) { mol[v] = m; q.push(v); } }); }
      m++;
    }
    const sol = new Uint8Array(n); for (let i = 0; i < n; i++) sol[i] = (!hasSolvent || !water[i]) && !virt[i] ? 1 : 0;
    const resAtoms = []; for (let i = 0; i < n; i++) (resAtoms[res[i]] = resAtoms[res[i]] || []).push(i);
    const hd = [];                                                // polar H: { h, d (donor N/O), skip (atoms within 2 bonds of the donor) }
    for (let i = 0; i < n; i++) {
      if (elem[i] !== 'H' || virt[i]) continue;
      const d = bonds[i].find(j => elem[j] === 'N' || elem[j] === 'O'); if (d === undefined) continue;
      const skip = new Set([d]); let layer = [d];
      for (let k = 0; k < 2; k++) { const nxt = []; layer.forEach(u => bonds[u].forEach(v => { if (!skip.has(v)) { skip.add(v); nxt.push(v); } })); layer = nxt; }
      hd.push({ h: i, d, skip });
    }
    const acc = []; for (let i = 0; i < n; i++) if ((elem[i] === 'N' || elem[i] === 'O') && !virt[i]) acc.push(i);
    e.topo = { n, elem, virt, water, res, resAtoms, bonds, mol, sol, hd, acc, hasSolvent,
      solute: Array.from(sol.keys()).filter(i => sol[i]), accSol: acc.filter(i => sol[i]) };
    const native = A.map(a => (a.properties && typeof a.properties.partialCharge === 'number' ? a.properties.partialCharge : null));   // mol2 / pqr charges
    if (native.some(v => v !== null && v !== 0) && !e.sets.length) e.sets.push({ name: e.name, q: native, note: 'Charges from the structure file' });
    return e.topo;
  }

  /* put the frame in a convenient place: solute centred at the origin, each water molecule in its nearest periodic image */
  function tidy(T, X) {
    const { x, y, z } = X; const S = T.solute; if (!S.length) return X;
    let cx = 0, cy = 0, cz = 0; S.forEach(i => { cx += x[i]; cy += y[i]; cz += z[i]; }); cx /= S.length; cy /= S.length; cz /= S.length;
    const c = X.cell && X.cell.ortho && X.cell.a > 0 ? X.cell : null;
    if (c && T.hasSolvent) {
      T.resAtoms.forEach(R => {
        const i0 = R[0]; if (T.sol[i0]) return;
        const sx = -Math.round((x[i0] - cx) / c.a) * c.a, sy = -Math.round((y[i0] - cy) / c.b) * c.b, sz = -Math.round((z[i0] - cz) / c.c) * c.c;
        if (sx || sy || sz) R.forEach(i => { x[i] += sx; y[i] += sy; z[i] += sz; });
      });
    }
    for (let i = 0; i < x.length; i++) { x[i] -= cx; y[i] -= cy; z[i] -= cz; }
    return X;
  }

  /* coordinates of frame i of an entry as typed arrays (null when that frame cannot be read as plain coordinates) */
  async function coordsOf(e, i) {
    if (e.traj) return tidy(e.topo, await MolTraj.readFrame(e.traj, i));
    const f = e.frames[i]; if (f.format !== 'xyz') return null;
    const L = f.data.split('\n'); const n = parseInt(L[0], 10); const x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n);
    for (let k = 0; k < n; k++) { const p = L[k + 2].trim().split(/\s+/); x[k] = +p[1]; y[k] = +p[2]; z[k] = +p[3]; }
    return { x, y, z, cell: null };
  }

  /* polar contacts H···N/O of one frame; returns [{ h, d, y, r, inter }] (indices) */
  function contactsOf(T, X, cut, only) {
    const { x, y, z } = X; const out = []; const c2 = cut * cut;
    const ok = only ? i => only[i] : () => true;
    for (const p of T.hd) {
      if (!ok(p.h)) continue;
      const list = T.hasSolvent && !T.sol[p.h] ? T.accSol : T.acc;    // in a solvated system only contacts that involve the solute
      for (const a of list) {
        if (p.skip.has(a) || !ok(a)) continue;
        const dx = x[p.h] - x[a], dy = y[p.h] - y[a], dz = z[p.h] - z[a]; const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 <= c2) out.push({ h: p.h, d: p.d, y: a, r: Math.sqrt(r2), inter: T.mol[p.h] !== T.mol[a] });
      }
    }
    return out.sort((p, q) => p.r - q.r);
  }

  /* ---------- charges ---------- */
  const allXYZ = () => { const X = new Float64Array(3 * atoms.length); atoms.forEach((a, i) => { X[3 * i] = a.x; X[3 * i + 1] = a.y; X[3 * i + 2] = a.z; }); return X; };
  function activeQ() {
    if (!cur || !cur.sets.length) return null;
    const v = $('q-set').value;
    if (v.startsWith('d')) { const k = +v.slice(1), a = cur.sets[0].q, b = cur.sets[k].q; return b.map((q, i) => (q === null || a[i] === null ? null : q - a[i])); }
    return (cur.sets[+v] || cur.sets[0]).q;
  }
  const qShow = () => (cur && cur.sets.length ? $('q-show').value : 'off');
  function qScale(Q) {
    const u = parseFloat($('q-scale').value); if (u > 0) return u;
    let m = 0; vis.forEach(a => { const q = Q[a.index]; if (q !== null && q !== undefined && Math.abs(q) > m) m = Math.abs(q); });
    return Math.max(0.05, m);
  }
  function qColour(q, s) {                                        // red (negative) – white – blue (positive)
    if (q === null || q === undefined) return 0xb0b0b0;
    const t = Math.max(-1, Math.min(1, q / s)); const w = Math.round(255 * (1 - Math.abs(t)));
    return t < 0 ? (255 << 16) | (w << 8) | w : (w << 16) | (w << 8) | 255;
  }
  function updateChargePanel() {
    const has = !!(cur && cur.sets.length); $('q-rows').hidden = !has;
    if (!has) { $('q-info').textContent = 'No charges. Drop a charge file while the structure is selected.'; return; }
    const sel = $('q-set'); const keep = sel.value; sel.innerHTML = '';
    cur.sets.forEach((s, k) => sel.add(new Option(s.name, String(k))));
    for (let k = 1; k < cur.sets.length; k++) sel.add(new Option(`Δ: ${cur.sets[k].name} − ${cur.sets[0].name}`, 'd' + k));
    sel.value = Array.from(sel.options).some(o => o.value === keep) ? keep : String(cur.sets.length - 1);
    const Q = activeQ(); const T = cur.topo; let n = 0, s = 0, sa = 0;
    const idx = T ? T.solute : Q.map((_, i) => i);
    idx.forEach(i => { const q = Q[i]; if (q !== null && q !== undefined) { n++; s += q; sa += Math.abs(q); } });
    const note = sel.value.startsWith('d') ? 'Difference between two sets' : (cur.sets[+sel.value].note || '');
    let dip = '';                                                 // dipole moment of a neutral molecule (charges, plus the atomic dipoles of a multipole set)
    if (!sel.value.startsWith('d') && n > 1 && Math.abs(s) < 0.02) {
      const ix = idx.filter(i => Q[i] !== null && Q[i] !== undefined); const S = cur.sets[+sel.value];
      const fx = Float64Array.from(ix, i => atoms[i].x), fy = Float64Array.from(ix, i => atoms[i].y), fz = Float64Array.from(ix, i => atoms[i].z), fq = Float64Array.from(ix, i => Q[i]);
      const D = m => (4.80321 * Math.hypot(m[0], m[1], m[2])).toFixed(2); const m0 = MolMultipole.dipole(fx, fy, fz, fq, null);
      if (S && S.mp) {
        const L = MolMultipole.lab(S.mp, allXYZ()); const fd = new Float64Array(3 * ix.length); ix.forEach((i, k) => { for (let c = 0; c < 3; c++) fd[3 * k + c] = L.d[3 * i + c]; });
        dip = `\n|μ| = ${D(MolMultipole.dipole(fx, fy, fz, fq, fd))} D with the atomic dipoles (charges only: ${D(m0)} D)`;
      } else dip = `\n|μ| = ${D(m0)} D`;
    }
    $('q-info').textContent = `Σq = ${signed(s, 3)} e   Σ|q| = ${sa.toFixed(3)} e   (${n} ${T && T.hasSolvent ? 'solute ' : ''}atoms)` + dip + (note ? '\n' + note : '');
    const sc = qScale(Q); $('q-legend').textContent = qShow() === 'colour' || qShow() === 'both' ? `red −${sc.toFixed(2)} … white 0 … blue +${sc.toFixed(2)} e` : '';
  }
  function chargeTable() {
    if (!cur || !cur.sets.length) return '';
    const rows = [['atom', 'element'].concat(cur.sets.map(s => s.name)).join('\t')];
    atoms.forEach((a, i) => { if (cur.sets.some(s => s.q[i] !== null && s.q[i] !== undefined)) rows.push([i + 1, a.elem].concat(cur.sets.map(s => (s.q[i] === null || s.q[i] === undefined ? '' : s.q[i].toFixed(5)))).join('\t')); });
    return rows.join('\n');
  }

  /* ---------- electrostatic potential on the molecular surface (atomic charges, or the full atomic multipoles) ---------- */
  const VDW = { H: 1.2, C: 1.7, N: 1.55, O: 1.52, F: 1.47, P: 1.8, S: 1.8, Cl: 1.75, Br: 1.85 };
  const COULOMB = { kcal: 332.0637, kJ: 1389.354 };                // energy of a +1 e probe charge, per (e / Å)
  let surfTok = 0;
  function drawSurface() {
    const tok = ++surfTok; const leg = $('esp-legend'); leg.hidden = true; leg.innerHTML = '';
    const on = !!cur && $('surf').value !== 'off'; $('surf-rows').hidden = !on; $('surf-info').textContent = '';
    if (!on) return;
    const Q = activeQ(); const T = cur.topo;
    if (!Q) { $('surf-info').textContent = 'The surface is coloured by the potential of the atomic charges: load charges first (see Charges).'; return; }
    const idx = T.solute.filter(i => Q[i] !== null && Q[i] !== undefined); const n = idx.length;
    if (!n) return;
    if (n > 500) { $('surf-info').textContent = 'ESP surface is limited to 500 charged atoms.'; return; }
    const unit = $('surf-unit').value; const K = COULOMB[unit]; const type = $('surf').value; const probe = type === 'SAS' ? 1.4 : 0;
    const ax = new Float64Array(n), ay = new Float64Array(n), az = new Float64Array(n), aq = new Float64Array(n), ar = new Float64Array(n);
    idx.forEach((i, k) => { const a = atoms[i]; ax[k] = a.x; ay[k] = a.y; az[k] = a.z; aq[k] = Q[i]; ar[k] = (VDW[a.elem] || 1.7) + probe; });
    const sel = $('q-set').value;                                    // a multipole set brings dipoles and quadrupoles; a difference of two sets subtracts them too
    const comps = sel.startsWith('d') ? [[cur.sets[+sel.slice(1)], 1], [cur.sets[0], -1]] : [[cur.sets[+sel] || cur.sets[0], 1]];
    const hasMP = comps.some(c => c[0] && c[0].mp); $('surf-src-row').hidden = !hasMP; const full = hasMP && $('surf-src').value !== 'mono';
    let aD = null, aQ = null;
    if (full) {
      const X = allXYZ(); aD = new Float64Array(3 * n); aQ = new Float64Array(6 * n);
      comps.forEach(([S, sg]) => {
        if (!S || !S.mp) return; const L = MolMultipole.lab(S.mp, X);
        idx.forEach((i, k) => { for (let c = 0; c < 3; c++) aD[3 * k + c] += sg * L.d[3 * i + c]; for (let c = 0; c < 6; c++) aQ[6 * k + c] += sg * L.Q[6 * i + c]; });
      });
    }
    const V0 = MolMultipole.potential(ax, ay, az, aq, aD, aQ, 0.6); const V = (x, y, z) => K * V0(x, y, z);
    let lo = Infinity, hi = -Infinity; const NP = 90;               // range: the potential on the surface itself (points on each atomic sphere that are outside the others)
    for (let k = 0; k < n; k++) for (let p = 0; p < NP; p++) {
      const ph = Math.acos(1 - 2 * (p + 0.5) / NP), th = Math.PI * (1 + Math.sqrt(5)) * p;
      const x = ax[k] + ar[k] * Math.sin(ph) * Math.cos(th), y = ay[k] + ar[k] * Math.sin(ph) * Math.sin(th), z = az[k] + ar[k] * Math.cos(ph);
      let inside = false;
      for (let j = 0; j < n && !inside; j++) { if (j === k) continue; const dx = x - ax[j], dy = y - ay[j], dz = z - az[j]; if (dx * dx + dy * dy + dz * dz < ar[j] * ar[j]) inside = true; }
      if (inside) continue; const v = V(x, y, z); if (v < lo) lo = v; if (v > hi) hi = v;
    }
    const rwb = $('surf-grad').value === 'rwb';
    if (rwb) { const m = Math.max(Math.abs(lo), Math.abs(hi)); lo = -m; hi = m; }      // white = 0
    const umin = parseFloat($('surf-min').value), umax = parseFloat($('surf-max').value);
    if (!isNaN(umin)) lo = umin; if (!isNaN(umax)) hi = umax; if (!(hi > lo)) hi = lo + 1;
    const mg = 2.4 + probe; let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let k = 0; k < n; k++) { x0 = Math.min(x0, ax[k]); x1 = Math.max(x1, ax[k]); y0 = Math.min(y0, ay[k]); y1 = Math.max(y1, ay[k]); z0 = Math.min(z0, az[k]); z1 = Math.max(z1, az[k]); }
    x0 -= mg; y0 -= mg; z0 -= mg; x1 += mg; y1 += mg; z1 += mg;
    let h = 0.2, nx, ny, nz;                                         // 3Dmol reads the nearest grid point, so the grid has to be fine
    for (;;) { nx = Math.ceil((x1 - x0) / h) + 1; ny = Math.ceil((y1 - y0) / h) + 1; nz = Math.ceil((z1 - z0) / h) + 1; if (nx * ny * nz * n <= (full ? 2e7 : 4e7)) break; h *= 1.12; }
    const data = new Float32Array(nx * ny * nz);
    for (let ix = 0; ix < nx; ix++) for (let iy = 0; iy < ny; iy++) { const o = (ix * ny + iy) * nz; const x = x0 + ix * h, y = y0 + iy * h; for (let iz = 0; iz < nz; iz++) data[o + iz] = V(x, y, z0 + iz * h); }
    const vd = new $3Dmol.VolumeData('', 'none'); vd.origin = { x: x0, y: y0, z: z0 }; vd.unit = { x: h, y: h, z: h }; vd.size = { x: nx, y: ny, z: nz }; vd.data = data;
    const grad = rwb ? new $3Dmol.Gradient.RWB(lo, hi) : new $3Dmol.Gradient.ROYGB(lo, hi);
    const pr = viewer.addSurface($3Dmol.SurfaceType[type], { opacity: parseFloat($('surf-op').value) || 0.8, voldata: vd, volscheme: grad }, { index: idx });
    Promise.resolve(pr).then(id => { if (tok !== surfTok) { try { viewer.removeSurface(id); } catch (e) { /* already gone */ } } viewer.render(); });
    const ttl = document.createElement('div'); ttl.className = 'ttl'; ttl.textContent = `ESP (${unit === 'kJ' ? 'kJ/mol' : 'kcal/mol'})`; leg.append(ttl);
    for (let k = 0; k <= 10; k++) {                                  // legend: positive (blue) on top, as in GaussView-style pictures
      const v = hi - (hi - lo) * k / 10; const row = document.createElement('div'); row.className = 'rowl';
      const sw = document.createElement('span'); sw.className = 'sw'; const col = '#' + grad.valueToHex(v).toString(16).padStart(6, '0'); sw.style.background = col;
      const tx = document.createElement('span'); tx.className = 'v'; tx.textContent = (v < 0 ? '−' : '') + Math.abs(v).toFixed(1); row.append(sw, tx); leg.append(row);
    }
    leg.hidden = false; leg.dataset.unit = unit;
    const note = sel.startsWith('d') ? '' : ((cur.sets[+sel] || {}).note || '');
    $('surf-info').textContent = (full ? 'Potential of the atomic multipoles (charges + dipoles + quadrupoles, rotated from their local frames; energy of a +1 e probe)' : 'Potential of the atomic point charges (energy of a +1 e probe)')
      + ' on the ' + (type === 'SAS' ? 'solvent-accessible' : 'van der Waals') + ' surface.'
      + (sel.startsWith('d') ? ' Showing the difference between two charge sets.' : '')
      + (!full && (hasMP || /multipole|monopole/i.test(note)) ? ' Multipole model: only the monopoles are used — atomic dipoles and quadrupoles are not included, so this is not the full model potential.' : '');
  }

  /* picture of the view with the ESP legend drawn in */
  function savePNG() {
    if (!cur) return; const name = cur.name.replace(/\.[^.]+$/, '') + '.png'; const uri = viewer.pngURI(); const leg = $('esp-legend');
    if (leg.hidden) { download(name, uri); return; }
    const img = new Image();
    img.onload = () => {
      const cv = document.createElement('canvas'); cv.width = img.width; cv.height = img.height; const g = cv.getContext('2d'); g.drawImage(img, 0, 0);
      const k = img.width / $('viewer').clientWidth; const rows = Array.from(leg.querySelectorAll('.rowl')); const x = 14 * k, rh = 17 * k; let y = 14 * k;
      g.fillStyle = 'rgba(255,255,255,0.85)'; g.fillRect(8 * k, 8 * k, 110 * k, (26 + rows.length * 17) * k);
      g.textBaseline = 'middle'; g.fillStyle = '#1f283d'; g.font = `600 ${13 * k}px "Segoe UI", Arial, sans-serif`; g.fillText(leg.querySelector('.ttl').textContent, x, y + 6 * k); y += 18 * k;
      g.font = `700 ${12 * k}px Consolas, "Courier New", monospace`;
      rows.forEach(r => { g.fillStyle = r.querySelector('.sw').style.background; g.fillRect(x, y, 18 * k, rh); g.fillStyle = '#1f283d'; g.textAlign = 'right'; g.fillText(r.querySelector('.v').textContent, x + 96 * k, y + rh / 2); g.textAlign = 'left'; y += rh; });
      download(name, cv.toDataURL('image/png'));
    };
    img.src = uri;
  }

  /* ---------- drawing ---------- */
  function styleSpec(base, Q, sc) {
    const out = {};
    Object.keys(base).forEach(k => { out[k] = Object.assign({}, base[k]); if (Q) { delete out[k].colorscheme; out[k].colorfunc = a => qColour(Q[a.index], sc); } });
    return out;
  }
  function applyStyle() {
    const s = $('style').value; const T = cur.topo; const show = qShow(); const Q = show === 'colour' || show === 'both' ? activeQ() : null; const sc = Q ? qScale(Q) : 1;
    const all = vis.length === atoms.length; const sel = all ? {} : { index: vis.map(a => a.index) };
    viewer.setStyle({}, {});
    if (s === 'cartoon') {
      viewer.setStyle(sel, { cartoon: { color: 'spectrum' } });
      viewer.setStyle(Object.assign({ hetflag: true }, sel), styleSpec(STYLE.ballstick, Q, sc));
    } else {
      viewer.setStyle(sel, styleSpec(STYLE[s], Q, sc));
      if (T.hasSolvent) { const w = vis.filter(a => T.water[a.index]).map(a => a.index); if (w.length) viewer.setStyle({ index: w }, styleSpec(WATER[s], Q, sc)); }
    }
    if ($('hide-h').checked) viewer.setStyle({ elem: 'H' }, {});
  }

  function label(text, pos, color, size, bg) {
    viewer.addLabel(text, { position: pos, fontSize: size || lsize() + 2, fontColor: color || (dark() ? 'white' : 'black'), showBackground: bg !== false,
      backgroundColor: dark() ? 'black' : 'white', backgroundOpacity: 0.65, inFront: true, alignment: 'center' });
  }
  function dashed(a, b, color) {
    viewer.addCylinder({ start: { x: a.x, y: a.y, z: a.z }, end: { x: b.x, y: b.y, z: b.z }, radius: 0.045, dashed: true, color, fromCap: 1, toCap: 1 });
  }

  const JMOL = { H: 0xffffff, C: 0x909090, N: 0x3050f8, O: 0xff0d0d, F: 0x90e050, P: 0xff8000, S: 0xffff30, Cl: 0x1ff01f, Br: 0xa62929 };
  const lum = c => (0.299 * ((c >> 16) & 255) + 0.587 * ((c >> 8) & 255) + 0.114 * (c & 255)) / 255;
  function drawLabels() {
    const m = $('labels').value; const show = qShow(); const Q = show === 'labels' || show === 'both' ? activeQ() : null;
    if (m === 'none' && !Q) return;
    const hideH = $('hide-h').checked;
    if (vis.length > 400) { $('hover').textContent = 'Labels are switched off above 400 atoms shown.'; return; }
    const bg = $('label-bg').checked; const QC = show === 'colour' || show === 'both' ? activeQ() : null; const sc = QC ? qScale(QC) : 1;      // without a background the text is white on dark atoms, black on light ones
    vis.forEach(a => {
      if (hideH && a.elem === 'H') return;
      const parts = [];
      if (m === 'index') parts.push(String(a.index + 1)); else if (m === 'both') parts.push(tag(a)); else if (m === 'elem' && a.elem !== 'C' && a.elem !== 'H') parts.push(a.elem);
      if (Q && Q[a.index] !== null && Q[a.index] !== undefined) parts.push(signed(Q[a.index], 3));
      if (parts.length) label(parts.join(' '), a, bg ? null : (lum(QC ? qColour(QC[a.index], sc) : (JMOL[a.elem] || 0x909090)) < 0.55 ? 'white' : 'black'), lsize(), bg);
    });
  }

  function frameX() {                                               // coordinates of the frame on screen as arrays
    const n = atoms.length; const x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n);
    atoms.forEach((a, i) => { x[i] = a.x; y[i] = a.y; z[i] = a.z; }); return { x, y, z, cell: null };
  }
  function drawContacts() {
    const box = $('contact-list');
    if (!$('contacts').checked || !atoms.length) { box.textContent = ''; return; }
    const cut = parseFloat($('cutoff').value) || 2.5; const shown = new Uint8Array(atoms.length); vis.forEach(a => { shown[a.index] = 1; });
    const C = contactsOf(cur.topo, frameX(), cut, shown); const draw = C.slice(0, 150);
    draw.forEach(c => { const col = c.inter ? INTER : INTRA; dashed(atoms[c.h], atoms[c.y], col); label(c.r.toFixed(2), mid(atoms[c.h], atoms[c.y]), col); });
    const ni = C.filter(c => c.inter).length;
    const head = C.length ? `${ni} between molecules (blue), ${C.length - ni} intramolecular (red)\n` : '';
    box.textContent = C.length ? head + draw.map(c => `${c.inter ? 'inter' : 'intra'} ${tag(atoms[c.d])}–${tag(atoms[c.h])}···${tag(atoms[c.y])}  ${c.r.toFixed(2)} Å  ${angle(atoms[c.d], atoms[c.h], atoms[c.y]).toFixed(0)}°`).join('\n') : `No contact below ${cut.toFixed(1)} Å.`;
  }

  /* value of one measurement; P(i) gives the position {x, y, z} of atom i (or undefined) */
  function measureValue(m, P) {
    const A = m.idx.map(P); if (A.some(p => !p)) return null;
    return m.type === 'distance' ? dist(A[0], A[1]) : m.type === 'angle' ? angle(A[0], A[1], A[2]) : dihedral(A[0], A[1], A[2], A[3]);
  }
  const fmtValue = (m, v) => v === null ? '—' : m.type === 'distance' ? v.toFixed(3) + ' Å' : v.toFixed(1) + '°';
  const measureName = m => m.idx.map(i => atoms[i] ? tag(atoms[i]) : '?').join('–');

  function drawMeasure() {
    const ul = $('measure-list'); ul.innerHTML = '';
    const P = picked.map(i => atoms[i]).filter(Boolean);
    P.forEach(a => viewer.addSphere({ center: { x: a.x, y: a.y, z: a.z }, radius: 0.45, color: '#f5b400', opacity: 0.55 }));
    for (let i = 0; i + 1 < P.length; i++) dashed(P[i], P[i + 1], '#f5b400');
    const need = NEED[mode];
    $('measure-hint').textContent = !cur ? 'Open a structure first.' :
      P.length ? `Picked ${P.map(tag).join(' – ')} (${P.length} of ${need}).` : `Click ${need} atoms in the structure.`;
    if (!cur) return;
    cur.measures.forEach((m, k) => {
      const A = m.idx.map(i => atoms[i]); const v = measureValue(m, i => atoms[i]); const c = MCOL[m.type];
      if (v !== null) {
        for (let i = 0; i + 1 < A.length; i++) dashed(A[i], A[i + 1], c);
        const at = m.type === 'distance' ? mid(A[0], A[1]) : m.type === 'angle' ? { x: (A[0].x + 2 * A[1].x + A[2].x) / 4, y: (A[0].y + 2 * A[1].y + A[2].y) / 4, z: (A[0].z + 2 * A[1].z + A[2].z) / 4 } : mid(A[1], A[2]);
        label(fmtValue(m, v), at, c);
      }
      const li = document.createElement('li');
      const dotEl = document.createElement('span'); dotEl.className = 'dot'; dotEl.style.background = c;
      const what = document.createElement('span'); what.className = 'what'; what.textContent = measureName(m); what.title = m.type + ' ' + measureName(m);
      const val = document.createElement('span'); val.className = 'val'; val.textContent = fmtValue(m, v);
      const x = document.createElement('button'); x.className = 'x'; x.textContent = '×'; x.title = 'Remove';
      x.onclick = () => { cur.measures.splice(k, 1); render(false); };
      li.append(dotEl, what, val, x); ul.append(li);
    });
  }

  /* ---------- per-frame data: time, energy, tables ---------- */
  /* frame -> row of the energy table. The step counter in a DCD header can carry an offset (equilibration, restarts),
     so when the header steps do not fall inside the table, the last frame is lined up with the last row instead. */
  function energyRow(e, i) {
    const E = e.energy; if (!E) return -1;
    const near = (arr, v) => { let lo = 0, hi = arr.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (arr[m] < v) lo = m; else hi = m; } return Math.abs(arr[lo] - v) <= Math.abs(arr[hi] - v) ? lo : hi; };
    if (e.traj && E.step) {
      const t = e.traj; const first = E.step[0], last = E.step[E.n - 1];
      const fits = MolTraj.frameStep(t, t.nframes - 1) <= last * 1.0001 && MolTraj.frameStep(t, 0) >= first - t.nsavc;
      return near(E.step, fits ? MolTraj.frameStep(t, i) : last - (t.nframes - 1 - i) * t.nsavc);
    }
    if (e.traj && E.t) return near(E.t, MolTraj.frameTime(e.traj, i));
    const n = nF(e); return n > 1 ? Math.round(i * (E.n - 1) / (n - 1)) : 0;
  }
  function frameTime(e, i) {                                        // ps, or null when unknown
    const r = energyRow(e, i);
    if (r >= 0 && e.energy.t) return e.energy.t[r];
    return e.traj ? MolTraj.frameTime(e.traj, i) : null;
  }
  function frameEnergy(e, i) {                                      // { v, unit } or null
    const r = energyRow(e, i); if (r >= 0) return { v: e.energy.e[r], unit: e.energy.eName };
    const f = e.traj ? null : e.frames[i]; return f && f.energy !== null && f.energy !== undefined ? { v: f.energy, unit: 'Eh' } : null;
  }
  const fmtTime = ps => (ps === null ? '' : ps >= 1000 ? (ps / 1000).toFixed(3) + ' ns' : ps.toFixed(1) + ' ps');

  async function measureTable() {
    if (!cur || !cur.measures.length) return '';
    const e = cur; const n = nF(e);
    const head = e.measures.map(m => m.type + ' ' + measureName(m) + (m.type === 'distance' ? ' (A)' : ' (deg)'));
    if (n < 2) return head.map((h, k) => h + '\t' + (measureValue(e.measures[k], i => atoms[i]) ?? NaN).toFixed(4)).join('\n');
    const rows = [['frame'].concat(e.traj ? ['time_ps'] : [], head).join('\t')];
    for (let i = 0; i < n; i++) {
      const X = await coordsOf(e, i); const P = X ? k => (k < X.x.length ? { x: X.x[k], y: X.y[k], z: X.z[k] } : undefined) : null;
      rows.push([i + 1].concat(e.traj ? [(frameTime(e, i) ?? NaN).toFixed(2)] : [], e.measures.map(m => (P ? (measureValue(m, P) ?? NaN).toFixed(4) : ''))).join('\t'));
      if (i % 50 === 0) $('measure-hint').textContent = `Reading frame ${i + 1} of ${n}…`;
    }
    return rows.join('\n');
  }

  async function scanContacts() {
    if (!cur || nF(cur) < 2) return;
    const e = cur; const n = nF(e); const cut = parseFloat($('cutoff').value) || 2.5;
    const inter = new Float32Array(n), intra = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const X = await coordsOf(e, i);
      if (X) { const C = contactsOf(e.topo, X, cut, null); let a = 0; C.forEach(c => { if (c.inter) a++; }); inter[i] = a; intra[i] = C.length - a; } else { inter[i] = NaN; intra[i] = NaN; }
      if (i % 25 === 0) { $('scan-status').textContent = `Scanning frame ${i + 1} of ${n}…`; await new Promise(r => setTimeout(r)); }
    }
    e.scan = { cut, inter, intra };
    const mean = a => a.reduce((s, v) => s + (isNaN(v) ? 0 : v), 0) / n; const frac = a => a.reduce((s, v) => s + (v > 0 ? 1 : 0), 0) / n * 100;
    $('scan-status').textContent = `H···N/O ≤ ${cut.toFixed(1)} Å over ${n} frames — between molecules: mean ${mean(inter).toFixed(2)}; intramolecular: mean ${mean(intra).toFixed(2)}, present in ${frac(intra).toFixed(0)} % of frames.`;
    if (cur === e) { updatePlotChoices('contacts'); drawPlot(); }
  }
  function scanTable() {
    if (!cur) return ''; const e = cur; const n = nF(e); const E = frameEnergy(e, 0);
    const head = ['frame', 'time_ps'].concat(E ? ['energy_' + E.unit] : [], e.scan ? ['contacts_between_molecules', 'contacts_intramolecular'] : []);
    const rows = [head.join('\t')];
    for (let i = 0; i < n; i++) { const t = frameTime(e, i); const en = frameEnergy(e, i); rows.push([i + 1, t === null ? '' : t.toFixed(2)].concat(E ? [en ? en.v : ''] : [], e.scan ? [e.scan.inter[i], e.scan.intra[i]] : []).join('\t')); }
    return rows.join('\n');
  }

  /* ---------- plot under the viewer ---------- */
  function plotSeries(e, what) {                                     // -> { x, ys: [{ y, color, name, step }], xName, yName, fx(i) = x of frame i }
    const n = nF(e);
    if (what === 'contacts' && e.scan) {
      const x = Float64Array.from({ length: n }, (_, i) => i + 1);
      return { x, ys: [{ y: e.scan.inter, color: INTER, name: 'between molecules' }, { y: e.scan.intra, color: INTRA, name: 'intramolecular' }], xName: 'frame', yName: 'H···N/O contacts', fx: i => i + 1, dots: true };
    }
    if (e.energy) {
      const E = e.energy; const useT = !!E.t; const x = useT ? E.t : Float64Array.from({ length: E.n }, (_, i) => i + 1);
      return { x, ys: [{ y: E.e, color: '#1f4e79', name: E.eName }], xName: useT ? (E.tName || 'time') : 'row', yName: E.eName, fx: i => { const r = energyRow(e, i); return useT ? E.t[r] : r + 1; }, smooth: E.n > 600 };
    }
    const en = e.frames.map(f => f.energy); if (e.traj || en.filter(v => v !== null && v !== undefined).length < 2) return null;
    const min = Math.min(...en.filter(v => v !== null && v !== undefined));
    return { x: Float64Array.from({ length: n }, (_, i) => i + 1), ys: [{ y: Float64Array.from(en, v => (v === null || v === undefined ? NaN : (v - min) * HARTREE)), color: '#1f4e79', name: 'E − E(min)' }], xName: 'frame', yName: 'E − E(min) (kJ/mol)', fx: i => i + 1, dots: true };
  }
  function updatePlotChoices(prefer) {
    const sel = $('plot-what'); const keep = prefer || sel.value; sel.innerHTML = '';
    if (cur) {
      if (cur.energy || (!cur.traj && cur.frames.filter(f => f.energy !== null && f.energy !== undefined).length > 1)) sel.add(new Option('Energy', 'energy'));
      if (cur.scan) sel.add(new Option('Contacts', 'contacts'));
    }
    if (Array.from(sel.options).some(o => o.value === keep)) sel.value = keep;
    const show = sel.options.length > 0; if ($('plot-wrap').hidden === show) { $('plot-wrap').hidden = !show; viewer.resize(); viewer.render(); }
  }
  function niceTicks(lo, hi, n) {
    const span = hi - lo || 1; const raw = span / n; const p = Math.pow(10, Math.floor(Math.log10(raw))); const f = raw / p;
    const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
  }
  function drawPlot() {
    plotMap = null; const wrap = $('plot-wrap'); if (wrap.hidden || !cur) return;
    const S = plotSeries(cur, $('plot-what').value); const cv = $('plot'); const dpr = window.devicePixelRatio || 1;
    const W = cv.clientWidth, H = cv.clientHeight; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    if (!S) return;
    const L = 68, R = 12, Tm = 8, B = 30; const pw = W - L - R, ph = H - Tm - B;
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < S.x.length; i++) { if (S.x[i] < x0) x0 = S.x[i]; if (S.x[i] > x1) x1 = S.x[i]; }
    S.ys.forEach(s => { for (let i = 0; i < s.y.length; i++) { const v = s.y[i]; if (isNaN(v)) continue; if (v < y0) y0 = v; if (v > y1) y1 = v; } });
    if (!(x1 > x0)) x1 = x0 + 1; if (!(y1 > y0)) { y0 -= 1; y1 += 1; } const pad = (y1 - y0) * 0.06; y0 -= pad; y1 += pad;
    const X = v => L + (v - x0) / (x1 - x0) * pw, Y = v => Tm + (1 - (v - y0) / (y1 - y0)) * ph;
    g.font = '11px "Segoe UI", Arial, sans-serif'; g.fillStyle = '#6b7280'; g.strokeStyle = '#e5e7eb'; g.lineWidth = 1;
    g.textAlign = 'right'; g.textBaseline = 'middle';
    niceTicks(y0, y1, 4).forEach(v => { g.beginPath(); g.moveTo(L, Y(v)); g.lineTo(L + pw, Y(v)); g.stroke(); g.fillText(Math.abs(v) >= 1e5 ? v.toExponential(3) : String(+v.toPrecision(7)), L - 6, Y(v)); });
    g.textAlign = 'center'; g.textBaseline = 'top';
    niceTicks(x0, x1, 6).forEach(v => { g.beginPath(); g.moveTo(X(v), Tm + ph); g.lineTo(X(v), Tm + ph + 4); g.stroke(); g.fillText(String(+v.toPrecision(6)), X(v), Tm + ph + 6); });
    g.fillText(S.xName, L + pw / 2, H - 13);
    g.strokeStyle = '#9ca3af'; g.strokeRect(L, Tm, pw, ph);
    S.ys.forEach(s => {
      const n = s.y.length; g.strokeStyle = s.color; g.fillStyle = s.color;
      if (S.smooth) {                                                 // long noisy trace: min–max envelope per pixel + block average
        const cols = Math.max(1, Math.floor(pw)); const lo = new Float64Array(cols).fill(Infinity), hi = new Float64Array(cols).fill(-Infinity), sm = new Float64Array(cols), ct = new Float64Array(cols);
        for (let i = 0; i < n; i++) { const c = Math.min(cols - 1, Math.floor((S.x[i] - x0) / (x1 - x0) * cols)); const v = s.y[i]; if (v < lo[c]) lo[c] = v; if (v > hi[c]) hi[c] = v; sm[c] += v; ct[c]++; }
        g.globalAlpha = 0.25; g.beginPath(); for (let c = 0; c < cols; c++) if (ct[c]) { g.moveTo(L + c + 0.5, Y(lo[c])); g.lineTo(L + c + 0.5, Y(hi[c])); } g.stroke(); g.globalAlpha = 1;
        const w = 12; g.lineWidth = 1.6; g.beginPath(); let started = false;
        for (let c = 0; c < cols; c++) { let a = 0, b = 0; for (let k = Math.max(0, c - w); k <= Math.min(cols - 1, c + w); k++) { a += sm[k]; b += ct[k]; } if (!b) continue; const yy = Y(a / b); if (started) g.lineTo(L + c + 0.5, yy); else { g.moveTo(L + c + 0.5, yy); started = true; } }
        g.stroke(); g.lineWidth = 1;
      } else {
        g.lineWidth = 1.4; g.beginPath(); let started = false;
        for (let i = 0; i < n; i++) { const v = s.y[i]; if (isNaN(v)) { started = false; continue; } if (started) g.lineTo(X(S.x[i]), Y(v)); else { g.moveTo(X(S.x[i]), Y(v)); started = true; } }
        g.stroke(); g.lineWidth = 1;
        if (S.dots && n <= 400) for (let i = 0; i < n; i++) { if (isNaN(s.y[i])) continue; g.beginPath(); g.arc(X(S.x[i]), Y(s.y[i]), 2.2, 0, 2 * Math.PI); g.fill(); }
      }
    });
    const n = nF(cur); const fx = S.fx(cur.i);                         // cursor at the frame on screen
    g.strokeStyle = '#f5b400'; g.lineWidth = 2; g.beginPath(); g.moveTo(X(fx), Tm); g.lineTo(X(fx), Tm + ph); g.stroke(); g.lineWidth = 1;
    plotMap = { n, fx: S.fx, x0, x1, L, pw };
    const bits = [`frame ${cur.i + 1} / ${n}`]; const t = frameTime(cur, cur.i); if (t !== null) bits.push('t = ' + fmtTime(t));
    const en = frameEnergy(cur, cur.i); if (en) bits.push(`E = ${en.unit === 'Eh' ? en.v.toFixed(6) : en.v.toFixed(2)} ${en.unit}`);
    if (cur.scan) bits.push(`contacts ${cur.scan.inter[cur.i]} between / ${cur.scan.intra[cur.i]} intra`);
    if (S.ys.length > 1) bits.push(S.ys.map(s => s.name).join(' · '));
    $('plot-read').textContent = bits.join('    ');
  }
  function plotClick(ev) {
    if (!plotMap || !cur) return; const r = $('plot').getBoundingClientRect(); const m = plotMap;
    const xv = m.x0 + (ev.clientX - r.left - m.L) / m.pw * (m.x1 - m.x0); let best = 0, bd = Infinity;
    for (let i = 0; i < m.n; i++) { const d = Math.abs(m.fx(i) - xv); if (d < bd) { bd = d; best = i; } }
    play(false); setFrame(best);
  }

  /* ---------- info ---------- */
  function formula() {
    const n = {}; const T = cur.topo; atoms.forEach((a, i) => { if (T.sol[i]) n[a.elem] = (n[a.elem] || 0) + 1; });
    const keys = Object.keys(n).sort((a, b) => (a === 'C' ? -2 : a === 'H' ? -1 : 0) - (b === 'C' ? -2 : b === 'H' ? -1 : 0) || a.localeCompare(b));
    return keys.map(k => k + (n[k] > 1 ? n[k] : '')).join(' ');
  }
  function updateInfo() {
    if (!cur) { $('info').textContent = 'No structure loaded.'; return; }
    const f = cur.traj ? cur.frames[0] : cur.frames[cur.i]; const T = cur.topo;
    const rows = [['File', cur.name], ['Format', cur.kind], ['Atoms', atoms.length + (vis.length !== atoms.length ? ` (${vis.length} shown)` : '')],
      [T.hasSolvent ? 'Solute' : 'Formula', T.solute.length <= 2000 ? formula() : '—']];
    if (f.comment) rows.push(['Comment', f.comment]);
    if (cur.traj) rows.push(['Trajectory', `${cur.traj.name} — ${cur.traj.nframes} frames, ${fmtTime(cur.traj.nsavc * cur.traj.dtps)} apart` + (cur.xyz && cur.xyz.cell ? `, box ${cur.xyz.cell.a.toFixed(1)} Å` : '')]);
    if (cur.energy) rows.push(['Energy table', `${cur.energy.name} — ${cur.energy.n} rows (${cur.energy.eName})`]);
    $('info').innerHTML = '';
    rows.forEach(([k, v]) => { const d = document.createElement('div'); const b = document.createElement('b'); b.textContent = k + ': '; d.append(b, String(v)); $('info').append(d); });
  }
  function updateTrajPanel() {
    const T = cur && cur.topo; const n = nF(cur); const show = !!cur && (n > 1 || (T && T.hasSolvent)); $('traj-panel').hidden = !show; if (!show) return;
    $('solvent-rows').hidden = !T.hasSolvent; $('scan').disabled = n < 2; $('copy-scan').disabled = n < 2;
    $('traj-info').textContent = (n > 1 ? `${n} frames. ` : '') + (T.hasSolvent ? `${T.solute.length} solute atoms + ${T.resAtoms.length - new Set(T.solute.map(i => T.res[i])).size} waters.` : '');
    if (!cur.scan) $('scan-status').textContent = n > 1 ? 'Scan contacts counts the H···N/O contacts (cut-off above) in every frame.' : '';
  }

  /* ---------- render ---------- */
  function computeVisible() {
    const T = cur.topo; const all = atoms.filter((a, i) => !T.virt[i]);
    const m = T.hasSolvent ? $('solvent').value : 'all';
    if (m === 'all') return all;
    if (m === 'solute') return atoms.filter((a, i) => T.sol[i]);
    const R = parseFloat($('near-r').value) || 3.5; const R2 = R * R; const S = T.solute.map(i => atoms[i]); const keep = new Uint8Array(T.resAtoms.length);
    atoms.forEach((a, i) => {
      if (T.sol[i] || T.virt[i] || keep[T.res[i]]) return;
      for (const s of S) { const dx = a.x - s.x, dy = a.y - s.y, dz = a.z - s.z; if (dx * dx + dy * dy + dz * dz <= R2) { keep[T.res[i]] = 1; break; } }
    });
    return atoms.filter((a, i) => !T.virt[i] && (T.sol[i] || keep[T.res[i]]));
  }

  function render(resetView) {
    viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels(); viewer.removeAllSurfaces();
    atoms = []; vis = [];
    if (cur) {
      const f = cur.traj ? cur.frames[0] : cur.frames[cur.i];
      const model = viewer.addModel(f.data, f.format, { keepH: true });     // 3Dmol drops PDB hydrogens unless asked
      atoms = model.selectedAtoms({});
      const T = cur.topo || buildTopo(cur, atoms);
      if (cur.traj && cur.xyz && cur.xyz.x.length === atoms.length) { const X = cur.xyz; atoms.forEach((a, i) => { a.x = X.x[i]; a.y = X.y[i]; a.z = X.z[i]; }); }
      else if (!cur.traj && f.charges && !cur.sets.some(s => s.fromFrame)) cur.sets.unshift({ name: cur.name, q: f.charges, note: 'Charge column of the structure file' + ((f.comment || '').includes('| q:') ? ':' + f.comment.split('| q:')[1] : ''), fromFrame: true });
      vis = computeVisible();
      updateChargePanel();
      applyStyle();
      const Q = activeQ();
      viewer.setClickable({}, true, a => { pick(a.index); });
      viewer.setHoverable({}, true,
        a => { $('hover').textContent = `${tag(a)}${a.resn ? ' ' + a.resn : ''}   x ${a.x.toFixed(3)}  y ${a.y.toFixed(3)}  z ${a.z.toFixed(3)}` + (Q && Q[a.index] !== null && Q[a.index] !== undefined ? `   q ${signed(Q[a.index], 4)}` : ''); },
        () => { $('hover').textContent = ''; });
      drawLabels(); drawContacts(); drawMeasure(); drawSurface();
      if (resetView) viewer.zoomTo(T.hasSolvent ? { index: T.solute } : {});
    } else { $('contact-list').textContent = ''; drawMeasure(); updateChargePanel(); drawSurface(); }
    viewer.render(); updateInfo(); updateFrameBar(); updateTrajPanel(); drawPlot();
  }

  function pick(i) {
    if (picked.includes(i)) picked = picked.filter(x => x !== i); else picked.push(i);
    if (cur && picked.length === NEED[mode]) { cur.measures.push({ type: mode, idx: picked.slice() }); picked = []; }     // enough atoms: keep the measurement
    render(false);
  }
  function setMode(m) {
    mode = m; picked = [];
    document.querySelectorAll('#mode button').forEach(b => b.classList.toggle('on', b.dataset.mode === m));
    render(false);
  }

  /* ---------- frames ---------- */
  function updateFrameBar() {
    const n = nF(cur);
    $('frame-bar').hidden = n < 2;
    if (n < 2) return;
    $('fr-slider').max = n - 1; $('fr-slider').value = cur.i;
    const t = frameTime(cur, cur.i); const e = frameEnergy(cur, cur.i);
    $('fr-label').textContent = `${cur.i + 1} / ${n}` + (t !== null ? `   ${fmtTime(t)}` : '') + (e && $('plot-wrap').hidden ? `   E = ${e.unit === 'Eh' ? e.v.toFixed(6) + ' Eh' : e.v.toFixed(1)}` : '');      // with the plot open the energy is in its read-out
  }
  async function setFrame(i, resetView) {
    if (!cur) return; const e = cur; const n = nF(e); const k = ((i % n) + n) % n; e.i = k;
    if (e.traj) {
      try { const X = await coordsOf(e, k); if (cur !== e || e.i !== k) return; e.xyz = X; }
      catch (err) { play(false); alert(err.message); return; }
    }
    render(!!resetView);
  }
  function play(on) {
    clearInterval(timer); timer = null; $('fr-play').innerHTML = '&#9654;';
    if (on && nF(cur) > 1) { timer = setInterval(() => setFrame(cur.i + 1), 250); $('fr-play').innerHTML = '&#10074;&#10074;'; }
  }

  /* ---------- files ---------- */
  function refreshList() {
    const ul = $('file-list'); ul.innerHTML = ''; $('empty-hint').hidden = files.length > 0;
    files.forEach(f => {
      const li = document.createElement('li'); if (f === cur) li.className = 'active';
      const nm = document.createElement('span'); nm.className = 'name'; nm.textContent = f.name; nm.title = f.name;
      const tg = document.createElement('span'); tg.className = 'tag'; tg.textContent = nF(f) > 1 ? nF(f) + ' frames' : '';
      const x = document.createElement('button'); x.className = 'x'; x.textContent = '×'; x.title = 'Close';
      x.onclick = ev => { ev.stopPropagation(); files.splice(files.indexOf(f), 1); if (cur === f) { cur = files[0] || null; picked = []; play(false); updatePlotChoices(); render(true); } refreshList(); };
      li.onclick = () => select(f);
      li.append(nm, tg, x); ul.append(li);
    });
  }
  function select(f) { play(false); cur = f; picked = []; refreshList(); updatePlotChoices(); render(true); }
  function addText(name, text) {
    const r = MolParsers.parse(name, text);
    const entry = { name, kind: r.kind, frames: r.frames, i: r.frames.length - 1, measures: [], sets: [], traj: null, xyz: null, energy: null, scan: null, topo: null };     // multi-frame files open on the last frame (final geometry)
    files.push(entry); select(entry);
  }
  const needCur = what => { if (!cur) throw new Error(`Open the structure first, then add the ${what}.`); return cur; };
  async function attachTraj(file) {
    const e = needCur('trajectory'); const t = await MolTraj.openDCD(file);
    if (!e.topo) throw new Error('The structure is not ready yet');
    if (t.natoms !== e.topo.n) throw new Error(`${file.name} has ${t.natoms} atoms, ${e.name} has ${e.topo.n}. Use the topology written with the trajectory (e.g. top.pdb).`);
    e.traj = t; e.scan = null; e.i = 0; refreshList(); updatePlotChoices(); await setFrame(0, true);
  }
  function attachEnergy(name, text) {
    const e = needCur('energy table'); const E = MolTraj.parseEnergyCSV(text); E.name = name; e.energy = E;
    updatePlotChoices('energy'); render(false);
  }
  function attachCharges(name, text) {
    const e = needCur('charges'); const f = e.traj ? e.frames[0] : e.frames[e.i];
    const r = MolCharges.parse(name, text, atoms.map(a => ({ elem: a.elem, atom: (a.atom || '').trim(), resn: (a.resn || '').trim(), resi: a.resi, chain: a.chain })), f.types || null);
    e.sets.push({ name, q: r.q, note: r.note, mp: r.mp || null }); $('q-set').value = ''; render(false);
  }
  const extOf = n => (n.split('.').pop() || '').toLowerCase();
  const isEnergy = (name, text) => extOf(name) === 'csv' && /potential|energy|epot|step|time/i.test(text.split('\n')[0] || '');
  async function addFiles(list) {
    const rank = f => (extOf(f.name) === 'dcd' ? 1 : extOf(f.name) === 'csv' || MolCharges.isChargeFile(f.name) ? 2 : 0);     // structures first, then what attaches to them
    for (const f of Array.from(list).sort((a, b) => rank(a) - rank(b))) {
      try {
        if (extOf(f.name) === 'dcd') { await attachTraj(f); continue; }
        const text = await f.text();
        if (isEnergy(f.name, text)) attachEnergy(f.name, text);
        else if (extOf(f.name) === 'csv' || MolCharges.isChargeFile(f.name)) attachCharges(f.name, text);
        else addText(f.name, text);
      } catch (err) { alert(f.name + '\n' + err.message); }
    }
  }

  /* ---------- events ---------- */
  function copyText(t, btn) {
    if (!t) return; const old = btn.textContent;
    const done = () => { btn.textContent = 'Copied'; setTimeout(() => { btn.textContent = old; }, 1200); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(t).then(done, () => window.prompt('Copy:', t));
    else window.prompt('Copy:', t);
  }
  $('file-input').onchange = e => { addFiles(e.target.files); e.target.value = ''; };
  document.querySelectorAll('[data-sample]').forEach(b => { b.onclick = () => fetch(b.dataset.sample).then(r => { if (!r.ok) throw new Error(r.status); return r.text(); }).then(t => addText(b.dataset.sample.split('/').pop(), t)).catch(e => alert('Could not load sample: ' + e.message)); });
  let depth = 0;
  window.addEventListener('dragenter', e => { e.preventDefault(); depth++; document.body.classList.add('dragging'); });
  window.addEventListener('dragleave', e => { e.preventDefault(); if (--depth <= 0) { depth = 0; document.body.classList.remove('dragging'); } });
  window.addEventListener('dragover', e => e.preventDefault());
  window.addEventListener('drop', e => { e.preventDefault(); depth = 0; document.body.classList.remove('dragging'); if (e.dataTransfer && e.dataTransfer.files.length) addFiles(e.dataTransfer.files); });

  let pend = null; const soon = () => { clearTimeout(pend); pend = setTimeout(() => render(false), 120); };
  ['style', 'labels', 'label-size', 'label-bg', 'hide-h', 'contacts', 'cutoff', 'solvent', 'near-r', 'q-set', 'q-show', 'q-scale', 'surf', 'surf-src', 'surf-grad', 'surf-unit', 'surf-op', 'surf-min', 'surf-max'].forEach(id => {
    const el = $(id); if (!el) return;                                   // an older cached page may not have the newer controls
    el.onchange = soon;
    if (el.tagName === 'INPUT' && (el.type === 'number' || el.type === 'range')) el.oninput = soon;     // react while typing / dragging, not only on Enter
  });
  $('dark').onchange = () => { viewer.setBackgroundColor(dark() ? 'black' : 'white'); render(false); };
  $('spin').onchange = e => viewer.spin(e.target.checked ? 'y' : false);
  $('reset-view').onclick = () => { viewer.zoomTo(cur && cur.topo && cur.topo.hasSolvent ? { index: cur.topo.solute } : {}); viewer.render(); };
  $('clear-measure').onclick = () => { picked = []; if (cur) cur.measures = []; render(false); };
  $('undo-pick').onclick = () => { picked.pop(); render(false); };
  document.querySelectorAll('#mode button').forEach(b => { b.onclick = () => setMode(b.dataset.mode); });
  $('copy-measure').onclick = async () => { const t = await measureTable(); render(false); copyText(t, $('copy-measure')); };
  $('copy-q').onclick = () => copyText(chargeTable(), $('copy-q'));
  $('clear-q').onclick = () => { if (cur) cur.sets = []; render(false); };
  $('scan').onclick = () => { scanContacts().catch(err => alert(err.message)); };
  $('copy-scan').onclick = () => copyText(scanTable(), $('copy-scan'));
  $('plot-what').onchange = () => drawPlot();
  $('plot').onclick = plotClick;
  $('fr-prev').onclick = () => { play(false); setFrame(cur.i - 1); };
  $('fr-next').onclick = () => { play(false); setFrame(cur.i + 1); };
  $('fr-play').onclick = () => play(!timer);
  $('fr-slider').oninput = e => { play(false); setFrame(parseInt(e.target.value, 10)); };
  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT' && e.target.type === 'number') return;
    if (e.key === 'Escape') { picked = []; render(false); }
    if (nF(cur) > 1 && e.key === 'ArrowRight') { play(false); setFrame(cur.i + 1); }
    if (nF(cur) > 1 && e.key === 'ArrowLeft') { play(false); setFrame(cur.i - 1); }
  });
  window.addEventListener('resize', () => { viewer.resize(); viewer.render(); drawPlot(); });

  function download(name, href) { const a = document.createElement('a'); a.href = href; a.download = name; document.body.append(a); a.click(); a.remove(); }
  $('save-png').onclick = savePNG;
  $('save-xyz').onclick = () => {                                     // what is on screen (solute + shown waters)
    if (!cur) return; const many = nF(cur) > 1;
    const text = MolParsers.xyzText(vis.map(a => ({ el: a.elem, x: a.x, y: a.y, z: a.z })), cur.name + (many ? ' frame ' + (cur.i + 1) : ''));
    download(cur.name.replace(/\.[^.]+$/, '') + (many ? '_' + (cur.i + 1) : '') + '.xyz', URL.createObjectURL(new Blob([text], { type: 'text/plain' })));
  };

  render(true);
  window.__molviewer = { addText, addFiles, attachCharges, attachEnergy, viewer, pick, setMode, setFrame, measureTable, scanContacts, scanTable, chargeTable,
    state: () => ({ files: files.map(f => f.name), atoms: atoms.length, shown: vis.length, picked, mode, measures: cur ? cur.measures : [], frame: cur ? cur.i : null, frames: nF(cur),
      sets: cur ? cur.sets.map(s => s.name) : [], traj: cur && cur.traj ? { n: cur.traj.nframes, natoms: cur.traj.natoms, dt: cur.traj.nsavc * cur.traj.dtps, istart: cur.traj.istart, cell: cur.xyz && cur.xyz.cell } : null,
      energy: cur && cur.energy ? cur.energy.n : 0, scan: cur && cur.scan ? { inter: Array.from(cur.scan.inter), intra: Array.from(cur.scan.intra) } : null }) };   // for debugging / automated checks
})();
