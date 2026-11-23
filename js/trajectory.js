/* Web Mol Viewer — trajectories and per-frame data.
 *   openDCD(file)        -> trajectory handle (header only; frames are read on demand, so large files are fine)
 *   readFrame(traj, i)   -> { x, y, z (Float32Array, Å), cell: { a, b, c, ortho } | null }
 *   parseEnergyCSV(text) -> { step, t, e, tName, eName }   (step / t may be null)
 * DCD = CHARMM / NAMD / OpenMM binary trajectory. The topology (atom names, bonds) comes from the structure file.
 */
'use strict';

const MolTraj = (() => {
  const AKMA_PS = 0.04888821;                                  // CHARMM time unit in ps

  async function openDCD(file) {
    const buf = await file.slice(0, Math.min(file.size, 65536)).arrayBuffer(); const dv = new DataView(buf);
    if (buf.byteLength < 104) throw new Error('File too short to be a DCD trajectory');
    let le = true;
    if (dv.getInt32(0, true) !== 84) { if (dv.getInt32(0, false) === 84) le = false; else throw new Error('Not a DCD trajectory (bad header)'); }
    const magic = String.fromCharCode(dv.getUint8(4), dv.getUint8(5), dv.getUint8(6), dv.getUint8(7));
    if (magic !== 'CORD') throw new Error('Not a coordinate DCD (found "' + magic + '")');
    const ic = k => dv.getInt32(8 + 4 * k, le);
    const charmm = ic(19) !== 0;
    if (ic(8) !== 0) throw new Error('DCD files with fixed atoms are not supported');
    const hasCell = charmm && ic(10) !== 0, has4d = charmm && ic(11) !== 0;
    const delta = charmm ? dv.getFloat32(8 + 36, le) : dv.getFloat64(8 + 36, le);
    let off = 92;
    const tlen = dv.getInt32(off, le); off += 8 + tlen;        // title block
    if (off + 12 > buf.byteLength) throw new Error('DCD title block is too long');
    const natoms = dv.getInt32(off + 4, le); off += 12;
    if (!(natoms > 0)) throw new Error('DCD: bad atom count');
    const frameSize = (hasCell ? 56 : 0) + (has4d ? 4 : 3) * (4 * natoms + 8);
    const nframes = Math.floor((file.size - off) / frameSize);
    if (nframes < 1) throw new Error('DCD contains no complete frame');
    return { file, name: file.name, le, natoms, nframes, hasCell, frameSize, headerSize: off,
      istart: ic(1), nsavc: Math.max(1, ic(2)), dtps: delta * AKMA_PS };
  }

  const frameStep = (t, i) => t.istart + i * t.nsavc;
  const frameTime = (t, i) => frameStep(t, i) * t.dtps;        // ps

  async function readFrame(t, i) {
    const start = t.headerSize + i * t.frameSize;
    const buf = await t.file.slice(start, start + t.frameSize).arrayBuffer();
    if (buf.byteLength < t.frameSize) throw new Error('Frame ' + (i + 1) + ' is incomplete');
    const dv = new DataView(buf); let o = 0, cell = null;
    if (t.hasCell) {
      const v = []; for (let k = 0; k < 6; k++) v.push(dv.getFloat64(4 + 8 * k, t.le));
      const ang = [v[1], v[3], v[4]];                           // degrees, or cosines in old NAMD files
      const ortho = ang.every(g => Math.abs(g - 90) < 0.01 || Math.abs(g) < 1e-4);
      cell = { a: v[0], b: v[2], c: v[5], ortho }; o = 56;
    }
    const N = t.natoms;
    const block = () => {
      let arr;
      if (t.le) arr = new Float32Array(buf.slice(o + 4, o + 4 + 4 * N));
      else { arr = new Float32Array(N); for (let k = 0; k < N; k++) arr[k] = dv.getFloat32(o + 4 + 4 * k, false); }
      o += 4 * N + 8; return arr;
    };
    const x = block(), y = block(), z = block();
    return { x, y, z, cell };
  }

  /* CSV / whitespace table with a header row; picks a step column, a time column and an energy column. */
  function parseEnergyCSV(text) {
    const L = text.replace(/\r/g, '').split('\n').filter(l => l.trim() && !/^\s*[#@]/.test(l));
    if (L.length < 2) throw new Error('Energy table is empty');
    const sep = L[0].includes(',') ? ',' : L[0].includes('\t') ? '\t' : /\s+/;
    const split = l => l.trim().split(sep).map(s => s.trim().replace(/^"|"$/g, ''));
    const head = split(L[0]); const hasHead = head.some(h => isNaN(Number(h)));
    const names = hasHead ? head : head.map((_, k) => 'col' + (k + 1));
    const rows = L.slice(hasHead ? 1 : 0);
    const find = re => names.findIndex(n => re.test(n));
    const iStep = find(/^#?\s*"?step/i), iTime = find(/time/i);
    let iE = find(/potential/i); if (iE < 0) iE = find(/energy|epot|^e$/i);
    if (iE < 0) iE = names.length > 1 ? (iStep === 0 || iTime === 0 ? (iTime === 1 || iStep === 1 ? Math.min(2, names.length - 1) : 1) : names.length - 1) : 0;
    const n = rows.length; const e = new Float64Array(n);
    const step = iStep >= 0 ? new Float64Array(n) : null, t = iTime >= 0 ? new Float64Array(n) : null;
    let m = 0;
    for (const l of rows) {
      const p = split(l); const v = Number(p[iE]); if (isNaN(v)) continue;
      e[m] = v; if (step) step[m] = Number(p[iStep]); if (t) t[m] = Number(p[iTime]); m++;
    }
    if (!m) throw new Error('No numeric energy values found (column "' + names[iE] + '")');
    return { n: m, e: e.subarray(0, m), step: step ? step.subarray(0, m) : null, t: t ? t.subarray(0, m) : null,
      eName: names[iE], tName: iTime >= 0 ? names[iTime] : null };
  }

  return { openDCD, readFrame, frameStep, frameTime, parseEnergyCSV };
})();
