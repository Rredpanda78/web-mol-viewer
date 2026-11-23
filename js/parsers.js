/* Parsers that turn various text formats into "frames" the viewer can show.
   A frame is { data: <text>, format: <3Dmol format name>, comment: <string>, energy: <number|null> }.
   Formats 3Dmol.js already understands (pdb, sdf, mol2, gro, cif, cube) are passed through as a single frame. */
'use strict';

const MolParsers = (() => {
  const SYM = [null, 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
    'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'];
  const isNum = s => s !== undefined && s !== '' && !isNaN(Number(s));
  const cleanEl = s => { const m = /^[A-Za-z]{1,2}/.exec(s || ''); if (!m) return 'X'; const e = m[0]; return e[0].toUpperCase() + (e[1] || '').toLowerCase(); };
  const xyzText = (atoms, comment) => atoms.length + '\n' + (comment || '').replace(/\n/g, ' ') + '\n' +
    atoms.map(a => a.el.padEnd(2) + ' ' + a.x.toFixed(6).padStart(12) + ' ' + a.y.toFixed(6).padStart(12) + ' ' + a.z.toFixed(6).padStart(12)).join('\n') + '\n';

  /* Standard (possibly multi-frame) XYZ. */
  function xyzFrames(text) {
    const L = text.replace(/\r/g, '').split('\n'); const frames = []; let i = 0;
    while (i < L.length) {
      const n = parseInt(L[i].trim(), 10);
      if (!(n > 0)) { i++; continue; }
      const atoms = [];
      for (let k = i + 2; k < i + 2 + n && k < L.length; k++) {
        const p = L[k].trim().split(/\s+/);
        if (p.length < 4 || !isNum(p[1]) || !isNum(p[2]) || !isNum(p[3])) break;
        atoms.push({ el: isNum(p[0]) ? (SYM[parseInt(p[0], 10)] || 'X') : cleanEl(p[0]), x: +p[1], y: +p[2], z: +p[3] });
      }
      if (atoms.length !== n) throw new Error('XYZ frame ' + (frames.length + 1) + ': expected ' + n + ' atoms, found ' + atoms.length);
      const comment = (L[i + 1] || '').trim();
      frames.push({ data: xyzText(atoms, comment), format: 'xyz', comment, energy: null });
      i += 2 + n;
    }
    if (!frames.length) throw new Error('No XYZ frames found');
    return frames;
  }

  /* Tinker XYZ / ARC: "N title" then "index element x y z type bonded...". */
  function looksLikeTinker(text) {
    const L = text.replace(/\r/g, '').split('\n').filter(l => l.trim());
    if (L.length < 2) return false;
    const h = L[0].trim().split(/\s+/), p = L[1].trim().split(/\s+/);
    return isNum(h[0]) && p.length >= 6 && p[0] === '1' && !isNum(p[1]) && isNum(p[2]) && isNum(p[3]) && isNum(p[4]) && isNum(p[5]);
  }
  function tinkerFrames(text) {
    const L = text.replace(/\r/g, '').split('\n').filter(l => l.trim()); const frames = []; let i = 0;
    while (i < L.length) {
      const h = L[i].trim().split(/\s+/); const n = parseInt(h[0], 10);
      if (!(n > 0)) break;
      let k = i + 1;
      if (k < L.length && L[k].trim().split(/\s+/).length === 6 && L[k].trim().split(/\s+/).every(isNum)) k++;   // periodic box line
      const atoms = [];
      for (let j = k; j < k + n && j < L.length; j++) {
        const p = L[j].trim().split(/\s+/);
        atoms.push({ el: cleanEl(p[1]), x: +p[2], y: +p[3], z: +p[4] });
      }
      if (atoms.length !== n) throw new Error('Tinker XYZ: expected ' + n + ' atoms');
      const comment = h.slice(1).join(' ');
      frames.push({ data: xyzText(atoms, comment), format: 'xyz', comment, energy: null });
      i = k + n;
    }
    if (!frames.length) throw new Error('No Tinker frames found');
    return frames;
  }

  /* Gaussian output: every geometry of the run becomes a frame (optimisation steps, scan points). */
  function gaussianLogFrames(text) {
    const L = text.replace(/\r/g, '').split('\n');
    const blocks = { 'Input orientation:': [], 'Standard orientation:': [], 'Z-Matrix orientation:': [] };
    const energies = [];
    for (let i = 0; i < L.length; i++) {
      const t = L[i].trim();
      if (blocks[t] !== undefined) {
        const atoms = []; let k = i + 5;
        for (; k < L.length && !/^\s*-{5,}/.test(L[k]); k++) {
          const p = L[k].trim().split(/\s+/);
          if (p.length >= 6) atoms.push({ el: SYM[parseInt(p[1], 10)] || 'X', x: +p[3], y: +p[4], z: +p[5] });
        }
        if (atoms.length) blocks[t].push({ atoms, line: i });
        i = k;
      } else if (t.startsWith('SCF Done:')) {
        const m = /=\s*(-?\d+\.\d+)/.exec(t); if (m) energies.push({ e: +m[1], line: i });
      }
    }
    const use = blocks['Input orientation:'].length ? blocks['Input orientation:'] :
      blocks['Standard orientation:'].length ? blocks['Standard orientation:'] : blocks['Z-Matrix orientation:'];
    if (!use.length) throw new Error('No geometry found in this Gaussian output');
    const ok = /Normal termination/.test(L.slice(-12).join('\n'));
    return use.map((b, idx) => {
      const next = idx + 1 < use.length ? use[idx + 1].line : Infinity;
      const en = energies.find(x => x.line > b.line && x.line < next);
      const comment = 'step ' + (idx + 1) + '/' + use.length + (en ? '  E = ' + en.e.toFixed(8) + ' Eh' : '') + (idx === use.length - 1 && ok ? '  (final)' : '');
      return { data: xyzText(b.atoms, comment), format: 'xyz', comment, energy: en ? en.e : null };
    });
  }

  /* Gaussian input (.gjf / .com) with Cartesian coordinates. */
  function gaussianInputFrames(text) {
    const L = text.replace(/\r/g, '').split('\n'); let i = 0;
    while (i < L.length && !L[i].trim().startsWith('#')) i++;
    if (i === L.length) throw new Error('No route line (#) found');
    while (i < L.length && L[i].trim()) i++;          // end of route
    while (i < L.length && !L[i].trim()) i++;         // blank
    const title = L[i] || '';
    while (i < L.length && L[i].trim()) i++;          // end of title
    while (i < L.length && !L[i].trim()) i++;         // blank
    i++;                                              // charge / multiplicity
    const atoms = [];
    for (; i < L.length && L[i].trim(); i++) {
      const p = L[i].trim().split(/[\s,]+/); const c = p.slice(1).filter(isNum).slice(-3);
      if (c.length === 3) atoms.push({ el: isNum(p[0]) ? (SYM[parseInt(p[0], 10)] || 'X') : cleanEl(p[0].split('(')[0].split('-')[0]), x: +c[0], y: +c[1], z: +c[2] });
    }
    if (!atoms.length) throw new Error('No Cartesian coordinates found (Z-matrix input is not supported yet)');
    return [{ data: xyzText(atoms, title.trim()), format: 'xyz', comment: title.trim(), energy: null }];
  }

  const PASS = { pdb: 'pdb', ent: 'pdb', sdf: 'sdf', mol: 'sdf', mol2: 'mol2', cif: 'cif', gro: 'gro', cube: 'cube', cub: 'cube' };

  /* Entry point: file name + text -> { kind, frames }. */
  function parse(name, text) {
    const ext = (name.split('.').pop() || '').toLowerCase();
    if (ext === 'log' || ext === 'out') return { kind: 'Gaussian output', frames: gaussianLogFrames(text) };
    if (ext === 'gjf' || ext === 'com') return { kind: 'Gaussian input', frames: gaussianInputFrames(text) };
    if (ext === 'txyz' || ext === 'arc' || (ext === 'xyz' && looksLikeTinker(text))) return { kind: 'Tinker XYZ', frames: tinkerFrames(text) };
    if (ext === 'xyz') return { kind: 'XYZ', frames: xyzFrames(text) };
    if (PASS[ext]) return { kind: ext.toUpperCase(), frames: [{ data: text, format: PASS[ext], comment: '', energy: null }] };
    if (/^\s*\d+\s*$/.test(text.split('\n')[0] || '')) return { kind: 'XYZ', frames: xyzFrames(text) };     // unknown extension that looks like XYZ
    throw new Error('Unsupported file type: .' + ext);
  }

  return { parse, xyzText };
})();
