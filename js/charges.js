/* Web Mol Viewer — read atomic charges for a structure that is already open.
 *   MolCharges.parse(name, text, atoms, types) -> { q: [number | null per atom], note, mp }
 *     mp = atomic dipoles / quadrupoles with their local frames (see multipole.js), or null; OpenMM AmoebaMultipoleForce XML only
 *     atoms = [{ elem, atom (name), resn, resi, chain }] in structure order; types = Tinker atom types (or null)
 * Supported: OpenMM force-field XML (Multipole c0, NonbondedForce charge, residue-template charge),
 *            GROMACS .itp/.top ([ atoms ] section), Tinker .key/.prm (multipole / charge lines; needs a Tinker xyz),
 *            charge tables: "El x y z q" (Multiwfn .chg, DDEC6 xyz) or one value per line (last numeric column).
 */
'use strict';

const MolCharges = (() => {
  const isNum = s => s !== undefined && s !== '' && !isNaN(Number(s));
  const el2 = s => { const m = /^[A-Za-z]{1,2}/.exec(s || ''); return m ? m[0][0].toUpperCase() + (m[0][1] || '').toLowerCase() : ''; };

  function checkElements(q, atoms, elems, label) {              // warn when the charge file lists other elements than the structure
    let bad = 0;
    elems.forEach((e, i) => { if (e && atoms[i] && atoms[i].elem && e.toUpperCase() !== atoms[i].elem.toUpperCase()) bad++; });
    return bad ? `Warning: ${bad} atom(s) have a different element in ${label} — the atom order probably differs.` : '';
  }

  function fromOpenMM(text, atoms) {
    const doc = new DOMParser().parseFromString(text, 'application/xml');
    if (doc.querySelector('parsererror')) throw new Error('XML could not be parsed');
    const typeEl = {}, typeQ = {}, classOf = {}, classQ = {};
    doc.querySelectorAll('AtomTypes > Type').forEach(t => { typeEl[t.getAttribute('name')] = t.getAttribute('element') || ''; classOf[t.getAttribute('name')] = t.getAttribute('class'); });
    doc.querySelectorAll('NonbondedForce > Atom').forEach(a => {
      const c = a.getAttribute('charge'); if (c === null) return;
      if (a.getAttribute('type') !== null) typeQ[a.getAttribute('type')] = +c; else if (a.getAttribute('class') !== null) classQ[a.getAttribute('class')] = +c;
    });
    const typeMP = {};
    doc.querySelectorAll('AmoebaMultipoleForce > Multipole').forEach(m => {                                                    // multipole monopole wins
      const g = k => +(m.getAttribute(k) || 0); const ty = m.getAttribute('type'); typeQ[ty] = g('c0');
      typeMP[ty] = { k: ['kz', 'kx', 'ky'].map(a => (m.getAttribute(a) || '').trim()), d: [g('d1'), g('d2'), g('d3')],
        q: [g('q11'), g('q21'), g('q31'), g('q21'), g('q22'), g('q32'), g('q31'), g('q32'), g('q33')] };
    });
    const tmpl = [];
    doc.querySelectorAll('Residues > Residue').forEach(r => {
      const A = Array.from(r.querySelectorAll('Atom')).map(a => {
        const ty = a.getAttribute('type'); const own = a.getAttribute('charge');
        const q = ty in typeQ ? typeQ[ty] : classOf[ty] in classQ ? classQ[classOf[ty]] : own !== null ? +own : null;
        return { name: a.getAttribute('name'), el: typeEl[ty] || '', q, ty };
      });
      const at = {}; A.forEach((a, k) => { at[a.name] = k; }); const adj = A.map(() => []);
      r.querySelectorAll('Bond').forEach(b => {
        const n1 = b.getAttribute('atomName1'), n2 = b.getAttribute('atomName2');
        const u = n1 !== null ? at[n1] : +b.getAttribute('from'), v = n2 !== null ? at[n2] : +b.getAttribute('to');
        if (adj[u] && adj[v]) { adj[u].push(v); adj[v].push(u); }
      });
      tmpl.push({ name: r.getAttribute('name'), atoms: A, adj });
    });
    if (!tmpl.length) throw new Error('No <Residue> templates in this XML');
    const q = new Array(atoms.length).fill(null); const elems = new Array(atoms.length).fill('');
    const byName = {}; tmpl.forEach(t => { byName[t.name] = t; });
    /* template of a residue: same name, else a name the (truncated, PDB) residue name is a prefix of, else the same atom names */
    const pick = (resn, R) => {
      if (!resn) return null;
      if (byName[resn]) return byName[resn];
      const pre = tmpl.filter(t => t.name.startsWith(resn) || resn.startsWith(t.name));
      if (pre.length === 1 && pre[0].atoms.length === R.length) return pre[0];
      const names = R.map(a => a.atom).sort().join(' ');
      return tmpl.find(t => t.atoms.length === R.length && t.atoms.map(x => x.name).sort().join(' ') === names) || null;
    };
    /* dipoles and quadrupoles: unit change (nm -> Å) and the atoms that define each local frame.
       The force field names the frame atoms by type; in the residue it is the nearest bonded atom of that type. */
    const N = atoms.length; let mpAtoms = 0, mpLost = 0;
    const mp = Object.keys(typeMP).length ? { d: new Float64Array(3 * N), Q: new Float64Array(9 * N), ax: new Int32Array(3 * N).fill(-1), type: new Uint8Array(N).fill(MolMultipole.NOAXIS) } : null;
    const place = (t, map) => {                                 // map[template atom] = atom of the structure
      if (!mp) return;
      t.atoms.forEach((x, ti) => {
        const k = map[ti], M = typeMP[x.ty]; if (k === undefined || !M || !(M.d.some(v => v) || M.q.some(v => v))) return;
        const dist = new Array(t.atoms.length).fill(Infinity); dist[ti] = 0; const queue = [ti];
        for (const u of queue) for (const v of t.adj[u]) if (dist[v] === Infinity) { dist[v] = dist[u] + 1; queue.push(v); }
        const sgn = M.k.map(c => (c === '' || c === '0' ? 0 : c[0] === '-' ? -1 : 1)); const want = M.k.map(c => (c[0] === '-' ? c.slice(1) : c));
        const ax = [-1, -1, -1]; const taken = [ti]; let ok = true;
        sgn.forEach((sg, j) => {
          if (!sg) return; let best = -1;
          t.atoms.forEach((y, c) => { if (y.ty === want[j] && !taken.includes(c) && map[c] !== undefined && (best < 0 || dist[c] < dist[best])) best = c; });
          if (best < 0) { ok = false; return; }
          taken.push(best); ax[j] = map[best];
        });
        if (!ok) { mpLost++; return; }
        mp.type[k] = MolMultipole.axisType(sgn[0], sgn[1], sgn[2]); mp.ax.set(ax, 3 * k);
        mp.d.set(M.d.map(v => 10 * v), 3 * k); mp.Q.set(M.q.map(v => 100 * v), 9 * k); mpAtoms++;
      });
    };
    const hasRes = atoms.some(a => a.resn);
    let used = [];
    if (hasRes) {                                               // residue by residue: match atom names, else order inside the residue
      let i = 0;
      while (i < atoms.length) {
        let j = i; const a0 = atoms[i];
        while (j < atoms.length && atoms[j].resn === a0.resn && atoms[j].resi === a0.resi && atoms[j].chain === a0.chain) j++;
        const t = pick(a0.resn, atoms.slice(i, j));
        if (t) {
          if (!used.includes(t.name)) used.push(t.name);
          const names = {}; t.atoms.forEach(x => { names[x.name] = x; });
          const allNamed = atoms.slice(i, j).every(a => names[a.atom]);
          const map = {};
          for (let k = i; k < j; k++) { const x = allNamed ? names[atoms[k].atom] : t.atoms[k - i]; if (x) { q[k] = x.q; elems[k] = x.el; map[t.atoms.indexOf(x)] = k; } }
          place(t, map);
        }
        i = j;
      }
    }
    if (!used.length) {                                         // no residue matched (xyz, other residue names): template with the same number of atoms, else the first
      const t = tmpl.find(x => x.atoms.length === atoms.length) || tmpl[0];
      if (t.atoms.length !== atoms.length) throw new Error(`Template ${t.name} has ${t.atoms.length} atoms, the structure has ${atoms.length}`);
      const map = {}; t.atoms.forEach((x, k) => { q[k] = x.q; elems[k] = x.el; map[k] = k; }); place(t, map); used = [t.name];
    }
    const tag = doc.querySelector('AmoebaMultipoleForce > Multipole') ? ' (multipole monopoles' + (mpAtoms ? `; ${mpAtoms} atoms also carry dipoles / quadrupoles` : '') + ')' : '';
    const lost = mpLost ? `Warning: the local frame of ${mpLost} atom(s) could not be resolved; their dipoles / quadrupoles are left out` : '';
    return { q, mp: mp && mpAtoms ? mp : null, note: ['OpenMM XML' + tag + ', residue ' + used.join(', '), checkElements(q, atoms, elems, 'the XML'), lost].filter(Boolean).join('. ') };
  }

  function fromGromacs(text, atoms) {
    const q = [], elems = []; let sec = '', mols = 0;
    for (const raw of text.replace(/\r/g, '').split('\n')) {
      const l = raw.split(';')[0].trim(); if (!l) continue;
      const m = /^\[\s*(\w+)\s*\]/.exec(l);
      if (m) { sec = m[1].toLowerCase(); if (sec === 'atoms') mols++; continue; }
      if (sec !== 'atoms' || mols > 1) continue;                // first molecule only
      const p = l.split(/\s+/);
      if (p.length >= 7 && isNum(p[0]) && isNum(p[6])) { q.push(+p[6]); elems.push(el2(p[4])); }
    }
    if (!q.length) throw new Error('No [ atoms ] section with charges found');
    if (q.length !== atoms.length) throw new Error(`[ atoms ] lists ${q.length} atoms, the structure has ${atoms.length}`);
    return { q, note: ['GROMACS [ atoms ]', checkElements(q, atoms, elems, 'the topology')].filter(Boolean).join('. ') };
  }

  function fromTinker(text, atoms, types) {
    if (!types) throw new Error('A Tinker key file needs the atom types of a Tinker .xyz/.txyz structure');
    const tq = {};
    for (const raw of text.replace(/\r/g, '').split('\n')) {
      const p = raw.trim().split(/\s+/); const k = (p[0] || '').toLowerCase();
      if (k === 'multipole') { const nums = p.slice(1).filter(isNum); if (nums.length >= 2) tq[Math.abs(parseInt(nums[0], 10))] = +nums[nums.length - 1]; }
      else if (k === 'charge' && isNum(p[1]) && isNum(p[2])) tq[Math.abs(parseInt(p[1], 10))] = +p[2];
    }
    if (!Object.keys(tq).length) throw new Error('No multipole / charge lines found');
    const q = types.map(t => (t in tq ? tq[t] : null));
    const miss = q.filter(v => v === null).length;
    return { q, note: 'Tinker key (monopoles)' + (miss ? `. ${miss} atom type(s) not found in the key` : '') };
  }

  function fromTable(text, atoms) {
    const L = text.replace(/\r/g, '').split('\n').map(l => l.trim()).filter(l => l && !/^[#;!]/.test(l));
    const split = l => l.split(/[\s,]+/).filter(Boolean);
    const xyzq = L.map(split).filter(p => p.length >= 5 && !isNum(p[0]) && isNum(p[1]) && isNum(p[2]) && isNum(p[3]) && isNum(p[4]));
    if (xyzq.length >= atoms.length && atoms.length) {           // "El x y z q" (Multiwfn .chg, DDEC6 net_atomic_charges.xyz): first block of N lines
      const P = xyzq.slice(0, atoms.length);
      return { q: P.map(p => +p[4]), note: ['Table El x y z q', checkElements(null, atoms, P.map(p => el2(p[0])), 'the table')].filter(Boolean).join('. ') };
    }
    const last = L.map(split).filter(p => isNum(p[p.length - 1]));
    let P = last;
    if (P.length !== atoms.length) P = last.filter(p => p.length >= 2);          // drop a lone count line
    if (P.length !== atoms.length) throw new Error(`Found ${last.length} numeric lines, the structure has ${atoms.length} atoms`);
    const elems = P.map(p => (p.length >= 2 && !isNum(p[0]) ? el2(p[0]) : p.length >= 3 && !isNum(p[1]) ? el2(p[1]) : ''));
    return { q: P.map(p => +p[p.length - 1]), note: ['Table (last column)', checkElements(null, atoms, elems, 'the table')].filter(Boolean).join('. ') };
  }

  const EXT = ['xml', 'itp', 'top', 'key', 'prm', 'chg', 'txt', 'dat', 'charges', 'q'];
  const isChargeFile = name => EXT.includes((name.split('.').pop() || '').toLowerCase());

  function parse(name, text, atoms, types) {
    const ext = (name.split('.').pop() || '').toLowerCase();
    if (ext === 'xml' || /^\s*<(\?xml|ForceField)/.test(text)) return fromOpenMM(text, atoms);
    if (ext === 'itp' || ext === 'top') return fromGromacs(text, atoms);
    if (ext === 'key' || ext === 'prm') return fromTinker(text, atoms, types);
    return fromTable(text, atoms);
  }

  return { parse, isChargeFile, fromTable };
})();
