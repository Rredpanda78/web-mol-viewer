/* Web Mol Viewer — UI and viewer logic (3Dmol.js does the drawing). */
'use strict';

(() => {
  const $ = id => document.getElementById(id);
  const viewer = $3Dmol.createViewer($('viewer'), { backgroundColor: 'white' });
  const files = [];            // { name, kind, frames, i }
  let cur = null;              // active entry
  let atoms = [];              // atoms of the frame on screen
  let picked = [];             // atom indices picked so far for the measurement being made
  let mode = 'distance';       // measuring tool: distance (2 atoms), angle (3), dihedral (4)
  let timer = null;
  const NEED = { distance: 2, angle: 3, dihedral: 4 };
  const MCOL = { distance: '#1f4e79', angle: '#2e7d32', dihedral: '#7030a0' };

  const STYLE = {
    ballstick: { stick: { radius: 0.14, colorscheme: 'Jmol' }, sphere: { scale: 0.27, colorscheme: 'Jmol' } },
    stick: { stick: { radius: 0.18, colorscheme: 'Jmol' } },
    line: { line: { linewidth: 2, colorscheme: 'Jmol' } },
    sphere: { sphere: { colorscheme: 'Jmol' } },
  };

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

  /* ---------- drawing ---------- */
  function applyStyle() {
    const s = $('style').value;
    if (s === 'cartoon') {
      viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
      viewer.setStyle({ hetflag: true }, STYLE.ballstick);
    } else viewer.setStyle({}, STYLE[s]);
    if ($('hide-h').checked) viewer.setStyle({ elem: 'H' }, {});
  }

  function label(text, pos, color, size) {
    viewer.addLabel(text, { position: pos, fontSize: size || 13, fontColor: color || (dark() ? 'white' : 'black'), showBackground: true,
      backgroundColor: dark() ? 'black' : 'white', backgroundOpacity: 0.65, inFront: true, alignment: 'center' });
  }
  function dashed(a, b, color) {
    viewer.addCylinder({ start: { x: a.x, y: a.y, z: a.z }, end: { x: b.x, y: b.y, z: b.z }, radius: 0.045, dashed: true, color, fromCap: 1, toCap: 1 });
  }

  function drawLabels() {
    const m = $('labels').value; if (m === 'none') return;
    const hideH = $('hide-h').checked;
    if (atoms.length > 400) { $('hover').textContent = 'Labels are switched off above 400 atoms.'; return; }
    atoms.forEach(a => {
      if (hideH && a.elem === 'H') return;
      if (m === 'elem' && (a.elem === 'C' || a.elem === 'H')) return;          // element labels: hetero atoms only
      label(m === 'elem' ? a.elem : m === 'index' ? String(a.index + 1) : tag(a), a, null, 11);
    });
  }

  /* polar H···N/O contacts: H bonded to N/O, acceptor N/O at least 3 bonds away (or in another molecule) */
  function findContacts(cut) {
    const out = [];
    atoms.forEach(h => {
      if (h.elem !== 'H' || !h.bonds || !h.bonds.length) return;
      const d = atoms[h.bonds[0]]; if (!d || (d.elem !== 'N' && d.elem !== 'O')) return;
      const near = new Set([d.index]); let layer = [d.index];
      for (let k = 0; k < 2; k++) {                          // atoms within 2 bonds of the donor are skipped
        const nxt = [];
        layer.forEach(i => (atoms[i].bonds || []).forEach(j => { if (!near.has(j)) { near.add(j); nxt.push(j); } }));
        layer = nxt;
      }
      atoms.forEach(y => {
        if ((y.elem !== 'N' && y.elem !== 'O') || near.has(y.index)) return;
        const r = dist(h, y);
        if (r <= cut) out.push({ h, d, y, r, ang: angle(d, h, y) });
      });
    });
    return out.sort((p, q) => p.r - q.r);
  }
  function drawContacts() {
    const box = $('contact-list');
    if (!$('contacts').checked || !atoms.length) { box.textContent = ''; return; }
    const cut = parseFloat($('cutoff').value) || 2.5; const C = findContacts(cut);
    C.forEach(c => { dashed(c.h, c.y, '#d62728'); label(c.r.toFixed(2), mid(c.h, c.y), '#d62728'); });
    box.textContent = C.length ? C.map(c => `${tag(c.d)}–${tag(c.h)}···${tag(c.y)}  ${c.r.toFixed(2)} Å  ${c.ang.toFixed(0)}°`).join('\n') : `No contact below ${cut.toFixed(1)} Å.`;
  }

  /* value of one measurement for a list of atoms ({x, y, z}); null if an atom is missing */
  function measureValue(m, A) {
    const P = m.idx.map(i => A[i]); if (P.some(p => !p)) return null;
    return m.type === 'distance' ? dist(P[0], P[1]) : m.type === 'angle' ? angle(P[0], P[1], P[2]) : dihedral(P[0], P[1], P[2], P[3]);
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
      const A = m.idx.map(i => atoms[i]); const v = measureValue(m, atoms); const c = MCOL[m.type];
      if (v !== null) {
        for (let i = 0; i + 1 < A.length; i++) dashed(A[i], A[i + 1], c);
        const at = m.type === 'distance' ? mid(A[0], A[1]) : m.type === 'angle' ? { x: (A[0].x + 2 * A[1].x + A[2].x) / 4, y: (A[0].y + 2 * A[1].y + A[2].y) / 4, z: (A[0].z + 2 * A[1].z + A[2].z) / 4 } : mid(A[1], A[2]);
        label(fmtValue(m, v), at, c);
      }
      const li = document.createElement('li');
      const dot = document.createElement('span'); dot.className = 'dot'; dot.style.background = c;
      const what = document.createElement('span'); what.className = 'what'; what.textContent = measureName(m); what.title = m.type + ' ' + measureName(m);
      const val = document.createElement('span'); val.className = 'val'; val.textContent = fmtValue(m, v);
      const x = document.createElement('button'); x.className = 'x'; x.textContent = '×'; x.title = 'Remove';
      x.onclick = () => { cur.measures.splice(k, 1); render(false); };
      li.append(dot, what, val, x); ul.append(li);
    });
  }

  /* coordinates of every frame (only for frames this page built itself, i.e. XYZ text) */
  function frameAtoms(f) {
    if (f.format !== 'xyz') return null;
    const L = f.data.split('\n'); const n = parseInt(L[0], 10);
    return L.slice(2, 2 + n).map(l => { const p = l.trim().split(/\s+/); return { x: +p[1], y: +p[2], z: +p[3] }; });
  }
  function measureTable() {
    if (!cur || !cur.measures.length) return '';
    const head = cur.measures.map(m => m.type + ' ' + measureName(m) + (m.type === 'distance' ? ' (A)' : ' (deg)'));
    if (cur.frames.length < 2) return head.map((h, k) => h + '\t' + (measureValue(cur.measures[k], atoms) ?? NaN).toFixed(4)).join('\n');
    const rows = [['frame'].concat(head).join('\t')];
    cur.frames.forEach((f, i) => {
      const A = i === cur.i ? atoms : frameAtoms(f);
      rows.push([i + 1].concat(cur.measures.map(m => A ? (measureValue(m, A) ?? NaN).toFixed(4) : '')).join('\t'));
    });
    return rows.join('\n');
  }

  function formula() {
    const n = {}; atoms.forEach(a => { n[a.elem] = (n[a.elem] || 0) + 1; });
    const keys = Object.keys(n).sort((a, b) => (a === 'C' ? -2 : a === 'H' ? -1 : 0) - (b === 'C' ? -2 : b === 'H' ? -1 : 0) || a.localeCompare(b));
    return keys.map(k => k + (n[k] > 1 ? n[k] : '')).join(' ');
  }
  function updateInfo() {
    if (!cur) { $('info').textContent = 'No structure loaded.'; return; }
    const f = cur.frames[cur.i];
    $('info').innerHTML = '';
    [['File', cur.name], ['Format', cur.kind], ['Atoms', atoms.length], ['Formula', atoms.length <= 2000 ? formula() : '—'],
      ['Comment', f.comment || '—']].forEach(([k, v]) => { const d = document.createElement('div'); const b = document.createElement('b'); b.textContent = k + ': '; d.append(b, String(v)); $('info').append(d); });
  }

  function render(resetView) {
    viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels();
    atoms = [];
    if (cur) {
      const f = cur.frames[cur.i];
      const model = viewer.addModel(f.data, f.format);
      atoms = model.selectedAtoms({});
      applyStyle();
      viewer.setClickable({}, true, a => { pick(a.index); });
      viewer.setHoverable({}, true,
        a => { $('hover').textContent = `${tag(a)}   x ${a.x.toFixed(3)}  y ${a.y.toFixed(3)}  z ${a.z.toFixed(3)}`; },
        () => { $('hover').textContent = ''; });
      drawLabels(); drawContacts(); drawMeasure();
      if (resetView) viewer.zoomTo();
    } else { $('contact-list').textContent = ''; drawMeasure(); }
    viewer.render(); updateInfo(); updateFrameBar();
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
    const n = cur ? cur.frames.length : 0;
    $('frame-bar').hidden = n < 2;
    if (n < 2) return;
    $('fr-slider').max = n - 1; $('fr-slider').value = cur.i;
    const e = cur.frames[cur.i].energy;
    $('fr-label').textContent = `${cur.i + 1} / ${n}` + (e !== null ? `   E = ${e.toFixed(6)} Eh` : '');
  }
  function setFrame(i) { if (!cur) return; const n = cur.frames.length; cur.i = ((i % n) + n) % n; render(false); }
  function play(on) {
    clearInterval(timer); timer = null; $('fr-play').innerHTML = '&#9654;';
    if (on && cur && cur.frames.length > 1) { timer = setInterval(() => setFrame(cur.i + 1), 250); $('fr-play').innerHTML = '&#10074;&#10074;'; }
  }

  /* ---------- files ---------- */
  function refreshList() {
    const ul = $('file-list'); ul.innerHTML = ''; $('empty-hint').hidden = files.length > 0;
    files.forEach(f => {
      const li = document.createElement('li'); if (f === cur) li.className = 'active';
      const nm = document.createElement('span'); nm.className = 'name'; nm.textContent = f.name; nm.title = f.name;
      const tg = document.createElement('span'); tg.className = 'tag'; tg.textContent = f.frames.length > 1 ? f.frames.length + ' frames' : '';
      const x = document.createElement('button'); x.className = 'x'; x.textContent = '×'; x.title = 'Close';
      x.onclick = ev => { ev.stopPropagation(); files.splice(files.indexOf(f), 1); if (cur === f) { cur = files[0] || null; picked = []; play(false); render(true); } refreshList(); };
      li.onclick = () => select(f);
      li.append(nm, tg, x); ul.append(li);
    });
  }
  function select(f) { play(false); cur = f; picked = []; refreshList(); render(true); }
  function addText(name, text) {
    try {
      const r = MolParsers.parse(name, text);
      const entry = { name, kind: r.kind, frames: r.frames, i: r.frames.length - 1, measures: [] };     // multi-frame files open on the last frame (final geometry)
      files.push(entry); select(entry);
    } catch (err) { alert(name + '\n' + err.message); }
  }
  function addFiles(list) { Array.from(list).forEach(f => { const rd = new FileReader(); rd.onload = () => addText(f.name, rd.result); rd.readAsText(f); }); }

  /* ---------- events ---------- */
  $('file-input').onchange = e => { addFiles(e.target.files); e.target.value = ''; };
  document.querySelectorAll('[data-sample]').forEach(b => { b.onclick = () => fetch(b.dataset.sample).then(r => { if (!r.ok) throw new Error(r.status); return r.text(); }).then(t => addText(b.dataset.sample.split('/').pop(), t)).catch(e => alert('Could not load sample: ' + e.message)); });
  let depth = 0;
  window.addEventListener('dragenter', e => { e.preventDefault(); depth++; document.body.classList.add('dragging'); });
  window.addEventListener('dragleave', e => { e.preventDefault(); if (--depth <= 0) { depth = 0; document.body.classList.remove('dragging'); } });
  window.addEventListener('dragover', e => e.preventDefault());
  window.addEventListener('drop', e => { e.preventDefault(); depth = 0; document.body.classList.remove('dragging'); if (e.dataTransfer && e.dataTransfer.files.length) addFiles(e.dataTransfer.files); });

  ['style', 'labels', 'hide-h', 'contacts', 'cutoff'].forEach(id => { $(id).onchange = () => render(false); });
  $('dark').onchange = () => { viewer.setBackgroundColor(dark() ? 'black' : 'white'); render(false); };
  $('spin').onchange = e => viewer.spin(e.target.checked ? 'y' : false);
  $('reset-view').onclick = () => { viewer.zoomTo(); viewer.render(); };
  $('clear-measure').onclick = () => { picked = []; if (cur) cur.measures = []; render(false); };
  $('undo-pick').onclick = () => { picked.pop(); render(false); };
  document.querySelectorAll('#mode button').forEach(b => { b.onclick = () => setMode(b.dataset.mode); });
  $('copy-measure').onclick = () => {
    const t = measureTable(); if (!t) return;
    const done = () => { $('copy-measure').textContent = 'Copied'; setTimeout(() => { $('copy-measure').textContent = 'Copy'; }, 1200); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(t).then(done, () => window.prompt('Copy the measurements:', t));
    else window.prompt('Copy the measurements:', t);
  };
  $('fr-prev').onclick = () => { play(false); setFrame(cur.i - 1); };
  $('fr-next').onclick = () => { play(false); setFrame(cur.i + 1); };
  $('fr-play').onclick = () => play(!timer);
  $('fr-slider').oninput = e => { play(false); setFrame(parseInt(e.target.value, 10)); };
  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT' && e.target.type === 'number') return;
    if (e.key === 'Escape') { picked = []; render(false); }
    if (cur && cur.frames.length > 1 && e.key === 'ArrowRight') { play(false); setFrame(cur.i + 1); }
    if (cur && cur.frames.length > 1 && e.key === 'ArrowLeft') { play(false); setFrame(cur.i - 1); }
  });
  window.addEventListener('resize', () => { viewer.resize(); viewer.render(); });

  function download(name, href) { const a = document.createElement('a'); a.href = href; a.download = name; document.body.append(a); a.click(); a.remove(); }
  $('save-png').onclick = () => { if (cur) download(cur.name.replace(/\.[^.]+$/, '') + '.png', viewer.pngURI()); };
  $('save-xyz').onclick = () => {
    if (!cur) return;
    const text = MolParsers.xyzText(atoms.map(a => ({ el: a.elem, x: a.x, y: a.y, z: a.z })), cur.name + (cur.frames.length > 1 ? ' frame ' + (cur.i + 1) : ''));
    download(cur.name.replace(/\.[^.]+$/, '') + (cur.frames.length > 1 ? '_' + (cur.i + 1) : '') + '.xyz', URL.createObjectURL(new Blob([text], { type: 'text/plain' })));
  };

  render(true);
  window.__molviewer = { addText, viewer, pick, setMode, measureTable,
    state: () => ({ files: files.map(f => f.name), atoms: atoms.length, picked, mode, measures: cur ? cur.measures : [], frame: cur ? cur.i : null }) };   // for debugging / automated checks
})();
