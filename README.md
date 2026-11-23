# Web Mol Viewer

A small molecule viewer that runs entirely in the browser. Drop a structure file on the page and look at it — nothing is installed and nothing is uploaded.

**Live page:** https://rredpanda78.github.io/web-mol-viewer/

## What it does

- **Open files** by drag-and-drop or the file picker; several files can be open at once.
- **Formats:** `.xyz` (multi-frame), Tinker `.xyz` / `.txyz` / `.arc`, `.pdb`, `.sdf` / `.mol`, `.mol2`, `.gro`, `.cif` (mmCIF), `.cube`, Gaussian output `.log` / `.out` (every optimisation step becomes a frame, with its SCF energy) and Gaussian input `.gjf` / `.com` (Cartesian).
- **Display:** ball and stick, stick, line, space filling, cartoon; element / atom-number labels; hide hydrogens; dark background; spin.
- **Measure:** pick the tool (Distance, Angle or Dihedral), then click 2, 3 or 4 atoms. Each measurement stays on screen and in a list, so you can keep several at once and remove them one by one. They follow the atoms when you change frame, and **Copy** puts them on the clipboard as a table with one row per frame (handy for scans and optimisations).
- **Polar contacts:** dashed lines and distances for H···N/O contacts below a cut-off (H bonded to N or O, acceptor at least three bonds away).
- **Frames:** slider, play / pause, ← / → keys.
- **Export:** PNG picture of the view, or the current frame as `.xyz`.

## Run it locally

It is a static page, so any web server works:

```
python -m http.server 8000
```

then open http://localhost:8000. The page loads [3Dmol.js](https://3dmol.csb.pitt.edu/) from cdnjs, so the first load needs an internet connection.

## Layout

```
index.html        page
css/style.css     styling
js/parsers.js     file formats -> frames (XYZ, Tinker, Gaussian)
js/app.js         viewer, measuring, contacts, frames
samples/          small idealised example molecules
```

To add a format, write a function in `js/parsers.js` that returns frames (`{ data, format, comment, energy }`) and register its extension in `parse()`.

## Ideas for later

- Overlay / align two structures and report the RMSD
- Plot energy against frame for optimisations and scans
- Orbitals and densities from `.cube` files as isosurfaces
- Edit a dihedral and export the new geometry
- Vibrational modes from Gaussian frequency jobs
- GROMACS / OpenMM trajectories

## 中文簡介

在瀏覽器裡看分子結構的小工具:把 `.xyz`、`.pdb`、Gaussian `.log` 等檔案拖進頁面就能看,不用安裝、檔案也不會上傳。可以選距離/鍵角/二面角工具後點原子量測(多筆同時保留,可複製成每幀一列的表)、顯示 H···N/O 接觸、播放多幀(最佳化每一步),並匯出 PNG 或 `.xyz`。

## License

MIT. Rendering is done by [3Dmol.js](https://github.com/3dmol/3Dmol.js) (BSD-3-Clause).
