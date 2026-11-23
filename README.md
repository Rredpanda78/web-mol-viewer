# Web Mol Viewer

A small molecule viewer that runs entirely in the browser. Drop a structure file on the page and look at it — nothing is installed and nothing is uploaded.

**Live page:** https://rredpanda78.github.io/web-mol-viewer/

## What it does

- **Open files** by drag-and-drop or the file picker; several files can be open at once.
- **Formats:** `.xyz` (multi-frame), Tinker `.xyz` / `.txyz` / `.arc`, `.pdb`, `.sdf` / `.mol`, `.mol2`, `.pqr`, `.gro`, `.cif` (mmCIF), `.cube`, Gaussian output `.log` / `.out` (every optimisation step becomes a frame, with its SCF energy) and Gaussian input `.gjf` / `.com` (Cartesian).
- **Display:** ball and stick, stick, line, space filling, cartoon; element / atom-number labels; hide hydrogens; dark background; spin.
- **Measure:** pick the tool (Distance, Angle or Dihedral), then click 2, 3 or 4 atoms. Each measurement stays on screen and in a list, so you can keep several at once and remove them one by one. They follow the atoms when you change frame, and **Copy** puts them on the clipboard as a table with one row per frame (also for trajectories).
- **Polar contacts:** dashed lines and distances for H···N/O contacts below a cut-off (H bonded to N or O, acceptor at least three bonds away). Contacts between different molecules are blue, intramolecular ones red.
- **Charges:** per-atom charges as labels and as a red – white – blue colouring. Several charge sets can be loaded for the same structure and shown as a difference (set 2 − set 1); **Copy table** gives all sets side by side.
- **Trajectories:** add a `.dcd` file (OpenMM, NAMD, CHARMM) to an open structure. Frames are read on demand, so long trajectories are fine. In a solvated system the view shows the solute and the water molecules near it (periodic images are handled), and the solute stays centred.
- **Energy and contact plots:** add an energy table (`.csv`) and the energy is plotted against time under the viewer; click the plot to jump to that frame. **Scan contacts** counts the H···N/O contacts in every frame (between molecules and intramolecular) and plots them; **Copy table** gives time, energy and counts per frame. Gaussian optimisations and scans get an energy-per-frame plot automatically.
- **Frames:** slider, play / pause, ← / → keys.
- **Export:** PNG picture of the view, or the atoms on screen as `.xyz`.

## Adding data to a structure

Open (or select) the structure first, then drop the extra file on the page.

| File | What it adds | Notes |
|---|---|---|
| `.dcd` | trajectory | must have the same atoms, in the same order, as the structure (use the topology written with the trajectory, e.g. `top.pdb`) |
| `.csv` with a header | energy per step | columns are found by name: `step`, `time…`, `potential…` / `energy…` |
| OpenMM force-field `.xml` | charges | multipole monopoles (`c0`), `NonbondedForce` charges or residue-template charges; residues are matched by name |
| GROMACS `.itp` / `.top` | charges | first `[ atoms ]` section |
| Tinker `.key` / `.prm` | charges | `multipole` / `charge` lines; the structure must be a Tinker `.xyz` / `.txyz` (atom types) |
| `.chg`, `.txt`, `.dat`, one-column `.csv` | charges | `El x y z q` lines, or one value per line (last numeric column) |

An `.xyz` file with a fifth column (`El x y z q`, e.g. DDEC6 `net_atomic_charges.xyz`) and `.mol2` / `.pqr` files bring their own charges.

When the charge file lists other elements than the structure, the page warns that the atom order probably differs.

## Run it locally

It is a static page, so any web server works:

```
python -m http.server 8000
```

then open http://localhost:8000. The page loads [3Dmol.js](https://3dmol.csb.pitt.edu/) from cdnjs, so the first load needs an internet connection.

## Layout

```
index.html         page
css/style.css      styling
js/parsers.js      structure formats -> frames (XYZ, Tinker, Gaussian)
js/trajectory.js   DCD reader (frame by frame) and energy tables
js/charges.js      charge files -> one value per atom
js/app.js          viewer, measuring, contacts, charges, frames, plot
samples/           small idealised example molecules
```

To add a structure format, write a function in `js/parsers.js` that returns frames (`{ data, format, comment, energy }`) and register its extension in `parse()`. To add a charge format, add a reader in `js/charges.js` and route its extension in `parse()` there.

## Ideas for later

- Overlay / align two structures and report the RMSD
- GROMACS `.xtc` / `.trr` trajectories
- Orbitals and densities from `.cube` files as isosurfaces
- Edit a dihedral and export the new geometry
- Vibrational modes from Gaussian frequency jobs
- Atomic dipoles as arrows

## 中文簡介

在瀏覽器裡看分子結構的小工具:把 `.xyz`、`.pdb`、Gaussian `.log` 等檔案拖進頁面就能看,不用安裝、檔案也不會上傳。可以量距離/鍵角/二面角、顯示 H···N/O 接觸(分子間藍色、分子內紅色)、顯示原子電荷(標籤 + 紅白藍著色,可比較兩組電荷的差)、播放 MD 軌跡(`.dcd`,溶液只顯示溶質與附近的水)、畫能量對時間的圖(點圖跳到該幀)、逐幀統計氫鍵數,並匯出 PNG 或 `.xyz`。

## License

MIT. Rendering is done by [3Dmol.js](https://github.com/3dmol/3Dmol.js) (BSD-3-Clause).
