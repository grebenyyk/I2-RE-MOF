# Interactive QTAIM figure kit (I₂@Y-MOF, 4×4×2) — Multiwfn GUI

Everything needed to pick your own sections, shift them by Å, add/remove
critical points and render publication maps — self-service, no calculations.

## What's here
| file | what |
|------|------|
| `val_<tag>.cube` | electron-density field to plot (113 MB, unit cell, renders in seconds) |
| `rho17_<tag>.cube` | the field the CP tables were made from (only needed to regenerate dots) |
| `CPs_<tag>.txt` | critical points for the dot overlay (iodine-involving BCPs, ρ ≥ 0.003 a.u.; st 11, vc 15, act 14, actmol 16) |
| `bcp_raw_<tag>.out`, `refine_<tag>.out` | critic2 CP sources (for regenerating CPs after tiling) |
| `Y-MOF<tag>.cif` | structure for VESTA (atom order = cube order; **iodines = atoms 73–80 = I1…I8**) |
| `Y-MOF<tag>_CPs.cif` | same + He pseudo-atoms at the CP positions — open this in **VESTA** to see the dots in the lattice and choose your section |
| `tile_cube.py`, `mwfn_cps_442.py`, `tile_cps.sh` | make a chain-spanning "big" cube + matching dots locally |
| `settings_kit.ini` | Multiwfn settings for this workflow (isilent=0 → graph windows pop up; iuserfunc=−3 → fn 100 = your cube) |
| `gui_442.sh` | run the GUI on the server itself via `ssh -X` |

## Setup (pick one)
**A. On your own PC (smoothest):** install Multiwfn from
<http://sobereva.com/multiwfn/download.html> (Windows 64-bit build; free).
Copy `settings_kit.ini` into the Multiwtn folder as `settings.ini` (keep a
backup of the original). Copy the whole kit folder over (cubes ≈ 850 MB, or
just the ones you need). Run `Multiwfn.exe` from that folder.

**B. Zero install, on the server:** from an X11-capable connection
(`ssh -X`, or MobaXterm/VcXsrv on Windows), run
`./gui_442.sh val_st.cube CPs_st.txt` — each rendered map opens a native window.

## Making a map (the core loop)
```
val_st.cube          <- your cube (path; ENTER on the next prompt)
4                    <- plane map
100                  <- user-defined function = the loaded cube
6                    <- gradient lines + contours
300,300              <- grid
0                    <- set plane extension  (input e.g. 10, in Bohr; 10 bohr = 5.3 A)
4                    <- plane through 3 atoms
74,76,1              <- the three atoms (indices = cube/CIF order)
                     <- graph window pops up; close it to continue
8                    <- show bonds   -> then choose color 3 (blue)
-8                   <- switch axis units to Angstrom
0                    <- save graph (dislin.png in the working dir; rename it)
q                    <- exit
```
Load the CP dots first if you want them (once per session, before step 4):
`2` (topology) → `-4` (modify CPs) → `5` (load CPs.txt; ENTER = current dir) → `0` → `-10` (back).
Telltale it worked: the map's post-menu gains options `4` and `6`.

## Choosing the section
- **VESTA picker:** open `Y-MOFst_CPs.cif` in VESTA — iodines with their I1…I8
  labels plus He dots at the CPs. Rotate, decide which three atoms define the
  plane you want, then use those indices in Multiwfn (I_k = cube atom 72+k; CP
  dots in the CIF are labelled CP001…).
- Try the planes used for the delivered masters:
  st → `73,74,155` (I1, I2, I3 of the next chain period; on locally tiled cube)
  vc → `73,74,77` (I1, I2, I5 — I₂ + both borderline links, chains edge-on)
  act → `73,74,75`   actmol → `73,74,75`
- On the unit-cell cube, chain bonds that cross the cell edge show cut — use
  `./tile_cps.sh st` (etc.) to build `val_st_big.cube` + matching dots for the
  final render (indices 1–80 unchanged; on that cube st's third atom is 155).

## Moving the section by Å (and rotating)
In the plane menu (before or instead of re-picking atoms):
- `0` — extension distance (Bohr; 1 Å = 1.8897 bohr)
- `-1` — translation + rotation of the map: input the translation vector in
  Bohr (e.g. shift 1 Å along x: `1.8897 0 0`) and a rotation angle in degrees
- `6` — expert mode: type the plane origin and translation vectors directly

## Adding / removing critical points
With CPs loaded: `2` → `-4` →
- `-1` — print CP list with coordinates in Å (find the index you dislike)
- `2` — delete CPs by index (e.g. `5,7,12`)
- `3` — add a CP: x y z in **Bohr**, then its type (2 = bond CP)
- `4` — save the edited set back to `CPs.txt`
Then re-run function 4 — the dots update immediately.

## Rendering quality
- `graphformat= png|pdf` in settings.ini (pdf = vector master)
- post-menu `11` (gradient-line parameters), `3` (contour settings, e.g. more
  levels), `4` (CP dot size/colour) are safe to experiment with; **avoid 17/18**
  (label threshold/style) — in the 2026.7.11 build they can switch labels off.
- Happy with a plane but want it rendered on the big tiled cube (or at 300 dpi)?
  Just tell me "atoms i,j,k, ext N, shift ..." — rendering it here takes a minute.
