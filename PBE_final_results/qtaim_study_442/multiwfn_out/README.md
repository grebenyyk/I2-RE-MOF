# QTAIM plane maps at the 4×4×2 k-mesh (Multiwfn)

Publication masters for the 442 CP study: electron-density contours + gradient-line
bond paths + critical-point dots, rendered by **Multiwfn** (2026.7.11, DISLIN build,
headless via xvfb) from the existing valence cubes of the density_442 campaign.
No new QE runs.

## Files
| file | panel | system | plane (cube atoms) |
|------|-------|--------|--------------------|
| `mwfn_442_st.png/.pdf`     | (a) | Y-MOFst(poly)  | I1, I2, I3(+a1 image) — atoms 73, 74, 315 |
| `mwfn_442_vc.png/.pdf`     | (b) | Y-MOFvc(poly)  | I1, I2, I5 — atoms 73, 74, 77 |
| `mwfn_442_act.png/.pdf`    | (c) | Y-MOFact(poly) | I1, I2, I3 — atoms 73, 74, 75 |
| `mwfn_442_actmol.png/.pdf` | (d) | Y-MOFact(mol.) | I1, I2, I3 — atoms 73, 74, 75 |
| `QTAIM_442_two_by_two.png/.pdf` | (a)–(d) composite, 2×2, letters | | |

Cube atom order reproduces the table labels: **atoms 73–80 = I1…I8** in all four
cells (every `BCP_II_full_442` d_pair value reproduced from the cube).

## What the maps show
Each map: black density contours (plot-0 valence field, as in the Γ-era figures),
gradient-line trajectories (bond paths), blue dots = (3,−1) CPs, blue sticks =
Multiwfn's empirical bonds, red element labels for atoms/CPs within 0.5 bohr
(0.26 Å) of the plane, axes in Å. Extension 10 bohr (7 for actmol).

## Corrected motif picture (important!)
Periodic-graph analysis of the table bonds (cell offsets resolved) shows that
**three of the four variants are infinite chains, not finite units**:

- **st** — infinite zigzag chain, period I1–I2–I3–I4 along a1: strong bonds
  2.88/3.01/2.90 Å within the period, weak I4···I1′ links (4.05 Å) and I1···I3′
  (4.41 Å) between periods. The slice shows one full period + the weak link.
- **vc** — TWO parallel zigzag chains (I4-I2-I3 and I7-I6-I8 periods, 2.90–3.06 Å)
  ~2.4 Å apart, with a molecular I2 (I1–I5, 2.795 Å) threaded BETWEEN them; each
  I end makes a borderline 3.40/3.43 Å contact to the opposite chain (I1···I2,
  I5···I6), so the whole assembly is one connected infinite component. The
  earlier reading "I2 bridging two cyclic I3 units" was an in-cell projection
  artifact: the three "ring" bonds connect different periodic images.
- **act** — infinite branched chain, period I1(branch)-I2-I3-I4 along a1:
  backbone 3.13 (I2–I3), 2.80 (I3–I4), 3.16 (I4–I2′) with a 2.80 branch (I1–I2)
  at each I2; second chain I5-I6-I7-I8 equivalent. Images stay within 0.43 Å of
  the slice plane, so the chain continuation is visible.
- **act(mol)** — the only finite motif: 4 discrete I2 molecules (2.70 Å) plus
  I···O contacts and one true vdW C···I CP (C1–I6, 3.901 Å).

## CP dots (imported, validated)
Multiwfn cannot find CPs on a cube (grid search returns 0) — dots are the
critic2 CPs converted to Multiwfn `CPs.txt` (Bohr Cartesian, type 2) by
`scripts/mwfn_cps_442.py`:
sources `bcp_raw_{tag}.out` ∪ `refine_{tag}.out` (seed-point refinement),
deduped <0.15 Å; BCPs only; iodine-involving (nearest-2 nuclei, same
classification as the tables); H contacts dropped (the path-verified artifacts);
near-core filters d2/d1 < 2 and d1 < 2.6 Å (PAW I-4d shell); **ρ ≥ 0.003 a.u.**
Per-cell counts vs the tables: st 11 (= 8 I-I + 1 I-O + 2 Y-I; absent: the two
newton-only weak links ρ 0.0036–0.0038 and O2–I6 ρ 0.0017 — all BONDED=0 in the
tables), vc 15, act 14, actmol 16 — vc/act/actmol exactly reproduce the union of
the II, IO, IC and MI tables. CPs are replicated over the tiled supercell
(matches `tile_cube.py` copy placement; originals keep indices).

## Cubes and tiling
Rendered from `val_{tag}.cube` (pp.x plot_num=0, smooth pseudo-valence — the Γ-era
figure field; the dots come from the plot-17-based critic2 analysis, as in Γ).
Tiled with `tile_cube.py` (originals first, so plane atom indices are unchanged):
st `[−1,+3]×[−1,+1]×[−1,+1]` (atoms: 73, 74, 315 = I3@+a1 in this ordering),
vc/act `[−2,+2]×[−1,+1]×[−1,+1]`, actmol `[−1,+1]³`. Black map boundaries mark
the edge of the periodic box (as in the Γ-era masters).

## Regenerate
```bash
cd qtaim/campaign/density_442/mwfn_out
bash ../../../../PBE_final_results/qtaim_study_442/scripts/run_cps_442.sh   # CP sets
./render_442.sh val_st_big2.cube 73,74,315 CPs_st.txt mwfn_442_st.png       # png
./render_442.sh val_st_big2.cube 73,74,315 CPs_st.txt mwfn_442_st.pdf       # vector
python3 .../scripts/compose_442.py                                          # composite
```
Requires `xvfb-run` + `libxm4`; `settings.ini` is written per-run from the
Γ-era template with the graphformat of the requested extension. Do **not** add
post-menu presses 17/18 (label style/threshold) — this build toggles labels off;
the settings.ini defaults (`iatmlabtype=1`, `disshowlabel=0.5`) are what we want.
