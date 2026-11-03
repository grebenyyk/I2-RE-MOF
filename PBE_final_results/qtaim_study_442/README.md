# qtaim_study_442 — CP/BCP analysis at the 4×4×2 tight footing (2026-10-31)

The `PBE_final_results/qtaim_study/` pipeline applied to the density_442 systems
(Y-only; act includes the molecular baseline). Canonical search: critic2 v1.2.1005
**`auto seed pair dist 4.0 nr 8`** on the plot-17 cubes (`density_442/*/rho17_*.cube`;
raw outputs in `../../qtaim/campaign/density_442/*/bcp_raw_*.out`), Newton recovery of
missing I–I BCPs from bond midpoints (flagged `src=newton`, trilinear interpolation —
critic2 rows use critic2's spline; the two footings differ up to ~20-40% near cores,
far less mid-bond).

**New vs the Γ-era filter:** besides the near-core filter (d1 < 2.6 Å, d2/d1 < 2), a
**shell filter** removes (3,−1) points within 1.35 Å of an iodine nucleus: the
`I.pbe-n-kjpaw` pseudopotential carries 4d¹⁰ in valence, and its shell structure creates
spurious near-core critical points at ~1.0 Å that would otherwise masquerade as short
I–I "bonds". (5 removed in act(mol), 1-2 elsewhere.)

## Tables (tables/, CSV; same columns as the Γ-era tables + extras)

- `BCP_II_full_442.csv` — all I–I contacts < 4.5 Å per cell, `BONDED(0.02)` = the
  ρ ≥ 0.02 e/bohr³ criterion, `rho_line_min` = straight-line profile minimum (see caveat).
- `BCP_IO_full_442.csv`, `BCP_MI_native_442.csv` — I···O and Y···I.

## Validation (2026-10-31, after the cell-matrix fix)

- **Every BCP is bond-path verified**: gradient paths traced uphill from each CP reach
  exactly its two assigned nuclei (39/39 I–I, 17/17 I–O, 5/5 Y···I, 1/1 C···I excluded
  rows re-checked; `path_ok` column = 1 everywhere in the final tables). Two earlier
  "C–I bonds" in vc were exposed by this check as a C–H···I feature and a ring-like CP —
  removed.
- **All I–I and I–O distances cross-checked against the CIF geometries** (39/39, 17/17
  exact matches; `scripts/validate_tables_vs_cif.py`).
- Cell-matrix bug found and fixed: `H.T * dims` broadcasting corrupted distances in cells
  with fully general lattice vectors (vc); st/act were unaffected by luck of axis-aligned
  grid vectors. Same family as the old fractional→Cartesian transpose bug.

## Headline numbers (ρ in e/bohr³; E = −V/2 in kJ/mol)

| system | chain structure | core I–I | core ρ | link I–I | link ρ | Y···I | I···O |
|---|---|---|---|---|---|---|---|
| st (poly_lS5) | 2×(4-core + links) | 2.87–2.91 Å | 0.042–0.056 | 3.00–3.01 Å | 0.042–0.043 | 2 found | none <4.5 Å |
| vc (poly_lS3) | **I₂(2.795 Å) bridging two I₃ rings** | I1–I5 **2.795** ρ 0.064 | I₂ core | 2.90–3.06 (triangles), 3.40/3.43 (bridges) | 0.020–0.053 | 2 | 4 (ρ 0.005–0.018) |
| act (poly_lS2) | 2×(I₂)₂ units | 2.80 Å | 0.062 | 3.13–3.15 Å | 0.033–0.034 | — | 3 (ρ 0.002–0.006) |
| act (mol_rd) | 4×I₂ | 2.70 Å (cores) | (see table) | none <0.02 | — | — | 3 (ρ 0.005–0.006) |

Motifs from the (path-verified) connectivity graphs:
- **st**: two linear I₄ chains (2.88/3.01 alternating) — classic polyiodide.
- **vc**: an **I₂ molecule (I1–I5, 2.795 Å, ρ 0.064) bridging two cyclic I₃ units**
  (I2–I3–I4 and I6–I7–I8, edges 2.90–3.06) via 3.40/3.43 Å contacts (ρ 0.020–0.021,
  borderline at the 0.02 criterion — flag these as weak links, not chain bonds).
- **act**: two branched (I₂)₂ units (2.80 cores + 3.13–3.16 links).
- **act(mol)**: four isolated I₂; only molecular cores pass 0.02.

Consistent with the energetics/decomposition: act cores are the strongest individual bonds
but fewest links; vc's I₂+2×I₃ aggregate is the most connected — the cooperative host.

## Value footing — unified on critic2 splines (2026-11-03)

All ρ/∇²ρ values come from critic2's spline of the plot-17 cube:
`critic2` = found by the `auto seed pair dist 4.0 nr 8` search; `critic2(seedpt)` = CPs the
seeded search missed, refined by critic2 1.4.309 (`tools/critic2_new`, built from upstream
master) seeded at our Newton-recovered positions (`auto seed point x0 ...`, one seed per
line; the old 1.2.1005 binary lacks this grammar). Only 2 rows remain trilinear — st's two
non-bonded 4.4 Å fringe contacts (BONDED=0, information only). The `src` column tracks
provenance; `path_ok` = the gradient-path verification (all rows = 1).

`rho_line_min` (straight-line profile minimum) is an information column only: it equals
ρ_BCP within 3% only for straight bond paths (vc 6/9, act(poly) 4/8, st 3/7); for curved
paths the line misses the saddle. The Γ-era "profile reproduces ρ_BCP to four significant
figures" statement does not transfer to 442 — use bond-path verification instead.

## Status / TODO
- CSV only; .docx versions can be generated with the old make_full_ii-style script on request.
- Multiwfn plane figures with CP overlays: pending (separate step).
- Optional refinement: re-run critic2 with explicit NEW seeds at the newton-recovered
  positions so all rows share the critic2-spline footing.
