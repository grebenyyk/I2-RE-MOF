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

## Headline numbers (ρ in e/bohr³; E = −V/2 in kJ/mol)

| system | chain structure | core I–I | core ρ | link I–I | link ρ | Y···I | I···O |
|---|---|---|---|---|---|---|---|
| st (poly_lS5) | 2×(4-core + links) | 2.87–2.91 Å | 0.042–0.056 | 3.00–3.01 Å | 0.042–0.043 | 2 found | none <4.5 Å |
| vc (poly_lS3) | continuous chain | 2.90–2.93 Å | 0.039–0.053 | 3.03–3.43 Å | 0.020–0.041 | 1 | 3 (ρ 0.008–0.018) |
| act (poly_lS2) | 2×(I₂)₂ units | 2.80 Å | 0.062 | 3.13–3.15 Å | 0.033–0.034 | — | 3 (ρ 0.002–0.006) |
| act (mol_rd) | 4×I₂ | 2.70 Å (cores) | (see table) | none <0.02 | — | — | 3 (ρ 0.005–0.006) |

Consistent with the energetics/decomposition: act cores are the *strongest individual*
bonds (ρ 0.062, E=−V/2 ≈ 75 kJ/mol) but the act chain has fewer/weak links; vc's chain is
the most connected (8 bonded contacts incl. 3.4 Å long links at ρ≈0.02) — the cooperative
host. act(mol) has no contact above the 0.02 criterion except its molecular cores.

## Straight-line profile caveat (important for the methods text)

`rho_line_min` equals ρ_BCP to <3% only for bonds with **straight bond paths**:
vc 6/8 bonded I–I within 3% (median 0.1%), act(poly) 4/8, st 2/8 (median 31-41%).
For curved paths the straight line misses the saddle and the line minimum is NOT the BCP —
the critic2 ρ/∇²ρ values are authoritative. The Γ-era "profile reproduces ρ_BCP to four
significant figures" statement held for the Γ geometries and does NOT transfer to 442 as-is.

## Status / TODO
- CSV only; .docx versions can be generated with the old make_full_ii-style script on request.
- Multiwfn plane figures with CP overlays: pending (separate step).
- Optional refinement: re-run critic2 with explicit NEW seeds at the newton-recovered
  positions so all rows share the critic2-spline footing.
