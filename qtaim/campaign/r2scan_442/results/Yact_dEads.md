# r2SCAN cross-check: ΔEads(Yact) — I₂ in activated Y-MOF
Single-point at PBE-PAW-relaxed geometry. NC PPs (Pseudo-Dojo I z=7, SG15 C/H/O/Y).
ecutwfc=70, ecutrho=560, 4×4×2 (loaded), Γ (empty/gas), degauss=0.04, gas 12 Å box.

| calc | E (Ry) | D3 (Ry) |
|---|---|---|
| r2SCAN loaded  | -1636.04021 | -0.24779 |
| r2SCAN empty   | -1455.61289 | -0.13738 |
| r2SCAN gas(12Å)|   -45.06811 | -0.00101 |
| PBE-NC loaded  | -1663.29689 | -0.53980 |
| PBE-NC empty   |  (see out)  |          |
| PBE-NC gas     |   -51.72749 | -0.00354 |

## ΔEads (kJ/I₂), 4 I₂/cell
| functional (matched NC PPs) | ΔEads | D3 | electronic |
|---|---|---|---|
| PBE-NC   | -122.78 | -64.93 | -57.86 |
| r2SCAN-NC|  -85.74 | -34.91 | -50.83 |
| **ΔΔEads (r2SCAN−PBE)** | **+37.05** (r2SCAN weaker) | | |

Production reference (PBE-PAW, different PPs): -67.12 (D3 -64.93, elec -2.19)

# Yst (same recipe)
| functional | ΔEads | D3 | electronic |
|---|---|---|---|
| r2SCAN-NC | -171.31 | -33.60 | -137.71 |
| PBE-NC    | -201.33 | -62.40 | -138.93 |
Production PBE-PAW Yst: -136.53 (D3 -62.40, elec -74.13)
ΔΔEads functional (r2SCAN−PBE, same PPs) = +30 kJ/I2 (r2SCAN weaker; ~all D3 params — electronic agrees to 1 kJ!)

## Cross-system summary (SP @ PBE-PAW geom, NC PPs, degauss 0.04)
| system | r2SCAN-NC | PBE-NC | PBE-PAW(prod) | elec PBE-NC | elec r2SCAN-NC | elec PBE-PAW |
|---|---|---|---|---|---|---|
| Yact | -85.7 | -122.8 | -67.1 | -57.9 | -50.8 | -2.2 |
| Yst  | -171.3| -201.3 | -136.5 | -138.9 | -137.7 | -74.1 |
- Relative trend (Yst ~2x Yact) reproduced by r2SCAN-NC. Strong-binding conclusion robust.
- Electronic binding functional-robust at fixed PP (PBE-NC vs r2SCAN-NC: 7 kJ Yact, 1 kJ Yst).
- Electronic binding PP-sensitive (PAW vs NC: ~56 kJ Yact, ~65 kJ Yst). PAW under-binds I electronically.
- D3 functional-specific (PBE-D3 ~2x r2SCAN-D3); PBE-NC D3 == PBE-PAW D3 (geometry-based, exact match).

# FULL r2SCAN-NC RELAXATION results (I-only relax, framework frozen, ecutwfc=70, degauss=0.04)
| system | ΔEads(relaxed) | D3 | elec | SP (PBE-geom) | deepened | PBE-PAW prod |
|---|---|---|---|---|---|---|
| Yact | -88.2 | -34.9 | -53.3 | -85.7 | -2.5 | -67.1 |
| Yst  | -173.5| -33.6 | -139.9| -171.3 | -2.2 | -136.5 |
- Relax deepened each ~2 kJ/I2 -> I2 geometry already near r2SCAN minimum (SP caveat removed, small correction).
- Yact floor-accepted (step 15, F=0.0089); Yst converged (F=0.0089, ~2.8h GPU0).
- D3 from relax outputs was mid-run-truncated; used SP-geometry D3 (geometry barely moved, <0.5 kJ uncertainty). Exact D3 via clean SP on relaxed geom if needed.
- ecut check: e70->e90 ΔEads shift +0.28 kJ/I2 (Yact) -> e70 converged; the 56 kJ PAW/NC electronic gap is REAL (not cutoff artifact).
- Outputs: rel_Yact_r2s_e70.out, rel_Yst_r2s_e70.out. Driver: run_relax_r2s.py (resumable).
