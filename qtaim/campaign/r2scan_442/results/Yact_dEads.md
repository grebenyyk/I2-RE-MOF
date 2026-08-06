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
