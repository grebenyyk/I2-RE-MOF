# Full vc-relax of I2@Y-MOFvc at 4×4×2 (2026-11-08)

Question: does the loaded vc cell, relaxed with **cell + ions free** at the
production 4×4×2 mesh, converge to something different from the Γ-converged
geometry that all 442 work has used (cell held frozen)?

Answer: **yes — the same motif in a meaningfully different cell, worth
−5.46 kJ per I₂.** The Γ geometry is not a 442 stationary point.

## Runs (this directory, real files)

| file | what |
|---|---|
| `vcrelax442_Yvc.in/.out` | vc-relax from the Γ geometry (rd_poly_Yvc.in structure), bfgs/bfgs, `cell_dofree='all'`, press 0 / press_conv_thr 0.5 kbar, degauss 0.04 relax footing, 4×4×2; **converged in 61 bfgs steps (~10.5 h, 1×L40S)** |
| `sp02_vcrelax442_Yvc.in/.out` | clean degauss-0.02 SP of the relaxed geometry (same footing as the adsorption table) |
| `launch_vcrelax.sh`, `resume.sh` | launch / restart-from-.save |

## Lattice: 442-relaxed vs Γ seed

| | a (Å) | b (Å) | c (Å) | α | β | γ | V (Å³) |
|---|---|---|---|---|---|---|---|
| Γ seed | 8.7975 | 8.3276 | 17.7257 | 90.68 | 103.44 | 99.93 | 1242.3 |
| 442 vc-relax | 8.8126 | 8.4240 | 17.7370 | 90.55 | 101.27 | 98.39 | **1276.6** |
| Δ | +0.17% | **+1.16%** | +0.06% | −0.13° | **−2.17°** | **−1.54°** | **+2.77%** |

Shape change, not size change: the expansion is carried by b and by β/γ shear;
a and c are stationary. Convergence path: a grew first (+0.67% at bfgs step 11),
then partially reversed while b kept growing — the polyiodide sublattice
redistributes within the expanding cell.

## Displacements (fractional, MIC, in the final cell)

- framework RMSD **0.105 Å** (max 0.23)
- iodine RMSD **0.249 Å** (max 0.42) — the chains do the adapting

## I–I motif: preserved topology, tightened bonding

Same connectivity and infinite-chain topology (two zigzag chains I2-I3-I4 /
I6-I7-I8 threaded by the molecular I1–I5; periodic images as in the motif
analysis). Bond lengths vs the frozen-cell 442 values:

| contact | frozen cell | vc-relaxed |
|---|---|---|
| I1–I5 (molecular I₂) | 2.795 | 2.798 |
| chain A: I2–I3 / I3–I4 / I2–I4 | 2.902 / 3.046 / 3.056 | 2.951 / 2.989 / 3.065 |
| chain B: I6–I7 / I7–I8 / I6–I8 | 2.925 / 3.058 / 3.030 | 2.951 / 2.991 / 3.062 |
| **borderline bridges I1···I2 / I5···I6** | 3.403 / 3.427 | **3.364 / 3.367** |
| weak cross-links (I1···I3, I5···I7, …) | 4.26–4.96 | 4.39–4.92 |

The two chains become near-degenerate (I2–I3 = I6–I7 = 2.951 Å), and the
b-axis growth pulls the threaded I₂ toward both chains: the borderline bridges
shorten by ~0.05 Å — the vc-relax **strengthens the defining interaction** of
this variant.

## Energy (clean 0.02 footing, 4×4×2)

| | total (Ry) | D3 (Ry) |
|---|---|---|
| frozen-cell winner (`sp_l442_vc_poly_lS3`) | −5767.31129446 | −0.61591957 |
| **vc-relaxed (this run)** | **−5767.32791737** | −0.60032776 |
| Δ (kJ/cell) | **−21.82** | +20.47 |
| Δ (kJ per I₂, 4 I₂/cell) | **−5.46** | +5.12 (elec: −10.57) |

The cell relaxation trades dispersion (weaker in the expanded cell) for a
2× larger electronic gain — consistent with the shortened borderline bridges.
At fixed empty/gas references, letting the vc cell relax buys **−5.46 kJ/I₂**
against the frozen-Γ-cell production footing (−99.06 → −104.5 kJ/I₂ equivalent
loaded-energy shift; the adsorption ΔEads itself also depends on the empty
reference, which we did not re-relax).

## Recommendation / status

Keep the frozen-Γ-cell campaign as the production protocol (uniform across all
variants; the difference is a geometry-only ±5.5 kJ/I₂ for vc), and cite this
number as the vc-variant cell-relaxation energy in the SI/limitations note:
*"a full 4×4×2 vc-relax of I2@Y-MOFvc changes the cell by +2.8% in volume
(predominantly β/γ shear and the b axis) and lowers the loaded energy by a
further 5.5 kJ per I₂, preserving the polyiodide motif."*

Caveat: the empty-cell reference (rd_empty_Yvc, Γ cell) was NOT re-relaxed at
442; per the empty-vc-relax study the empty vc cell sits in a metastable
compressed basin, so an empty-442-vc-relax could change the net adsorption
number in either direction. Not run.
