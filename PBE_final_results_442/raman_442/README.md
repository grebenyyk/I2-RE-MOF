# 442 Raman campaign — I-sublattice spectra of I₂@Y-MOF (2026-11-14)

Finite-difference phonons + bond-polarizability Raman for the four production
winner geometries (framework frozen, 8 iodines × 3 axes × ±0.01 Å = 48 displaced
SCFs per system; central differences; 4×4×2 mesh; QE 7.5 PBE+D3(BJ), ecut 50/400;
conv_thr 1e-9 geometry SCFs, 1e-8-9 displaced SCFs; factor 3634.77 cm⁻¹, validated
on gas I₂ = 214.5 vs exp 213).

## Deliverables
| file | what |
|---|---|
| `raman_modes_442.csv` | all 89 modes: system, smearing, mode#, ν, relative activity |
| `raman_spectra_442.png/.pdf` | four broadened spectra vs experimental band guides (orange 199–214 intact I₂; green 160–174 bridged I₂) |

## Results (per-system, smearing footing)
| system | smearing | main bands ν (activity) |
|---|---|---|
| st (zigzag chains) | 0.03 | 148.9/146.1 (100/99), 131.7/131.5 (82), 89.9/89.1 (27) |
| vc (chains + bridged I₂) | 0.02 | **167.8 (100, bridged I₁–I₅)**, 140.2/137.8 (95/93, chain primary), 100.5/99.4 (57, secondary) |
| act (branched chains) | 0.04 | 158.9/157.2 (100/99), 129.5/126.7 (93/92), 90.4/89.4 (49) |
| act(mol.) (molecular I₂) | 0.04 | **211.0/208.7/207.8/205.0 (≈100 each, intact I₂)**; nothing 40–200 |

## Assignments vs the experimental bands (202–203 and ~167 cm⁻¹)
- **202–203 → intact molecular I₂** (act-mol-like population; computed 205–211,
  smearing-insensitive since the system is gapped). The persistent band = the
  un-converted molecular fraction inside the pores.
- **~167 → chain-bridged I₂** (vc-type I₁–I₅): intact but σ*-populated by the
  borderline contacts — BCP density 0.064 a.u. vs 0.50 for intact I₂, Bader
  population +0.21 e, ν red-shifted 205→168.
- **st/act chain ladders have no bands in either window** — pure chain variants
  produce only the 89–158 envelope.
- Adduct test (see `campaign/raman_442/adduct_test/`): forcing I₂ into a
  carboxylate contact (2.3 Å seed) **dissociates** it (one I chemisorbed on O,
  one released; split state near-degenerate with intact, −11 kJ/cell) — the
  wall route cannot produce a 167-class *intact* I–I band. Its computed
  two-population spectrum (3 intact molecules + dispersed pair) is the
  quantitative check of the 202–203 assignment.

## Campaign history (why the footings differ)
Production smearing 0.02 deadlocked the SCFs of act and (partially) st — the
known act-floor behavior. act pair: FD sets at 0.04 (after 0.02/0.03 relax
floors). st: 0.02 set abandoned (4 displacements incl. I074_xn would not converge
below 0.05) → whole set regenerated at 0.03. Each Hessian is internally
consistent; cross-system comparisons carry the smearing footnote. FD linearity
validated: force-difference ratio at doubled amplitude = 2.0006 (ideal 2.0).
Nodes died twice mid-campaign; all queues resumed from `.done` markers.

## Reproduce
`campaign/raman_442/`: `gen_rerelax_442.py`, `gen_scf_loaded_442.py`,
`raman_driver_442.py` (gen/run/hessian/raman; RAMAN_GPUS, FD_DEGAUSS,
MAX_SCF_S env), `make_deliverables_442.py`, `adduct_test/`.
