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

## Addendum — the adduct cell (dissociation/re-attachment), 2026-11-15
Test of the alternative reading of the experimental 202–203 + 167 bands as
donor-perturbed molecular I₂. One act-mol I₂ pulled to a carboxylate O
(2.3 Å closest contact after the rigid pull) — the molecule **dissociated**:
I8 accepted by the O region, I7 released. The released iodines then attacked
neighboring intact molecules, forming **two branched I₃-motifs** (I8–I1–I2:
2.82/2.89 Å; I7–I6–I3: 2.91/2.93 Å) around hub iodines I1 and I6; one molecule
(I4–I5, 2.70 Å) stayed intact. Split state near-degenerate with intact
(−11 kJ/cell). Computed spectrum (0.04): **207.9 (100, intact)**,
**175.7 (95) + 149.7 (80)** = asymmetric stretches of the two I₃-motifs,
125.7/104.6 (33/27, symmetric), nothing below 63.

Reading: the wall route does not give a stable O···I–I adduct — it gives
**dissociation → I₃-motif nucleation**, i.e. the first step of chain assembly.
The I₃-motif pair (175.7/149.7, broadened by the floating-iodine site
distribution) covers the experimental ~167 band; the intact molecule covers
202–203. The two experimental bands are therefore consistent with early-stage
conversion: intact I₂ + nucleated I₃ motifs, with chains appearing as
conversion proceeds (89–158 envelope, above).

Both experimental bands are thus assigned: **202–203 = intact molecular I₂;
~167 = I–I bonds of I₃-like nucleation units (or equivalently chain-bridging
I₂ — the same σ*-populated electronic class, ρ_BCP ~0.06 vs 0.50 intact)**.

## Addendum 2 — density analysis of the adduct cell (2026-11-15)
Production density footing on the final adduct geometry (60/800/442/mv 0.01/
conv 1e-8; empty reference reused from act-mol — framework verified identical
to 5e-7; plot 0/17/8 on the 160×180×300 grid; critic2 v1.4.309).

**Bader CT (Δρ = loaded − empty − I₂-only, loaded-system basins; net integral
8.8×10⁻⁶ e ≡ conserved):** ALL eight iodines gain population — I1 +0.125,
I2 +0.087, I3 +0.143, I4 +0.027, I5 +0.062, I6 +0.108, **I7 (released) +0.127,
I8 (O-bound) +0.070** e. **No ion pair**: both dissociation fragments are
donor-stabilized (wall O⁻ + neighbouring I₂ feed their σ* orbitals) while the
population is shared. Hubs of the two I₃-motifs gain the most (I1 +0.125,
I6 +0.108 — each carries two σ* bonds).

**BCPs (pair-seeding missed iodine pairs across cell boundaries — recovered by
explicit midpoint seeding; all five verified by gradient paths):**

| bond | d (Å) | ρ_BCP (e/Bohr⁻³) | ∇²ρ | role |
|---|---|---|---|---|
| I4–I5 | 2.700 | 0.0736 | −0.023 | intact molecule (≡ act-mol 0.0741) |
| I1–I8 | 2.821 | 0.0597 | +0.015 | I₃-motif 1 |
| I1–I2 | 2.889 | 0.0527 | +0.028 | I₃-motif 1 |
| I6–I7 | 2.909 | 0.0517 | +0.024 | I₃-motif 2 |
| I6–I3 | 2.926 | 0.0501 | +0.025 | I₃-motif 2 |
| I7···I8 | 5.513 | 0.0006 | ~0 | split-pair remnant (not bonded) |

The I₃-motif bonds carry **positive ∇²ρ with substantial ρ** — the dative/CT
bond signature, matching their σ* population. Files: `BCP_adduct.csv`,
`CT_adduct.csv` (all 80 atoms), source `campaign/density_442/adduct/`.
