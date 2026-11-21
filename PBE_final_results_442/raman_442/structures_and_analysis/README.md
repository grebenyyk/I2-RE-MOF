# 442 Raman campaign — structures, energies, mode assignments

Naming convention (binding across ALL files here, Table S11, and the energy
table): structures are **Y-MOFst(poly)**, **Y-MOFvc(poly)**, **Y-MOFact(poly)**,
**Y-MOFact(mol)** and the test cell **Y-MOFadduct(test)**; iodines are
**I1…I8** = the last 8 atoms of the 80-atom cell in cell order (= cube order =
order of the production tables). Identical labels are used in the CIF atom
blocks, in the vibration-mode table, and in Table S11.

## Contents
| file | what |
|---|---|
| `Y-MOFst_poly.cif` … `Y-MOFadduct_test.cif` (5 files) | the exact 80-atom geometries on which the 4×4×2 FD Raman spectra were computed (framework frozen at the 4×4×2-relaxed positions; iodines re-optimized at each variant's Raman production footing) |
| `energy_decomposition_442.csv` | frozen-fragment decomposition (ΔE_assembly / ΔE_anchor / ΔE_def / ΔE_ads, kJ per I₂; from `campaign/decomp_442/`) |
| `vibration_modes_442.csv` | per-band assignment: band position, relative intensity, contributing iodine atoms, description |

## Energy decomposition (kJ per I₂)
| System | ΔE_assembly | ΔE_anchor | ΔE_def | ΔE_ads |
|---|---|---|---|---|
| Y-MOFst(poly) (seed lS5) | −11.58 | −123.59 | −1.97 | −137.14 |
| Y-MOFact(poly) (seed lS2) | −10.61 | −48.20 | −1.24 | −60.05 |
| Y-MOFvc(poly) (seed lS3) | −28.76 | −109.92 | +39.63 | −99.06 |

Identity ΔE_assembly + ΔE_anchor + ΔE_def = ΔE_ads closes exactly in all rows.
Y-MOFact(mol) (ΔE_ads = −67.1 kJ per I₂; the act molecular minimum) and the
Y-MOFadduct(test) cell (no adsorption reference) have no decomposition rows.

## Vibration bands (computed, 4×4×2; intensities = relative activities)
See `vibration_modes_442.csv`. Headline bands:

| structure | band (cm⁻¹) | act | assignment (short) |
|---|---|---|---|
| Y-MOFst(poly) | 148.9 / 146.1 | 100 / 99 | antisymmetric stretches of the two zigzag chains |
| Y-MOFst(poly) | 131.7 / 131.5 | 82 | chain stretches on the short I3–I4 / I7–I8 bonds |
| Y-MOFst(poly) | 89.9 / 89.1 | 27 | inter-period link elongations (I6–I7 / I2–I3) |
| Y-MOFvc(poly) | **167.8** | 100 | bridging I1–I5 molecule (2.79 Å) |
| Y-MOFvc(poly) | 140.2 / 137.8 | 95 / 93 | antisymmetric stretches of the two chain triangles |
| Y-MOFvc(poly) | 100.5 / 99.4 | 57 | secondary-bond antisymmetric stretches |
| Y-MOFact(poly) | 158.9 / 157.2 | 100 / 99 | branched-unit antisymmetric stretches (core vs links) |
| Y-MOFact(poly) | 129.5 / 126.7 | 93 / 92 | core-bond combinations |
| Y-MOFact(poly) | 90.4 / 89.4 | 49 | link antisymmetric stretches (I2 / I6 hubs) |
| Y-MOFact(mol) | 211.0–205.0 | ≈100 | the four intact-molecule I–I stretches |
| Y-MOFadduct(test) | 207.9 | 100 | intact control molecule I4–I5 |
| Y-MOFadduct(test) | 175.7 / 149.7 | 95 / 81 | antisymmetric stretches of the two branched I₃-type units |
| Y-MOFadduct(test) | 125.7 / 104.6 | 33 / 27 | their symmetric stretches |

Weak low-frequency modes (activity ≤ 5) are summarized in single rows
("collective librational/translational") rather than listed individually.

## Footings
FD displaced SCFs: 50/400 Ry, 4×4×2, conv 10⁻⁹; smearing — st 0.03, vc 0.02,
act / act(mol) / adduct 0.04 Ry (each Hessian internally consistent; see main
README). Geometries = the iodine sublattices re-optimized at these footings on
the 4×4×2-relaxed frameworks; framework identical to act-mol in the adduct cell
(verified 5×10⁻⁷).
