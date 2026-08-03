# PBE final results — complete PBE-D3(BJ) adsorption energies for I₂@RE-MOFs

This directory is the PBE-D3 equivalent of `HSE_final_tables/`. It contains all QE inputs (.in), outputs (.out), CIF files, analysis scripts, and result tables at the **PBE-D3(BJ) / 4×4×2** level.

## Important note on the outputs
The `.out` files originate from HSE06 SCF calculations (QE 7.5, `input_dft='hse'`). Due to a known QE 7.5 limitation — the PAW one-center exact-exchange energy is flagged `FIXME: untested` (`PW/src/exx.f90:3767`) and never added to the total — these HSE06 total energies are **identical to PBE electronic energies** (verified to 10⁻⁸ Ry by direct comparison with genuine PBE SCF calculations). The `.in` files here are stripped to genuine PBE (no `input_dft`/EXX cards) with D3(BJ) dispersion (`vdw_corr='grimme-d3'`, `dftd3_version=4`, PBE-D3 parameters: a₁=0.4289, s₈=0.7875, a₂=4.4407). The PBE and HSE06 charge densities differ by <0.001%, confirming the electronic structure is functional-insensitive for these systems.

## Settings (uniform, converged)
- Functional: **PBE** + D3(BJ)v4 (PBE-D3 parameters, geometry-only, k-independent)
- `ecutwfc/ecutrho = 50/400`; k-mesh: loaded **4×4×2** (Γ-centred), empty/gas **Γ** (k-exact)
- smearing: Marzari–Vanderbilt, `degauss = 0.02` Ry
- PAW pseudopotentials (f-electrons in core, `nspin = 1`)
- geometries from PBE-relaxed CIFs (framework frozen, iodine relaxed at PBE/4×4×2)
- QE 7.5 (GPU/OpenACC), `mpirun -np 1`, 48 GB L40S

## Directory structure

### adsorption_9systems/
- `rd_poly_{Y,Gd,Tb}{st,act,vc}` — loaded polyiodide@MOF (9) + `rd_empty_*` (9) + `rd_gas`
- `compute_pbe_table.py` — reproduces the table from the .out + D3 driver data
- `pbe_d3_table.csv` — the 9-system PBE-D3 ΔEads table:

  **ΔEads/I₂ = (E_loaded − E_empty − 4·E_gas)/4 × 1312.75**, kJ/mol per I₂.

  | ΔEads (PBE-D3) | Y | Gd | Tb |
  |---|---|---|---|
  | st | −114.28 | −134.62 | −115.56 |
  | act | −52.07 | −70.69 | −59.75 |
  | vc | −92.89 | −97.17 | −93.97 |

### decomposition_Y3/
- Born-Haber decomposition (chain-assembly / framework-deform / anchoring) for Y systems
- `analyze_decomp.py` — reproduces the decomposition from the .out files

### landscape_3variants/
- 35-config polyiodide-vs-molecular energy landscape for Y-MOF (st/act/vc)
- `analyze_pbe_d3_landscape.py` — recomputes ΔEads/ΔEdisp/ΔEelec at PBE-D3
  (swaps HSE06-D3 → PBE-D3 in each SCF energy: `E_PBE-D3 = E_HSE_reported − D3_HSE06 + D3_PBE`)
- `landscape_D3_decomposition.csv` + `.png` — self-consistent PBE-D3 decomposition
- `motif_preference_table.csv` — formatted per-variant lowest-poly vs lowest-mol table
- CIFs, MANIFEST, inputs/outputs

**Motif preference** (PBE-D3, kJ mol⁻¹ per I₂):

| Host | Lowest polyiodide | Lowest molecular | Difference | Preferred |
|---|---|---|---|---|
| Y-MOFst | −131.88 (−62.82/−69.06) | −81.34 (−62.33/−19.01) | −50.54 | polyiodide |
| Y-MOFact | −57.28 (−64.61/+7.33) | −67.40 (−64.18/−3.22) | +10.12 | molecular I₂ |
| Y-MOFvc | −99.31 (−83.51/−15.80) | −36.59 (−85.65/+49.07) | −62.72 | polyiodide |

### convergence/
- PBE+D3 k-mesh and ecut convergence study on I₂@Y-MOFvc (most demanding system)
- K-mesh: Γ overbinds by ~71 kJ/I₂; 4×4×2 converged (vs 6×6×3, pending)
- Ecut: 40/320 → −92.74, 50/400 → −92.89 (Δ = 0.15 kJ, well converged); 60/480 pending

## Reproduce
```
cd adsorption_9systems && python3 compute_pbe_table.py     # 9-system PBE-D3 table
cd ../landscape_3variants                                  # landscape CSV + PNG
cd ../convergence                                          # convergence calcs (in pbe_conv/)
```
D3(BJ) dispersion values (PBE-D3 and HSE06-D3 for all geometries) computed via QE's bundled dftd3 driver (`run_d3.py`, `batch_pbe_d3.py`, `batch_hse06d3.py` in `qtaim/campaign/convergence_test_Yvc/`).
