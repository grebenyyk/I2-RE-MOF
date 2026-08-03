# HSE final tables — exact QE inputs/outputs

This directory holds the **exact QE (pw.x) inputs and outputs** used to compute:
1. the **9-system adsorption table** (st/act/vc × Y/Gd/Tb),
2. the **chain-form/anchoring/deformation decomposition** for the three Y systems, and
3. the **polyiodide-vs-molecular energy landscape** for Y-MOF (35 configs across st/act/vc) with D3/HSE06 binding-energy decomposition.

## Settings (uniform, converged — locked by the Y-vc sweep)
- Functional: **HSE06** (`input_dft='hse'`, `exx_fraction=0.25`, `screening_parameter=0.106` Bohr⁻¹ = 0.2 Å⁻¹) + D3(BJ) v4
- `ecutwfc/ecutrho = 50/400`; `ecutfock = 100`
- k-mesh: **loaded 4×4×2**; empty/gas **Γ** (insulator/molecule, k-exact); `nqx1-3` default (= k; trivially converged)
- smearing: Marzari–Vanderbilt, `degauss = 0.02` (uniform everywhere)
- f-electrons in the pseudopotentials (Gd/Tb `*-spdn-*`), `nspin = 1`
- geometries: from the relaxed CIFs (`final_activated_results/{solvent-induced,activated,vc-relaxed}/{Metal}_{st,act,vc}_{loaded,empty}.cif`); gas I₂ in a 20 Å box at 2.71 Å
- QE 7.5 (GPU/OpenACC), `mpirun -np 1`, 48 GB L40S
- loaded metallic calcs print `! total energy` then may OOM/hang at the PAW-EXX post-step (benign — energy is captured before)

## adsorption_9systems/
- `rd_poly_{Y,Gd,Tb}{st,act,vc}` — loaded polyiodide@MOF (9)
- `rd_empty_{Y,Gd,Tb}{st,act,vc}` — empty-host reference Host⁰ (9)
- `rd_gas` — gas-phase I₂ reference
- `analyze_redo.py` — reproduces the table:

  **ΔEads/I₂ = (E_loaded − E_empty − 4·E_gas)/4 × 1312.75**, kJ/mol per I₂.

  | ΔEads | Y | Gd | Tb |
  |---|---|---|---|
  | st | −130.08 | −151.22 | −131.74 |
  | act | −69.64 | −89.12 | −77.66 |
  | vc | −111.98 | −117.18 | −113.71 |

- **ΔEads = ΔED3 + ΔEHSE06 decomposition (HSE06-D3, functional-matched)** — `compute_table1.py` /
  `generate_table1_artifacts.py` → `table1_HSE06D3.csv` + `.docx`. ΔED3 is the D3(BJ) dispersion
  evaluated at the **HSE06-D3(BJ)** level (QE maps `hse→hse06`; BJ params a1=0.383, s8=2.310, a2=5.685),
  via the bundled dftd3 driver (`d3_pbe_vs_hse06.json`). ΔEads is unchanged; only the split is corrected
  (the previous ΔED3 used PBE-D3, which understated dispersion by ~17 kJ/I₂). Values, kJ/mol per I₂:

  | ΔED3 / ΔEHSE06 | Y | Gd | Tb |
  |---|---|---|---|
  | st | −79.80 / −50.28 | −79.50 / −71.72 | −79.30 / −52.44 |
  | act | −82.04 / +12.40 | −80.72 / −8.40 | −81.42 / +3.76 |
  | vc | −102.55 / −9.43 | −104.51 / −12.67 | −101.12 / −12.59 |

  Corrected picture: st is anchoring-driven (ΔEHSE06 −50 to −72); **act is dispersion-dominated with
  ΔEHSE06 ≈ 0 to +12 (repulsive electronic — the act pore does not bind polyiodide electronically)**;
  vc is mostly dispersion (ΔEHSE06 ≈ −10), and rejects molecular I₂ electronically (not shown).


## decomposition_Y3/
- `dc_isubl_Y{st,act,vc}` — iodine sublattice I* at loaded geometry, no framework (3)
- `dc_hostld_Yvc` — framework at loaded geometry (Host*), no I, vc only (st/act: Host* = Host⁰ = `rd_empty`, deformation = 0)
- also includes the shared Y calcs it reads (`rd_poly_Y{st,act,vc}`, `rd_empty_Y{st,act,vc}`, `rd_gas`) so the dir is **self-contained**
- `analyze_decomp.py` — reproduces the Born–Haber decomposition:

  ```
  ΔE_chain-assembly      = (E*_I      − 4·E_gas)/4
  ΔE_framework-deform    = (E*_host   − E⁰_host)/4     [=0 for st/act]
  ΔE_anchoring           = (E_loaded  − E*_host − E*_I)/4
  ΔE_ads                 = chain-assembly + framework-deform + anchoring
  ```

  | System | ΔE_chain-assm | ΔE_deform | ΔE_anchoring | ΔE_ads |
  |---|---|---|---|---|
  | Y-st | +2.15 | 0.00 | −132.23 | −130.08 |
  | Y-act | −7.90 | 0.00 | −61.74 | −69.64 |
  | Y-vc | −31.36 | +42.77 | −123.40 | −111.98 |

## landscape_3variants/
The **polyiodide-vs-molecular energy landscape** for **Y-MOF** across the three pore variants (st/act/vc): 35 configurations — the authoritative CIF polyiodides plus perturbed/clustered polyiodide alternatives and void-placed/seeded molecular I₂ — all evaluated at HSE06-D3/4×4×2 on PBE/4×4×2-relaxed geometries (framework frozen). Establishes which motif each pore prefers and isolates the dispersion vs electronic contributions.
- `rd_poly/rd_mol/hmol/hse_lmol/hse_lpoly_*.in/.out` — the 35 config calcs
- `rd_empty_Y{st,act,vc}`, `rd_gas` — references
- `MANIFEST.md` — every config: E (Ry), ΔEads (kJ/I₂), and a **Lowest** column marking the 6 minima
- `*.cif` — final geometry of all 35 configs; the 6 lowest also copied as `<name>__LOWEST.cif`
- `landscape_energy_levels.png` — energy-level diagram (3 pore ticks × ΔEads; CIF polyiodides highlighted)
- `landscape_D3_decomposition.{csv,png}` — ΔEads split into **HSE06-D3** dispersion + HSE06 electronic. ΔEdisp is functional-matched (HSE06 BJ a1=0.383/s8=2.310/a2=5.685, via `analyze_d3_decomp_hse06.py` + `hse06d3_energies.json` from the bundled dftd3 driver); consistent with the 9-system Table 1 (st −79.8, act −82.0, vc −102.6). (An earlier PBE-D3 version understated dispersion by ~17 kJ/I₂; see `analyze_d3_decomp.py`.)
- `landscape_writeup.docx` — motivation / method / results / decomposition / justification
- `binding_decomposition.png`, `anchoring_vs_chargetransfer.png` — supporting figures (st/act/vc disp+elec; anchoring vs charge-transfer)

**Motif preference** (lowest polyiodide vs lowest molecular, kJ mol⁻¹ per I₂):

| Variant | lowest polyiodide | lowest molecular | polyiodide − molecular | preferred motif |
|---|---|---|---|---|
| Y-MOF st | −148.21 | −97.05 | +51.16 | polyiodide |
| Y-MOF act | −74.93 | −86.31 | −11.38 | molecular I₂ |
| Y-MOF vc | −118.39 | −55.11 | +63.29 | polyiodide |

**Binding-energy decomposition** ΔEads = ΔEdisp + ΔEelec (lowest configs):

| Variant | lowest polyiodide (ΔEdisp / ΔEelec) | lowest molecular (ΔEdisp / ΔEelec) | regime |
|---|---|---|---|
| st | −148.21 (−79.15 / −69.06) | −97.05 (−78.04 / −19.01) | anchoring-driven |
| act | −74.93 (−82.26 / +7.33) | −86.31 (−83.09 / −3.22) | dispersion-dominated (electronic ≈ 0) |
| vc | −118.39 (−102.59 / −15.80) | −55.11 (−104.18 / +49.07) | polyiodide weak electronic; molecular repulsive |

ΔEdisp (HSE06-D3) is a pore property (st ≈ act ≈ −80; vc ≈ −103), motif-independent; ΔEelec is the discriminator — ~zero-to-repulsive for the barely-charged act polyiodide (+7 to −3) and **repulsive (+49)** for vc molecular I₂.

## Reproduce
```
cd adsorption_9systems && python3 analyze_redo.py        # 9-system ΔEads table
cd adsorption_9systems && python3 compute_table1.py      # Table 1 split ΔEads=ΔED3+ΔEHSE06 (HSE06-D3)
# HSE06-D3(BJ) energies via the bundled dftd3 driver: run_d3.py / batch_9sys_d3.py / batch_hse06d3.py
#   (test.x in qe-7.5/dft-d3, func=hse06 version=4; validated: func=pbe reproduces the old PBE-D3 values)
cd ../landscape_3variants  # data + figures; build/analyze scripts live in
# qtaim/campaign/convergence_test_Yvc/ (analyze_d3_decomp_hse06.py = HSE06-D3 split, build_landscape_deliverables.py)
```
