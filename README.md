# I2-RE-MOF — I₂ sorption in rare-earth MOFs (Quantum ESPRESSO)

Periodic DFT study of molecular iodine captured in rare-earth (Y, Gd, Tb) carboxylate MOFs,
with a smaller Cd/U/Cu screening set. The repository holds the full calculation tree: every
Quantum ESPRESSO input and output, the Python builders/collectors that generate and analyse
them, the relaxed structures (CIF/POSCAR), and the result tables and figures that go into the
manuscript.

Two questions drive the work:

1. **How strongly does I₂ bind, and why?** Adsorption energies per I₂ are computed for the
   9-system matrix {Y, Gd, Tb} × {st, act, vc} and decomposed into a dispersion (D3) and an
   electronic part.
2. **What is the guest actually?** Real-space analysis (QTAIM bond critical points, ELF,
   charge-transfer densities, PDOS, Raman) is used to decide whether the adsorbate is discrete
   I₂ or a linker-templated polyiodide, and what anchors it to the framework.

## Nomenclature

| tag | meaning |
|-----|---------|
| `st` | solvent-templated / strained pore — framework geometry as in the guest-loaded crystal |
| `act`| activated pore — framework relaxed empty, then loaded |
| `vc` | variable-cell — cell and framework relaxed with the guest present |
| `poly` / `mol` | polyiodide-chain vs discrete-molecular I₂ guest motif |
| `05_Cd`, `06_Tb`, `07_Cd`, `08_Gd`, `10_Cd`, `12_Tb`, `13_Y`, `17_Cd` | per-system dirs from the initial screening; the number is the structure ID in `poscars_unitcell/` |

## Repository layout

```
poscars_unitcell/          starting I₂-loaded unit cells (11 structures: Y, Gd, Tb, Cd, U, Cu)
05_Cd … 17_Cd              per-system screening trees: empty_mof.in/out, loaded.in/out,
                           activate_empty/, activated_i2*/ variants, seeded/uniform-chain runs
build_*.py, collect_*.py   input generators and energy collectors for the screening stage
run.sh, run_conv.sh, run_tight.sh   pw.x launchers (MPI/GPU, --detach)

qtaim/                     real-space analysis: tight SCF + pp.x cubes (rho, val, ELF),
                           critic2 BCP/Bader runs, difference densities, PDOS
qtaim/campaign/            the production campaign: per-system dirs (Y_str, Y_act, Y_actc,
                           Y_deep, Tb_*, Gd_*), Raman driver, ELF/PDOS plotting, HSE checks,
                           relax_442/, landscape_442/, vc_relax_442/, r2scan_442/,
                           paper_draft.md, consolidated_tables.txt
gd_neb/                    constrained-coordinate scans and NEB attempts for the
                           2 I₂ → I₄ chain-formation step in Gd-MOF (+ a cluster NWChem route)
cshm/                      continuous-shape-measure analysis of the metal coordination
                           polyhedra (cosymlib + spglib) for pristine/stripped/DFT geometries
i2_gas_raman/, i2_gas_raman.py   gas-phase I₂ finite-difference Raman, pipeline validation
i2_vac_ref/, verify_yfinal/      isolated-I₂ reference SCF and a geometry spot-check
hse/                       early HSE06 single points (empty / filled / I₂)

PBE_final_results/         production PBE-D3(BJ) results at Γ-relaxed geometries
PBE_final_results_442/     same, with geometries relaxed on the 4×4×2 k-mesh (current best)
HSE_final_tables/          HSE06-D3 tables (adsorption, Born–Haber decomposition, landscape)
final_activated_results/   packaged relaxed structures (st / act / vc-relaxed CIFs) + figures
```

Each `*_final_*` directory carries its own `README.md` with the exact settings, file naming and
the analysis script that regenerates its tables — read those before reusing numbers.

## Methods

- **Code:** Quantum ESPRESSO 7.5 (`pw.x`, `pp.x`, `projwfc.x`), GPU/OpenACC build; NWChem only
  for the cluster-NEB side experiment in `gd_neb/cluster_neb/`.
- **Functional:** PBE + D3(BJ) (`vdw_corr='grimme-d3'`, `dftd3_version=4`) for production.
  HSE06 and r2SCAN cross-checks live in `HSE_final_tables/` and `qtaim/campaign/r2scan_442/`.
- **Basis/cutoffs:** PAW (PSlibrary 1.0.0), `ecutwfc/ecutrho = 50/400` Ry for energetics,
  60/800 for the real-space (QTAIM/ELF) single points; norm-conserving PseudoDojo sets are used
  in the r2SCAN cross-check.
- **k-points:** loaded cells 4×4×2 (Γ-centred); empty frameworks and gas I₂ at Γ (k-exact).
  Γ-relaxed geometries were found to underbind by 6–26 kJ/I₂, which is why the 442-relaxed set
  supersedes them.
- **Smearing/spin:** Marzari–Vanderbilt, `degauss = 0.02` Ry in production; `nspin = 1` with the
  4f/5f electrons frozen in the Gd/Tb pseudopotential core (see the caveat below).
- **Relaxation:** framework frozen (`if_pos 0 0 0`), iodine free; `bfgs` or `damp` depending on
  how soft the polyiodide modes are. `vc_relax` runs for the `vc` variants.
- **Analysis:** `critic2` (v1.2.1005) for BCPs and Bader basins on the PAW-reconstructed valence
  cube (`plot_num=17`); charge transfer from δρ = ρ(loaded) − ρ(empty) − ρ(I-only); ELF from
  `plot_num=8`; PDOS from `projwfc.x`; Raman from a finite-difference Hessian (Δ = 0.01 Å) plus a
  bond-polarizability model; shape measures with `cosymlib`/`spglib`; structure I/O with `ase`.

**Metal / f-electron caveat.** The PSlibrary `*-spdn-*` pseudopotentials used for Gd and Tb are
f-in-core: the 4f shell is frozen into a non-spin-polarized core and is not a variational valence
manifold, so Hubbard U cannot be applied to it and no f moment is represented. With the trivalent
ion's valence closed-shell, `nspin = 1` and no U is the consistent choice, and the inert f shell
cancels between empty and loaded — energy differences are sound, f magnetism is not described.
`build_inputs.py` retains a `METALS_FVALENCE` table with U values should f-in-valence potentials
ever be used.

## Key results (as recorded in the result directories)

Adsorption energy per I₂, ΔEads = (E_loaded − E_empty − 4·E_gas)/4, in kJ mol⁻¹:

| level / geometry | Y-st | Gd-st | Tb-st | Y-act | Gd-act | Tb-act | Y-vc | Gd-vc | Tb-vc |
|---|---|---|---|---|---|---|---|---|---|
| PBE-D3, Γ-relaxed (`PBE_final_results`) | −114.28 | −134.62 | −115.56 | −52.07 | −70.69 | −59.75 | −92.89 | −97.17 | −93.97 |
| PBE-D3, 4×4×2-relaxed (`PBE_final_results_442`) | −136.53 | −138.82 | −138.76 | −67.12 | −76.41 | −73.52 | −99.26 | −103.58 | −100.01 |
| HSE06-D3, Γ-relaxed (`HSE_final_tables`) | −130.08 | −151.22 | −131.74 | −69.64 | −89.12 | −77.66 | −111.98 | −117.18 | −113.71 |

- Binding is **pore-geometry driven, not metal driven**: the isostructural Y/Gd/Tb series tracks
  together, and the strained (`st`) pore binds roughly twice as strongly as the activated one.
- The **dispersion/electronic split** (`PBE_final_results_442/adsorption_9systems/adsorption_table_442.csv`)
  puts ΔEdisp at −60 to −85 kJ/I₂ for all systems; the electronic term is what discriminates:
  −74 to −78 for `st`, ≈ −2 to −14 for `act`, ≈ −16 to −19 for `vc`.
- The **motif landscape** (18 configurations per Y variant, 4×4×2-relaxed,
  `PBE_final_results_442/landscape_3variants/motif_preference_table_442.csv`) finds polyiodide
  chains lowest in `vc` and `st`, with molecular I₂ nearly degenerate in `act` (1.1 kJ/I₂).
- **QTAIM/ELF/Bader** (`qtaim/campaign/consolidated_tables.txt`, `paper_draft.md`) show a bonded
  polyiodide with 0.3–0.9 e per cell transferred from the linker C/O to iodine while the
  rare-earth cation stays a spectator (|δq_M| < 0.04 e).
- **Cross-checks:** r2SCAN/NC reproduces the PBE-PAW totals (Yst within 3 kJ, Yact within 14 kJ,
  `qtaim/campaign/r2scan_442/results/Yact_dEads.md`); HSE06 does not open a clean gap, so the
  charge transfer is not a PBE self-interaction artefact (`qtaim/campaign/hse_results_summary.txt`).

## Running / reproducing

Everything is driven by plain QE inputs, so any `pw.x` build works; `run.sh` is a convenience
launcher and hard-codes the original machine's paths (`/workspace/qe-7.5`, local OpenMPI and
NVHPC) — edit `QE_ROOT`/`MPI_ROOT` for your system. Python scripts need `numpy`, `ase`, plus
`spglib`/`cosymlib` for `cshm/` and `matplotlib` for the plotting scripts. Several older scripts
also embed the absolute base path `/workspace/dimitry/calc/qe/i2_mof`.

Build and run one system from a loaded POSCAR:

```bash
python build_inputs.py poscars_unitcell/POSCAR_13_2020_Y_pts_3_ASR_1_I2loaded_1000Pa_UC 13_Y 13_Y
NP=2 OMP_NUM_THREADS=1 ./run.sh --detach 13_Y/empty_mof.in
NP=2 OMP_NUM_THREADS=1 ./run.sh --detach 13_Y/loaded.in
python collect_energies.py        # E_int and E_per_I2 from */empty_mof.out + */loaded.out
```

Regenerate the published tables from the archived outputs (no new SCF needed):

```bash
python3 PBE_final_results/adsorption_9systems/compute_pbe_table.py        # PBE-D3, Γ geometries
python3 PBE_final_results_442/adsorption_9systems/compute_442_table.py    # PBE-D3, 4×4×2 geometries
python3 HSE_final_tables/adsorption_9systems/analyze_redo.py              # HSE06-D3 table
python3 HSE_final_tables/decomposition_Y3/analyze_decomp.py               # Born–Haber decomposition
python3 PBE_final_results_442/landscape_3variants/analyze_l442.py         # motif landscape ranking
python3 cshm/cshm_i2mof.py                                               # coordination shape measures
```

## Where results live

- Tables and CIFs for the manuscript: `PBE_final_results_442/` (current), `PBE_final_results/`,
  `HSE_final_tables/`.
- Relaxed structures by pore state: `final_activated_results/{solvent-induced,activated,vc-relaxed}/`.
- Figures (ELF slices, PDOS, Raman, landscape levels, BCP plots): `qtaim/campaign/` and
  `qtaim/campaign/bcp_plots/`.
- Manuscript text in progress: `qtaim/campaign/paper_draft.md`,
  `qtaim/campaign/final_activated_results/raman_discussion.md`.

## Repository hygiene

`.gitignore` is a **default-deny allowlist**: everything is ignored, directories are traversed,
and only known-good text/code/result/figure extensions (`.in`, `.out`, `.py`, `.sh`, `.md`,
`.cif`, `POSCAR*`, images, …) can be committed. The working tree is hundreds of GB of QE
wavefunctions, charge densities, `.save/` dirs and cubes; none of that is tracked, and new
calculations are covered automatically without editing the ignore file.
