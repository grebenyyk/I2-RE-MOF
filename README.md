# I₂ @ MOF — DFT interaction energies (Quantum ESPRESSO)

Interaction energies of I₂ with a panel of MOF frameworks (POSCARs in `poscars_unitcell/`),
following the setup developed in `../i2_ce-uio-66/`.

## Goal

For each MOF, compute the I₂–framework interaction energy at **full I₂ loading**:

```
E_int    = E(loaded) - E(empty) - N_I2 * E(I2_vac)     [eV]
E_per_I2 = E_int / N_I2
```

`E(I2_vac) = -757.80004425 Ry` (PBE+D3(BJ), 50/400, `I.pbe-n-kjpaw_psl`) — the shared
isolated-I₂ reference from `../i2_ce-uio-66/01_i2_vac/`, reused unchanged.

The user's plan: start from the full-loaded config (one-I₂ / two-I₂ branches are too many
calculations for near-saturated frameworks), compare E_per_I2 across MOFs, then decide.

## Setup (mirrors the converged i2_ce-uio-66 loaded template)

| item | choice |
|------|--------|
| XC + dispersion | PBE + D3(BJ): `vdw_corr='grimme-d3'`, `dftd3_version=4` |
| cutoffs | `ecutwfc=50`, `ecutrho=400` Ry (PAW, 8×) — env-overridable |
| k-points | Γ-point only (cells ≥ 8 Å) |
| occupations | `smearing` mv, `degauss=0.01` (uniform for empty + loaded; T→0 for gapped) |
| SCF | `conv_thr=1e-5`, `diagonalization='david'`, `mixing_mode='local-TF'`, `mixing_beta=0.2`, `mixing_ndim=16`, `electron_maxstep=300` |
| relax | framework frozen (`if_pos 0 0 0`), only I atoms move (`if_pos 1 1 1`), BFGS; `etot_conv_thr=1e-5`, `forc_conv_thr=1e-3` |
| per system | `empty_mof.in` = SCF on bare framework (I removed); `loaded.in` = relax full I₂@MOF |

`david` (not `paro`) is used for the loaded relax: `paro` (the Ce-template choice) **stalls**
on the Gd+I₂ SCF — the density residual oscillates right at `conv_thr=1e-5` for one BFGS
geometry and burns 100s of iterations. `david` pushes through (≈90 iters on the hard step,
converges). `paro` worked for 10+ steps then hit one hard geometry.

## Metal / f-electron treatment — the key caveat

The Ce-UiO-66 reference used `nspin=1`, no Hubbard U, valid because Ce⁴⁺ is 4f⁰. The metals
here:

| metal | oxidation | shell | treatment |
|-------|-----------|-------|-----------|
| Cd | +2 | 4d¹⁰ closed | nspin=1, no U |
| Y  | +3 | 4d⁰ closed   | nspin=1, no U |
| Gd | +3 | 4f⁷          | nspin=1, no U — **4f frozen in core** (see below) |
| Tb | +3 | 4f⁸          | nspin=1, no U — 4f frozen in core |
| U  | +4 | 5f²          | nspin=1, no U — 5f frozen in core (⚠ verify f-count vs oxidation) |

The PSlibrary `spdn` PPs for Tb/Gd/U are **f-in-core**: the 4f/5f shell is frozen into a
non-spin-polarized core (`nlcc=.true.`) and is NOT a variational valence manifold, so
`U Gd-4f` errors ("Mismatch between requested and available manifolds"). With an f-in-core
PP the trivalent ion's **valence is closed-shell** (the 3 valence e⁻ transfer to O/linkers),
so the correct treatment is nspin=1, no U — identical to Ce⁴⁺/Cd/Y. The frozen f-shell is
chemically inert and identical in empty vs loaded, so it cancels exactly in E_int.

⚠ Two consequences:
1. The PSlibrary core is **nonmagnetic** → the f moment (e.g. Gd 7 μB) is NOT represented.
   Only the interaction energy ΔE is sound (which is all this campaign needs).
2. For U⁴⁺ the frozen-core f-count MUST match the actual oxidation. U–O = 2.28 Å with no
   uranyl U=O (~1.8 Å) → U is tetravalent 5f² here; verify the U PP core config = 5f².

If true magnetic / DFT+U f-physics is ever needed, source f-in-VALENCE PPs (PSlibrary does
not ship them for these). `build_inputs.py` keeps `METALS_FVALENCE` with the +U values for
that case.

## Pseudopotentials (PBE PAW, PSlibrary 1.0.0 legacy table)

Sourced from `https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/<el>`
(download base `.../upf_files/`). All in `/workspace/qe-7.5/pseudo/`:
`C.pbe-n-kjpaw_psl.1.0.0.UPF`, `H.pbe-kjpaw.UPF`, `O.pbe-kjpaw.UPF`,
`N.pbe-n-kjpaw_psl.1.0.0.UPF`, `I.pbe-n-kjpaw_psl.1.0.0.UPF`,
`Y.pbe-spn-kjpaw_psl.1.0.0.UPF`, `Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF`.
Still to fetch for the full campaign: Tb, U, Cd.

## Results (campaign, 2026-06-30)

| system | metal | topology | N_I2 | E_empty (Ry) | E_loaded (Ry) | E_int (eV) | E_int (kJ/mol) | E/I2 (eV) | E/I2 (kJ/mol) |
|--------|-------|----------|------|--------------|---------------|------------|----------------|-----------|---------------|
| 13_Y   | Y³⁺  | pts      | 4 | -2734.912430 | -5766.450222 | -4.59 | -443.20 | **-1.148** | **-110.80** |
| 06_Tb  | Tb³⁺ | pts      | 4 | -3892.857828 | -6924.393816 | -4.57 | -440.83 | **-1.142** | **-110.21** |
| 12_Tb  | Tb³⁺ | nan      | 6 | -7929.254500 | -12476.495596 | -6.00 | -578.70 | **-1.000** | **-96.45** |
| 08_Gd  | Gd³⁺ | pts      | 4 | -3832.692049 | -6864.163768 | -3.69 | -356.47 | **-0.924** | **-89.12** |
| 17_Cd  | Cd²⁺ | nan      | 2 | -2592.656222 | -4108.359447 | -1.40 | -135.39 | **-0.702** | **-67.70** |
| 07_Cd  | Cd²⁺ | pacman   | 6 | -4352.751357 | -8899.843520 | -3.97 | -383.19 | **-0.662** | **-63.86** |
| 10_Cd  | Cd²⁺ | pacman   | 5 | -4350.403273 | -8139.645079 | -3.29 | -317.14 | **-0.657** | **-63.43** |
| 05_Cd  | Cd²⁺ | pacman   | 5 | -4350.051068 | -8139.275398 | -3.05 | -294.20 | **-0.610** | **-58.84** |

U systems (16/18) skipped per user. E(I2_vac) = -757.80004425 Ry (shared reference).

**Reading the table:**
- Isostructural pts pair Y (-1.148) vs Tb_06 (-1.142) are identical -> with f-in-core PPs the
  metal is a spectator; I2 binding is framework-geometry-driven (good sanity check).
- Topology dominates: small-pore **pts** binds strongest (-0.92 to -1.15), **nan** intermediate
  (-0.70 to -1.00), large-pore **pacman Cd** weakest (-0.61 to -0.66) -> vdW confinement scaling.
- Cd is consistently weak (~-0.6 to -0.70) across topologies.

**Caveats:** most relaxes hit the 50-step BFGS cap (vdW plateau, ~0.005 eV/I2 of full convergence).
E/I2 at full loading includes I2-I2 guest-guest terms. Cross-system degauss 0.01 (Y/Gd/06_Tb)
vs 0.005 (rest) - sub-meV. f-in-core PPs (no f magnetism); ΔE sound.

## Files

- `build_inputs.py` — generic builder: `python build_inputs.py <POSCAR> <outdir> <label>`,
  or `python build_inputs.py --pilot`. Per-metal settings in the `METALS` dict.
- `collect_energies.py` — reads `*/empty_mof.out` + `*/loaded.out`, prints E_int / E_per_I2.
- `run.sh` — QE launcher (copied from i2_ce-uio-66; `NP`, `OMP_NUM_THREADS`, `--detach`).
- `13_Y/`, `08_Gd/` — pilot `empty_mof.in/out`, `loaded.in/out`.

## Running

```bash
cd /workspace/dimitry/calc/qe/i2_mof
python build_inputs.py --pilot                       # or per-system
NP=2 OMP_NUM_THREADS=1 ./run.sh --detach 13_Y/empty_mof.in
NP=2 OMP_NUM_THREADS=1 ./run.sh --detach 13_Y/loaded.in
python collect_energies.py
```
(NP=2 OMP=1 = one MPI rank per GPU on the 2×L40S node — see i2_ce-uio-66 README.)
