# density_442 act pilot — real-space density at the 4x4x2 tight footing (2026-10-30)

First real-space study on the **442-relaxed winner geometries** (fixes both the Gamma-density
and the Gamma-geometry approximations of the old Y_actc-era study). Settings: PBE-D3(BJ),
ecut 60/800, degauss 0.01, conv_thr 1e-8, **K 4x4x2**, PAW plot 17 rho / plot 0 val / plot 8 ELF,
FFT 160x180x300 (0.056 A). Two systems = the act speciation contrast:
- `act_poly_lS2` (lowest polyiodide; cores 2.80 + links 3.13-3.16 A)
- `act_mol_rd` (molecular I2; 4x 2.70 A; act overall minimum)

## Charge transfer (Bader-basin delta-rho, e/cell, plot-17 basins, old verified 4th-column convention)

| system | C | H | O | metal | I | Sigma |
|---|---|---|---|---|---|---|
| poly_lS2 | -0.084 | -0.018 | -0.132 | -0.015 | **+0.249** | +0.0000 |
| mol_rd   | -0.009 | -0.003 | -0.042 | -0.015 | **+0.069** | +0.0005 |

**The polyiodide carries 3.6x more charge transfer than molecular I2 in the same host**;
the linker C/O donates in both cases and Y stays a spectator (-0.015) in both.
(old Gamma-geometry Y_actc: I +0.495 on a different basin/geometry.)

## 2x grid sensitivity check (320x360x600 Fourier-interpolated, rho17+deltarho, ct rerun)
Per-species delta-rho identical to 4 decimals (all d = +0.0000). The 0.056 A FFT grid is
converged for basin integration; the methods-text spacing note is safe.

## Topology
bcp (auto, OMP=1): poly_lS2 316 bond CPs; mol_rd 242. Full CP lists in bcp_<tag>.out.
NOTE: critic2 needs OMP_NUM_THREADS=1 on this node (libgomp team abort otherwise) and
CRITIC_HOME (tools/critic2_run.sh).

## Files
rs_{loaded,empty,i2only}_<tag>.{in,out} — tight SCFs (all converged);
pp_*.in/out + *.cube — cubes; ct/bcp _<tag>.cri/out — critic2; *_2x_* — sensitivity set;
deltarho_<tag>.cube = val(loaded) - val(empty) - val(i2only).
