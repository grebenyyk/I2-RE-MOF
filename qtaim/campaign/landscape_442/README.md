# Landscape 4×4×2 relaxation campaign (Y-only, 18 configs)

## Goal
Find the true lowest-energy iodine arrangement per pore variant (vc/act/st) at the
proper 4×4×2 relaxation. Supersedes the `landscape_3variants` Γ-relaxed-geometry
ranking, which undersamples the PES by 6–26 kJ/I₂.

## The 18 configurations
One representative per distinct I-I-fingerprint basin (dedup'd from 35 Γ-relaxed
seeds). Seeds in `PBE_final_results/landscape_3variants/`. Tag = `l442_<var>_<motif>_<snap>`:

| variant | tags |
|---|---|
| vc | vc_poly_S1, vc_poly_lS3, vc_poly_lS6, vc_mol_lS4, vc_mol_ref, vc_mol_lS5, vc_mol_lS3 |
| act | act_mol_rd, act_poly_lS2, act_poly_lS5, act_mol_lS2, act_mol_lS6 |
| st | st_poly_S1, st_poly_lS5, st_poly_lS6, st_poly_lS4, st_mol_lS6, st_mol_lS3 |

## Recipe (all 18, uniform)
- I-only relax: framework frozen (if_pos 0 0 0), 8 I free (1 1 1). Cell fixed per variant.
- ion_dynamics = damp; degauss = 0.04; conv_thr = 1e-6; etot_conv_thr = 1e-5; forc_conv_thr = 5e-3;
  ecutwfc=50; ecutrho=400; smearing=mv; K 4×4×2; davidson; mixing local-TF β=0.2 ndim=16; nstep=200.
- (Final ΔEads uses clean 0.02 SPs on the 3 winners — degauss during relaxation doesn't affect it.)

## Run / resume
```
bash resume.sh        # launch (idempotent, flock-guarded). Re-run on any node after a death.
tail -f run.log       # progress
```
- Driver `run_l442.py` is fully idempotent: derives state from `l442_<tag>.out`
  (converged marker). Partials reseeded from latest geometry. Memory-guarded 2 GPUs.
  Divergence guard: SCF >1 Ry past iter 50 → retry once at degauss=0.06 → else FAILED.
- A durable cron watches and relaunches the driver if it dies.

## Outputs per config
- `l442_<tag>.in` (seed input), `l442_<tag>_run.in` (reseeded, regenerated each run),
  `l442_<tag>.out` (output; overwritten each run — progress is in the reseed chain).
- Markers: `.DIVERGENT04` (retried at 0.06), `.FAILED`, `.retries` (count).

## End state
`python3 analyze_l442.py` when converged:
- ranks configs per variant by relaxed energy → lowest = true minimum;
- writes CIFs (P1, fractional normalized to [0,1)) for ALL converged geometries in `cifs/`;
- generates `sp02/sp_l442_<winner>.in` (clean 0.02 SP) per variant winner; run them, then
  ΔEads = (E_loaded − E_empty(variant) − 4·E_gas)/4, refs in `landscape_3variants/`.
