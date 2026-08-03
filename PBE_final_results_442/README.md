# PBE_final_results_442 — 4×4×2-relaxed campaign results

Mirrors `PBE_final_results/` (the Γ-relaxed-geometry production set) but with geometries
**relaxed at the 4×4×2 k-mesh**. The Γ-relaxed geometries undersample the PES by 6–26 kJ/I₂
(see `delta` columns); this directory holds the corrected, properly-relaxed results.

## Layout (mirrors PBE_final_results)
```
adsorption_9systems/   # 9 systems (Y/Gd/Tb × vc/act/st), 4×4×2-relaxed loaded + Γ empty/gas
  relax_442_<sys>.in/out  # the ACTUAL relaxation (calculation='relax', damp) that produced the geometry
  rd_poly_442_<sys>.in/out# clean single-point (calculation='scf', degauss=0.02, 4×4×2) at that relaxed
                          #   geometry — the energy used for ΔEads. (Mirrors original rd_poly_*, which
                          #   were also scf single-points on the Γ-relaxed geometry.)
  rd_empty_<sys>.out      # empty framework (Γ, 0.02) — reused from PBE_final_results (unchanged)
  rd_gas.out              # gas I2 (Γ, 0.02) — reused
  compute_442_table.py    # self-refreshing ΔEads table (pulls latest SP02, writes adsorption_table_442.csv)
  adsorption_table_442.csv
landscape_3variants/    # (pending) 18-config Y landscape, 4×4×2-relaxed — populated when the campaign finishes
```

### Why two inputs per system (relax + scf)?
The relaxation (`relax_442_*`, calculation='relax') is what we actually ran to find the
4×4×2-relaxed geometry. The act/st systems were relaxed at degauss=0.04 (cold-start SCF aid),
so their relaxation energy isn't at the production smearing. `rd_poly_442_*` (calculation='scf',
degauss=0.02) re-evaluates the energy at the relaxed geometry with the production smearing —
that clean energy is what enters ΔEads. (vc and Yact were relaxed at 0.02 already, so for them
the relaxation output and the scf single-point give the same energy.) This two-file split matches
the original `PBE_final_results` convention where `rd_poly_*` was the scf single-point and the
Γ-relaxation lived separately.

## Recipe (relaxation, in `qtaim/campaign/relax_442/`)
I-only relax (framework frozen via if_pos), 4×4×2, ecutwfc=50/rho=400, conv_thr=1e-6,
smearing=mv. Optimizer: bfgs for vc, damp for act/st (soft polyiodide modes). degauss=0.04
during the act/st relaxation (cold-start SCF aid); **final energies here are clean 0.02 SPs**
(SP02) on those relaxed geometries, so the reported ΔEads is at production degauss=0.02.
Empty frameworks + gas I2 are Γ-converged (insulating; validated) — reused unchanged.

## Status
- **adsorption_9systems:** 9 geometries relaxed (Tbst accepted at energy-converged floor F=0.013).
  Clean 0.02 energies: vc×3, Yact, Gdact done; Tbact/Yst/Gdst/Tbst = SP02 in progress → PENDING.
- **To refresh as SP02 completes:** `python3 adsorption_9systems/compute_442_table.py`
  (auto-copies converged SP02 outputs, recomputes the table).
- **landscape_3variants:** populated from `qtaim/campaign/landscape_442/` when the 18-config
  campaign finishes (its own durable cron + `resume.sh`).

## Key result (so far)
4×4×2 relaxation **deepens I₂ binding** vs Γ: vc +6, Yact +15, Gdact +6 kJ/I₂ (clean 0.02).
The Γ-relaxed production geometries underestimate adsorption, most severely in the open act/st pores.
