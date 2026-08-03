Attempt 1 (mixing_beta=0.1, conv_thr=1e-6) FAILED identically to the chain calc:
first SCF hit electron_maxstep=300 without converging (accuracy stuck ~2e-6 to 1e-2
oscillation), pw.x printed "convergence NOT achieved after 300 iterations: stopping"
and exited JOB DONE after completing only bfgs step 0. Outputs were inadvertently
overwritten during the regenerate step; see full_relax_chain/attempt1_beta0p1/relax.out
for the same failure mode. Attempt 2 uses mixing_beta=0.2 + conv_thr=1e-5 (proven on
this exact seed geometry by uniform_chain/loaded.out).
