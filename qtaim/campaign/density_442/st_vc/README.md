# density_442 st + vc — real-space density at the 4×4×2 tight footing (2026-10-30)

Rollout of the act-pilot pipeline (PBE-D3(BJ), 60/800, degauss 0.01, conv 1e-8, K 4×4×2,
FFT 150×144×320 for vc / 160×180×300 for st, ~0.056 Å) to the st and vc landscape winners.
NOTE: the st loaded SCF does not converge from atomic superposition at 0.01/1e-8 within 300
steps (known metallic-st case) — it was run seeded from the clean-0.02 SP charge density
(`rs_seed_st.in`, converged in 19 iterations, E = −5766.54334492 Ry; the converged density was
swapped into `tmp_rs_loaded_st_poly_lS5/` before cubing).

## Charge transfer, complete 4×4×2 real-space set (Bader-basin δρ, e/cell, 4th-column convention)

| system (config) | C | H | O | Y | I | Σ | bond CPs |
|---|---|---|---|---|---|---|---|
| vc (poly_lS3) | −0.502 | −0.086 | −0.460 | −0.045 | **+1.094** | +0.0000 | 312 |
| st (poly_lS5) | −0.275 | −0.098 | −0.377 | −0.021 | **+0.768** | −0.0036 | 301 |
| act (poly_lS2) | −0.084 | −0.018 | −0.132 | −0.015 | **+0.249** | +0.0000 | 316 |
| act (mol_rd)  | −0.009 | −0.003 | −0.042 | −0.015 | **+0.069** | +0.0005 | 242 |

**Reading (the host-state control story, now on one footing):**
- CT to iodine spans **+0.07 → +1.09 e/cell** across hosts: the more the framework can
  cooperate with the guest (vc = guest-relaxed, st = compressed/anchoring-rich, act = rigid
  open pore), the more charge the iodine sublattice accepts. The linker C/O are the donors in
  every case; Y stays a spectator (−0.015…−0.045).
- act molecular I₂ (+0.069) is the **null baseline** (polarization + weak I···O contacts);
  act polyiodide's +0.249 is 3.6× that; st/vc are 10–16× it. The chain is unambiguously the
  charged species, and its charge scales with the same host cooperation that drives the
  decomposition (vc: assembly −28.8/anchor −109.9; st: anchor −123.6; act: anchor −48.2).
- Consistent with the Raman signature: subvalent chain = partially reduced polyiodide.

## Notes
- vc cubes: pp.x wrote a 2-value line every 54 lines (buffer flush quirk, 21,600 lines);
  all 6,912,000 values present — normalized to 6/line before use (`README` in ../act for the
  act set; `density442_ct_summary.json` here has machine-readable numbers).
- critic2 needs `OMP_NUM_THREADS=1` + `tools/critic2_run.sh` (CRITIC_HOME).
- st loaded SCF required the seeded restart (see above); act set in `../act/`.
