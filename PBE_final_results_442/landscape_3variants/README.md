# landscape_3variants (4×4×2-relaxed) — COMPLETE

18 Y-config landscape from `qtaim/campaign/landscape_442/` (all 18 converged 2026-08-05,
clean degauss-0.02 SPs added 2026-10-30). Mirrors `PBE_final_results/landscape_3variants/`
but with 4×4×2-relaxed geometries: 18 config .in/.out, per-variant rankings + winner
geometries/CIFs, `motif_preference_table_442.csv` with clean 0.02 ΔEads for the SP'd configs.

## Clean 0.02 result (production footing)

| variant | best config | ΔEads (kJ/I₂) | motif | runner-up motif gap |
|---|---|---|---|---|
| st | poly_lS5 | **−137.14** | polyiodide | mol +20.8 |
| act | mol_rd | **−67.30** | molecular I₂ | poly +7.2 |
| vc | poly_lS3 | **−99.06** | polyiodide | mol ≥ +37 (relax-footing) |

Sanity: landscape winners reproduce the 9-system adsorption table (vc −99.06 vs −99.26,
st −137.14 vs −136.53; act SP of the production geometry −67.12 sits in the molecular basin).

**The act flip is real and explained.** On relaxation-run energies (degauss 0.04) act
ranked polyiodide first by 1.1 kJ; at clean 0.02 the ranking reverses: molecular wins by
7.2 kJ/I₂. Cause: the wide 0.04 relaxation smearing penalizes molecular configs (partially
populates framework/conduction states) far more than the intrinsically metallic polyiodide —
a biased footing, not physics. Fingerprint check: act_mol_rd = 4× isolated I–I 2.70 Å
(molecular); act_poly_lS2 = cores 2.80 + true 3.13–3.16 Å links (associated chain); the
production rd_poly_442_Yact geometry (2.73–2.77 + 3.05–3.37 Å links) SPs in the molecular
basin — quasi-molecular with stretched I₂ units.

**Implication for the paper:** at production footing DFT favors molecular I₂ in the activated
pore by ~7 kJ/I₂ while the experiment (Raman) sees polyiodide. Both motifs are within ~7 kJ;
the strained (st) pore prefers polyiodide by 21 kJ, so partial strain/solvation retention or
high-loading cooperative effects can tip the balance. State this openly; the Raman experiment
is the arbiter of the realized motif (see also the act_poly_lS5 confirmation SP, planned).

## Files
- `motif_preference_table_442.csv` — all 18 configs: relax energies (0.04) + clean-0.02
  energies, ΔE_vs_min and ΔEads for the 5 SP'd configs
- `l442_<tag>.{in,out}` + `sp_l442_<tag>.{in,out}` — relaxations and clean SPs
- `cifs/l442_<tag>_relaxed.cif` (18) + `_LOWEST` copies of the 3 relax-stage winners
- master campaign dir: `qtaim/campaign/landscape_442/` (driver `run_l442.py`, `resume.sh`,
  `analyze_l442.py`, sp02/)

## Caveats
- Remaining near-degenerate pairs unconfirmed at clean footing (relax-footing gaps ≤ 1.7 kJ,
  i.e. within the 0.04→0.02 footing shift): st poly_S1 (0.0), act poly_lS5 (1.7), vc poly_lS6
  (0.4), st poly_lS4 (0.7). Confirmation SPs queued in the follow-up campaign.
- Landscape ΔEads uses 4-I₂-per-cell orderings only; real saturation loadings may differ.
