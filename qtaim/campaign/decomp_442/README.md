# decomp_442 — frozen-fragment decomposition at the 4×4×2 clean footing (2026-10-30)

Rebuild of the Γ-era decomposition table (`HSE_final_tables/decomposition_Y3`, which is
PBE-electronic due to the QE 7.5 PAW-EXX issue) for the **lowest polyiodide of each variant
at its 4×4×2-relaxed geometry**, everything at production footing: PBE-D3(BJ), degauss 0.02,
loaded SPs and fragments at the 4×4×2 mesh, empty/gas references the production Γ SCFs
(k-exact; validated). Same convention as the old `analyze_decomp.py`: host⁰ = per-variant
`rd_empty_Y<var>`; only vc carries a real ΔE_def (its frozen guest-templated framework vs
that reference); st/act ΔE_def ≈ 0 (small Γ-vs-442 empty footing residual).

ΔE_assembly = (E_I* − 4·E_gas)/4;  ΔE_anchor = (E_loaded − E_host* − E_I*)/4;  ΔE_def = (E_host* − E_host⁰)/4

## Table (kJ/mol per I₂)

| variant | config | ΔE_assembly | ΔE_anchor | ΔE_def | Σ | ΔEads (SP02) |
|---|---|---|---|---|---|---|
| st | st_poly_lS5 | −11.58 | −123.59 | −1.97 | −137.14 | −137.14 ✓ |
| act | act_poly_lS2 | −10.61 | −48.20 | −1.24 | −60.05 | −60.05 ✓ |
| vc | vc_poly_lS3 | −28.76 | −109.92 | +39.63 | −99.06 | −99.06 ✓ |

Identity closes exactly in all three rows. vs the Γ-era table the structure is unchanged
(anchor dominates st/vc; act weak) with vc def +39.63 vs old +42.77.

**Reading:**
- **act anchors weakly** (−48.2): consistent with the clean-footing finding that act's true
  minimum is molecular I₂ (−67.3); the associated chain in act is ~7 kJ/I₂ above it.
- **vc stabilizes the chain cooperatively**: assembly −28.8 (vs ~−11 in st/act) *and*
  anchoring −109.9 — the guest-equilibrated framework lets both the sublattice and the
  chain–framework contact relax. The +39.6 def is what the framework pays for it, netting
  −99.1.
- **st = anchoring-dominated** (−123.6, 90% of ΔEads): the strained host's carboxylates
  anchor the chain at essentially no deformation cost (def ≈ −2).

## Winner confirmation SPs (clean 0.02; no flips anywhere)

| variant | ranking (kJ/I₂ above winner) |
|---|---|
| st | **poly_lS5 MIN**, poly_S1 +0.26, poly_lS4 +0.70 |
| act (polyiodide track) | **poly_lS2 MIN**, poly_lS5 +0.70; overall act minimum remains act_mol_rd (molecular, −67.30) |
| vc | **poly_lS3 MIN**, poly_lS6 +0.41 |

All relax-footing winners survive the clean footing; the `_LOWEST` stamps for st/vc are
confirmed, act's polyiodide `_LOWEST` is the lowest *polyiodide* (the overall act lowest is
molecular act_mol_rd).

## Files
`dc442_{isubl,hostld}_<tag>.{in,out}` — fragment SCFs (8-I sublattice / 72-atom frozen host);
`sp_l442_<tag>.{in,out}` — confirmation SPs; `run_dc442.sh`, `run_dc442.log` — driver.
Loaded energies from `qtaim/campaign/landscape_442/sp02/`; references from
`PBE_final_results_442/adsorption_9systems/`.
