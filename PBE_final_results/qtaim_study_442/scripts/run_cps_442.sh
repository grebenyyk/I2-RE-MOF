#!/usr/bin/env bash
# Build the four 442 CPs.txt sets (bcp_raw + refine merged) — run from mwfn_out.
set -euo pipefail
CAMP=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442
S=/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/scripts
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/mwfn_out

# st: chain wraps into +a1 -> 3x positive a1, negative a2/a3
python3 $S/mwfn_cps_442.py $CAMP/st_vc/rho17_st_poly_lS5.cube CPs_st.txt \
  $CAMP/st_vc/bcp_raw_st_poly_lS5.out $CAMP/st_vc/refine_st_poly_lS5.out \
  --shift0 0,1,2 --shift1 0,-1 --shift2 0,-1
# vc / act / actmol: aggregates in-cell -> 2x negative each axis
for s in "vc st_vc vc_poly_lS3" "act act act_poly_lS2" "actmol act act_mol_rd"; do
  set -- $s
  python3 $S/mwfn_cps_442.py $CAMP/$2/rho17_$3.cube CPs_$1.txt \
    $CAMP/$2/bcp_raw_$3.out $CAMP/$2/refine_$3.out \
    --shift0 0,-1 --shift1 0,-1 --shift2 0,-1
done
