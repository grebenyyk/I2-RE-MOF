#!/bin/bash
# Finisher: re-run the HSE single-points that didn't complete (driver died in node disconnect).
# 7 tags: 5 never-finished (3 not-started + 2 killed mid-SCF) + 2 prior fails (retry once).
# PBE is already done; this only runs HSE SPs, then regenerates deliverables.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
cd "$D"; mkdir -p _run
TAGS="hse_lpoly_Yvc_S5 hse_lpoly_Yst_S6 hse_lmol_Yst_S4 hse_lmol_Yst_S6 hse_lpoly_Yst_S5 hse_lmol_Yact_S5 hse_lmol_Yst_S5"
for t in $TAGS; do rm -f "_run/$t.done" "_run/$t.fail" "_run/${t}_relaunched.marker"; done
echo "======== HSE FINISHER $(date) ========"
echo "clearing stale markers for: $TAGS"
# split: GPU0 = the 3 polyiodide (sure wins) + 1 st-molecular; GPU1 = remaining molecular incl. 2 fail-retries
G0="hse_lpoly_Yvc_S5 hse_lpoly_Yst_S6 hse_lpoly_Yst_S5 hse_lmol_Yst_S4"
G1="hse_lmol_Yst_S6 hse_lmol_Yact_S5 hse_lmol_Yst_S5"
bash run_hse_conv.sh 0 $G0 > run_finish_gpu0.log 2>&1 &
H0=$!
bash run_hse_conv.sh 1 $G1 > run_finish_gpu1.log 2>&1 &
H1=$!
echo "launched GPU0 (pid $H0): $G0"
echo "launched GPU1 (pid $H1): $G1"
wait $H0; wait $H1
echo "[$(date +%H:%M)] all HSE finisher calcs resolved"
nhse=0; for t in $TAGS; do [ -f "_run/$t.done" ] && nhse=$((nhse+1)); done
echo "finisher HSE done: $nhse/7"
echo "[$(date +%H:%M)] regenerating deliverables + analysis"
python3 build_landscape_deliverables.py
python3 analyze_landscape_full.py
echo "======== FINISHER COMPLETE $(date) ($nhse/7) ========"
