#!/bin/bash
# Landscape batch-2 driver: 18 new seeds -> PBE relax (plateau-kill) -> HSE SP (poll-kill) -> regenerate.
# 2 GPUs, 2 concurrent. Self-contained; logs to run_landscape2.log.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
cd "$D"; mkdir -p _run
TAGS="lmol_Yact_S4 lmol_Yact_S5 lmol_Yact_S6 lpoly_Yact_S4 lpoly_Yact_S5 lpoly_Yact_S6 \
lmol_Yvc_S4 lmol_Yvc_S5 lmol_Yvc_S6 lpoly_Yvc_S4 lpoly_Yvc_S5 lpoly_Yvc_S6 \
lmol_Yst_S4 lmol_Yst_S5 lmol_Yst_S6 lpoly_Yst_S4 lpoly_Yst_S5 lpoly_Yst_S6"
G0="lmol_Yact_S4 lmol_Yact_S5 lmol_Yact_S6 lpoly_Yact_S4 lpoly_Yact_S5 lpoly_Yact_S6 lmol_Yvc_S4 lmol_Yvc_S5 lmol_Yvc_S6"
G1="lpoly_Yvc_S4 lpoly_Yvc_S5 lpoly_Yvc_S6 lmol_Yst_S4 lmol_Yst_S5 lmol_Yst_S6 lpoly_Yst_S4 lpoly_Yst_S5 lpoly_Yst_S6"
pbe_worker() {
  local gpu=$1; shift
  for tag in "$@"; do
    [ -f "_run/${tag}.done" ] && { echo "[$(date +%H:%M)] PBE skip $tag (done)"; continue; }
    [ ! -f "${tag}.in" ] && { echo "[$(date +%H:%M)] PBE skip $tag (no .in)"; continue; }
    rm -rf "tmp_${tag}"
    echo "[$(date +%H:%M)] PBE START $tag (GPU$gpu)"
    CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "${tag}.in" > "${tag}.out" 2>&1 &
    local pid=$!
    while kill -0 "$pid" 2>/dev/null; do
      sleep 30
      if grep -q 'JOB DONE' "${tag}.out" 2>/dev/null; then
        kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; touch "_run/${tag}.done"
        echo "[$(date +%H:%M)] PBE DONE $tag (JOB DONE)"; break
      fi
      local steps=$(grep -c 'Forces acting on atoms' "${tag}.out" 2>/dev/null)
      local force=$(grep 'Total force' "${tag}.out" 2>/dev/null | tail -1 | awk '{print $4}')
      if [ -n "$force" ] && python3 -c "exit(0 if (float('$force')<0.05 and ${steps:-0}>=8) or (float('$force')<0.2 and ${steps:-0}>=10) else 1)" 2>/dev/null; then
        kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; touch "_run/${tag}.done"
        echo "[$(date +%H:%M)] PBE PLATEAU $tag (force=$force, ${steps} steps)"; break
      fi
    done
    if [ ! -f "_run/${tag}.done" ]; then
      wait "$pid" 2>/dev/null
      grep -q 'JOB DONE' "${tag}.out" 2>/dev/null && touch "_run/${tag}.done" && echo "[$(date +%H:%M)] PBE DONE $tag (exited)"
    fi
  done
}
echo "======== LANDSCAPE BATCH 2 $(date) ========"
echo "[Phase 1] PBE relax (plateau-kill, 2 GPUs)"
pbe_worker 0 $G0 &
P0=$!; pbe_worker 1 $G1 &
P1=$!; wait $P0; wait $P1
nrel=0; for t in $TAGS; do [ -f "_run/${t}.done" ] && nrel=$((nrel+1)); done
echo "[$(date +%H:%M)] relaxed: $nrel/18"
echo "[Phase 2] build HSE inputs"
for t in $TAGS; do
  [ -f "_run/${t}.done" ] || continue
  python3 build_mol_hse_single.py "$t" "hse_${t}" 2>/dev/null && echo "  built hse_${t}.in"
done
echo "[Phase 3] HSE single-points (poll-kill, 2 GPUs)"
HG0=""; HG1=""; i=0
for t in $TAGS; do
  [ -f "_run/${t}.done" ] || continue
  if [ $((i%2)) -eq 0 ]; then HG0="$HG0 hse_${t}"; else HG1="$HG1 hse_${t}"; fi; i=$((i+1))
done
bash run_hse_conv.sh 0 $HG0 &
H0=$!; bash run_hse_conv.sh 1 $HG1 &
H1=$!; wait $H0; wait $H1
nhse=0; for t in $TAGS; do [ -f "_run/hse_${t}.done" ] && nhse=$((nhse+1)); done
echo "[$(date +%H:%M)] HSE done: $nhse/18"
echo "[Phase 4] regenerate deliverables + analysis"
python3 build_landscape_deliverables.py
python3 analyze_landscape_full.py
echo "======== BATCH 2 COMPLETE $(date) ($nhse/18 HSE) ========"
