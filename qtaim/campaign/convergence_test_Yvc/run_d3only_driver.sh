#!/bin/bash
# Build + run PBE+D3 Gamma single-points for all landscape configs + references,
# then decompose ΔEads = ΔEdisp + ΔEelec. Detached (survives disconnect).
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
cd "$D"
echo "[$(date +%H:%M)] building d3only inputs..."; python3 build_d3only.py
TAGS=$(ls d3only_*.in | sed 's,d3only_,,; s,\.in,,')
G0=""; G1=""; i=0
for t in $TAGS; do
  if [ $((i%2)) -eq 0 ]; then G0="$G0 $t"; else G1="$G1 $t"; fi; i=$((i+1))
done
worker() {
  local gpu=$1; shift
  for tag in "$@"; do
    local pre="d3only_$tag"
    [ -f "_run/${pre}.done" ] && continue
    rm -rf "tmp_$pre"
    CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "${pre}.in" > "${pre}.out" 2>&1 &
    local pid=$!; local n=0
    while kill -0 "$pid" 2>/dev/null; do
      sleep 20; n=$((n+1))
      if grep -qE 'DFT-D3 Dispersion +=' "${pre}.out" 2>/dev/null && grep -q 'JOB DONE\|^!.*total energy' "${pre}.out" 2>/dev/null; then
        kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; touch "_run/${pre}.done"; break
      fi
      if (( n*20 > 600 )); then kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
        grep -qE 'DFT-D3 Dispersion +=' "${pre}.out" 2>/dev/null && touch "_run/${pre}.done"
        break
      fi
    done
    [ ! -f "_run/${pre}.done" ] && { wait "$pid" 2>/dev/null; grep -qE 'DFT-D3 Dispersion +=' "${pre}.out" 2>/dev/null && touch "_run/${pre}.done"; }
    echo "[$(date +%H:%M)] $pre: $([ -f _run/${pre}.done ] && echo OK || echo MISS)"
  done
}
mkdir -p _run
echo "[$(date +%H:%M)] running d3only SPs on 2 GPUs ($i tags)"
worker 0 $G0 & W0=$!; worker 1 $G1 & W1=$!; wait $W0; wait $W1
nok=0; for t in $TAGS; do [ -f "_run/d3only_$t.done" ] && nok=$((nok+1)); done
echo "[$(date +%H:%M)] d3only done: $nok/$i"
echo "[$(date +%H:%M)] decomposing..."; python3 analyze_d3_decomp.py
echo "[$(date +%H:%M)] D3 DECOMPOSITION COMPLETE"
