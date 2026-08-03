#!/bin/bash
# Phase B + C-light HSE06/Gamma queue. Sequential on GPU0 (GPU1 busy w/ poly relax).
# KEY: HSE on the metallic loaded cells prints '!total energy' then HANGS at the
# PAW-EXX post-step. So we poll for the energy line, capture it, then kill pw.x to
# free the GPU and advance. Insulators (empty, gas-I2) exit normally on their own.
# Success = '!    total energy' present in .out (JOB DONE is NOT required/expected).
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
GPU=${1:-0}
mkdir -p "$D/_run"
cd "$D"
MAXWAIT=240   # polls x 15s = 60 min cap per calc (HSE SCF converges <24min)
run_one () {
  local tag="$1"
  [ -e "_run/$tag.done" ] && { echo "[$(date +%H:%M)] $tag already done, skip"; return; }
  rm -rf "tmp_$tag"
  echo "[$(date +%H:%M)] START $tag (GPU$GPU)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "$tag.in" > "$tag.out" 2>&1 &
  local pid=$!
  local n=0
  while (( n < MAXWAIT )); do
    sleep 15; n=$((n+1))
    if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
      sleep 20   # let an insulator finish naturally; energy already captured either way
      kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
      touch "_run/$tag.done"
      echo "[$(date +%H:%M)] DONE $tag  E=$(grep '^!.*total energy' "$tag.out" | tail -1 | awk '{print $NF}')"
      return
    fi
    if ! kill -0 "$pid" 2>/dev/null; then      # process exited on its own (insulator)
      wait "$pid" 2>/dev/null
      if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
        touch "_run/$tag.done"
        echo "[$(date +%H:%M)] DONE(natural) $tag  E=$(grep '^!.*total energy' "$tag.out" | tail -1 | awk '{print $NF}')"
      else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (exited, no energy)"; fi
      return
    fi
  done
  kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null   # 60min cap
  if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
    touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE(cap) $tag"
  else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (60min cap, no energy)"; fi
}
# references first (insulating, fast -> table assembles early), then loaded
for tag in hse_i2_vac hse_empty_yact hse_empty_gdact hse_empty_tbact \
           hse_mol_gdact hse_mol_tbact hse_poly_gdact hse_poly_tbact; do
  run_one "$tag"
done
echo "[$(date +%H:%M)] HSE B+C-light queue COMPLETE"
