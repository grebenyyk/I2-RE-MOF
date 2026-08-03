#!/bin/bash
# HSE06 for Yvc + Yst (metallic variants). Empties are Gamma/insulating (fast);
# loaded are 4x4x2 (expensive metallic). Poll-for-energy-then-kill (metallic HSE
# prints energy then may hang at PAW-EXX). Cap 6h per loaded calc. GPU0.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
GPU=${1:-0}
mkdir -p "$D/_run"; cd "$D"
run_one () {
  local tag="$1" cap="$2"
  [ -e "_run/$tag.done" ] && { echo "[$(date +%H:%M)] $tag done, skip"; return; }
  rm -rf "tmp_$tag"
  echo "[$(date +%H:%M)] START $tag (GPU$GPU, cap ${cap}s)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "$tag.in" > "$tag.out" 2>&1 &
  local pid=$! n=0
  while (( n*30 < cap )); do
    sleep 30; n=$((n+1))
    if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
      sleep 30; kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
      touch "_run/$tag.done"
      echo "[$(date +%H:%M)] DONE $tag  E=$(grep '^!.*total energy' "$tag.out"|tail -1|awk '{print $5}')"
      return
    fi
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" 2>/dev/null
      if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE(nat) $tag"
      else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (no energy)"; fi
      return
    fi
  done
  kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
  if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE(cap) $tag"
  else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (${cap}s cap, no energy)"; fi
}
# empties (Gamma, fast) first, then loaded (fine mesh, expensive)
run_one hse_empty_yvc  3600
run_one hse_empty_yst  3600
run_one hse_poly_yvc  21600   # 6h cap (metallic HSE 4x4x2)
run_one hse_poly_yst  21600
echo "[$(date +%H:%M)] HSE Yvc+Yst queue COMPLETE"
