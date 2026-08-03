#!/bin/bash
# vc-empty fixup: the 3 vc-relaxed empty HSE/Gamma calcs hit a mixing limit cycle
# (energy converged but residual stuck ~1e-3, no '!' line). Inputs already patched
# to mixing_beta=0.05/ndim24/conv_thr1e-5/maxstep300. Re-run on GPU0. Quick (~5min each).
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
GPU=${1:-0}
mkdir -p "$D/_run"; cd "$D"
for tag in hse_empty_yvc hse_empty_gdvc hse_empty_tbvc; do
  [ -e "_run/$tag.done" ] && { echo "[$(date +%H:%M)] $tag done, skip"; continue; }
  rm -rf "tmp_$tag" "_run/$tag.fail"
  echo "[$(date +%H:%M)] START $tag (GPU$GPU)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "$tag.in" > "$tag.out" 2>&1 &
  pid=$!; n=0
  while (( n*15 < 1800 )); do
    sleep 15; n=$((n+1))
    if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
      sleep 15; kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
      touch "_run/$tag.done"
      echo "[$(date +%H:%M)] DONE $tag  E=$(grep '^!.*total energy' "$tag.out"|tail -1|awk '{print $5}')"; break
    fi
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" 2>/dev/null
      if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE(nat) $tag"
      else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (no energy)"; fi
      break
    fi
  done
done
echo "[$(date +%H:%M)] vc-empty fixup COMPLETE"
