#!/bin/bash
# Y-vc HSE convergence sweep runner. usage: run_hse_conv.sh <GPU> <tag1> <tag2> ...
# Poll-for-energy-then-kill (metallic HSE prints energy then may hang at PAW-EXX).
# Cap 90 min/calc (4x4x2/nqx4 Fock is heavy). Success = '! total energy' present.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
GPU=$1; shift
mkdir -p "$D/_run"; cd "$D"
for tag in "$@"; do
  [ -e "_run/$tag.done" ] && { echo "[$(date +%H:%M)] $tag done, skip"; continue; }
  rm -rf "tmp_$tag"
  echo "[$(date +%H:%M)] START $tag (GPU$GPU)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "$tag.in" > "$tag.out" 2>&1 &
  pid=$!; n=0; ok=0
  while (( n*30 < 5400 )); do
    sleep 30; n=$((n+1))
    if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then
      sleep 30; kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; touch "_run/$tag.done"; ok=1
      echo "[$(date +%H:%M)] DONE $tag  E=$(grep '^!.*total energy' "$tag.out"|tail -1|awk '{print $5}')"; break
    fi
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" 2>/dev/null
      if grep -q '^!.*total energy' "$tag.out" 2>/dev/null; then touch "_run/$tag.done"; ok=1; echo "[$(date +%H:%M)] DONE(nat) $tag"
      else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (no energy)"; fi
      break
    fi
  done
  if (( ! ok )); then kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
    grep -q '^!.*total energy' "$tag.out" 2>/dev/null && { touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE(cap) $tag"; } || { touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag (cap)"; }
  fi
done
echo "[$(date +%H:%M)] GPU$GPU sweep batch done"
