#!/bin/bash
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
cd "$D"; GPU=$1; shift
for tag in "$@"; do
  [ -e "_run/$tag.done" ] && continue
  rm -rf "tmp_$tag"
  echo "[$(date +%H:%M)] START $tag (GPU$GPU)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "$tag.in" > "$tag.out" 2>&1
  if grep -q "JOB DONE" "$tag.out"; then touch "_run/$tag.done"; echo "[$(date +%H:%M)] DONE $tag"; else touch "_run/$tag.fail"; echo "[$(date +%H:%M)] FAIL $tag"; fi
done
