#!/bin/bash
# Sequential lane: run tight I-only re-relaxes for the listed systems on one GPU.
# usage: run_lane.sh <gpu> <tag1> <tag2>
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
GPU=$1; shift
for TAG in "$@"; do
  cd "$BASE/$TAG"
  echo "[$(date '+%H:%M:%S')] rerelax $TAG start (GPU$GPU)"
  CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
    -inp rerelax_${TAG}.in > rerelax_${TAG}.out 2>&1
  echo "[$(date '+%H:%M:%S')] rerelax $TAG done rc=$? : $(grep -c 'bfgs converged' rerelax_${TAG}.out 2>/dev/null)"
done
echo "[$(date '+%H:%M:%S')] LANE COMPLETE"
