#!/bin/bash
# act pair equilibrium SCFs at degauss 0.04 (from scratch), one per GPU; then FD gen+run.
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
GPU=$1; TAG=$2
cd "$BASE/$TAG"
echo "[$(date '+%H:%M:%S')] scf_loaded $TAG (0.04, scratch) start GPU$GPU"
CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp scf_loaded.in > scf_loaded.out 2>&1
echo "[$(date '+%H:%M:%S')] scf_loaded $TAG done conv=$(grep -c 'convergence has been achieved' scf_loaded.out)"
if grep -q "convergence has been achieved" scf_loaded.out; then
  cd "$BASE/$TAG"
  echo "[$(date '+%H:%M:%S')] FD gen $TAG"
  FD_DEGAUSS=0.04 python3 ../raman_driver_442.py gen $TAG > fd_gen.log 2>&1
  echo "[$(date '+%H:%M:%S')] FD queue $TAG start (GPU$GPU, FD_DEGAUSS=0.04)"
  cd "$BASE"
  FD_DEGAUSS=0.04 RAMAN_GPUS=$GPU nohup python3 raman_driver_442.py run $TAG > fd_run_$TAG.log 2>&1 &
  echo "[$(date '+%H:%M:%S')] FD queue $TAG launched"
fi
