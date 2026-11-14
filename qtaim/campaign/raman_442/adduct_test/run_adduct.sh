#!/bin/bash
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442/adduct_test"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$B"
echo "[$(date '+%H:%M:%S')] adduct I-only relax start GPU1"
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp adduct_relax.in > adduct_relax.out 2>&1
echo "[$(date '+%H:%M:%S')] adduct relax done: conv=$(grep -c 'bfgs converged' adduct_relax.out) fail=$(grep -c 'bfgs failed' adduct_relax.out)"
echo ADDUCT RELAX DONE
