#!/bin/bash
# Tight density SCFs for the adduct geometry (60/800/442/0.01/1e-8) — loaded on GPU0, i2only on GPU1.
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/adduct"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$B"
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp rs_loaded_adduct.in > rs_loaded_adduct.out 2>&1 &
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp rs_i2only_adduct.in > rs_i2only_adduct.out 2>&1 &
wait
echo "[$(date '+%H:%M:%S')] RS DONE loaded=$(grep -c 'convergence has been achieved' rs_loaded_adduct.out) i2only=$(grep -c 'convergence has been achieved' rs_i2only_adduct.out)"
