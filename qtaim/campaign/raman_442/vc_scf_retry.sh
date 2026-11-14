#!/bin/bash
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$BASE/vc_poly_lS3"
echo "[$(date '+%H:%M:%S')] vc scf_loaded retry start"
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp scf_loaded.in > scf_loaded.out 2>&1
echo "[$(date '+%H:%M:%S')] vc scf_loaded retry done conv=$(grep -c 'convergence has been achieved' scf_loaded.out)"
