#!/bin/bash
set -e
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/vc_relax_442
rm -rf tmp_vcr_$2* 2>/dev/null || true
CUDA_VISIBLE_DEVICES=$1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 -inp vcr_$2.in > vcr_$2.out 2>&1
