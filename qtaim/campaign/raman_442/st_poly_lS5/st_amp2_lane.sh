#!/bin/bash
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442/st_poly_lS5"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$B"
for T in I073_xp I073_xn; do
  mkdir -p tmp_$T
  CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
    -inp raman_amp_check/disp2_$T.in > raman_amp_check/$T.out 2>&1
  echo "[$T] conv=$(grep -c 'convergence has been achieved' raman_amp_check/$T.out)"
done
echo AMP2 DONE
