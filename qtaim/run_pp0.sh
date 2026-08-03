#!/usr/bin/env bash
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof/qtaim
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
for d in 13Y_opt 13Y_empty 13Y_i2only; do
  pref=$(grep -m1 prefix $d/scf.in | sed "s/.*= *'\(.*\)'.*/\1/")
  ( cd $d && CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 $QE/bin/pp.x -inp pp0.in > pp0.out 2>&1; echo "$d pp0 rc=$? cube=$(stat -c %s val_rho.cube 2>/dev/null || echo 0)" )
done
echo "PP0 DONE"
