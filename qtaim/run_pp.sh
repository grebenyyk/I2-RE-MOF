#!/usr/bin/env bash
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
for d in 13Y_rho400_G 13Y_rho600_G 13Y_rho800_G; do
  echo "=== $(date '+%H:%M:%S') pp.x $d ==="
  cd qtaim/$d
  mpirun --allow-run-as-root -np 2 $QE/bin/pp.x -inp pp.in > pp.out 2>&1
  echo "  rc=$? cube=$(ls -la ae_rho.cube 2>/dev/null | awk '{print $5}')"
  cd ../..
done
echo "=== ALL PP DONE $(date '+%H:%M:%S') ==="
