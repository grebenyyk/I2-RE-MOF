#!/usr/bin/env bash
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=1
NP="${NP:-2}"; NI="${NI:-2}"
mpirun --allow-run-as-root -np $NP $QE/bin/neb.x -ni $NI -inp neb.in > neb.out 2>&1
echo "neb rc=$?"
