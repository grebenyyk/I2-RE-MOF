#!/usr/bin/env bash
set -u
lbl="$1"; gpu="${2:-0}"
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb/cscan_${lbl}
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=1
CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp relax.in > relax.out 2>&1
echo "RC=$?" >> relax.out
