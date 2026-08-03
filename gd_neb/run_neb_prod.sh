#!/usr/bin/env bash
# Production NEB: molecular I2 -> chain on activated Gd framework.
# 2 GPUs, -ni 2 (2 image-groups -> 2 images computed concurrently, 1 GPU each).
set -u
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=0,1
mpirun --allow-run-as-root -np 2 $QE/bin/neb.x -ni 2 -inp neb.in > neb_prod.out 2>&1
echo "NEB_RC=$?" >> neb_prod.out
