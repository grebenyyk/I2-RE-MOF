#!/usr/bin/env bash
# usage: run_test.sh <image#> <gpu#>  — standalone SCF of a NEB image at conv_thr=1e-4, mixing 0.1
set -u
img="$1"; gpu="$2"
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=1
CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp test_img${img}.in > test_img${img}.out 2>&1
echo "PW_RC=$? img=$img" >> test_img${img}.out
