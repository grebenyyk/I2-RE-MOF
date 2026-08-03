#!/bin/bash
# Launch an HSE pw.x SCF on one GPU with the correct NVHPC/OpenMPI env.
# usage: run_hse.sh <gpu> <variant_dir> <inp> <out>
export MPI=/workspace/openmpi
export NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
export CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
GPU=$1; DIR=$2; INP=$3; OUT=$4
cd "$DIR" || exit 1
rm -rf tmp_hse
CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp "$INP" > "$OUT" 2>&1
