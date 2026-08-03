#!/bin/bash
# Run projwfc.x (CPU) to produce PDOS from a converged SCF save dir.
# usage: run_proj.sh <variant_dir> <inp> <out>
export MPI=/workspace/openmpi
export NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
export CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
DIR=$1; INP=$2; OUT=$3
cd "$DIR" || exit 1
/workspace/qe-7.5/bin/projwfc.x -inp "$INP" > "$OUT" 2>&1
