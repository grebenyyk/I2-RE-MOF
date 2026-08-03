#!/usr/bin/env bash
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_vs_hse_442
MPI=/workspace/openmpi; NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$NVHPC/compilers/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64:${LD_LIBRARY_PATH:-}
QE=/workspace/qe-7.5
cd "$D"
# loaded 4x4x2 on GPU0
setsid bash -c "CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp pbe_poly_Yvc.in > pbe_poly_Yvc.out 2>&1" </dev/null >/dev/null 2>&1 &
# empty (Gamma) + gas (Gamma) on GPU1, sequential (fast)
setsid bash -c "CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp pbe_empty_Yvc.in > pbe_empty_Yvc.out 2>&1 && CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp pbe_gas.in > pbe_gas.out 2>&1" </dev/null >/dev/null 2>&1 &
sleep 8; nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader
