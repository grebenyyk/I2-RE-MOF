#!/usr/bin/env bash
# Launch HSE06 60/800/Gamma spot-check calcs (loaded + empty) on the two L40S.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck
MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
QE=/workspace/qe-7.5
cd "$D"
# loaded (nat=80, metallic) on GPU0
setsid bash -c "CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp hse800_vc_loaded.in > loaded.out 2>&1" < /dev/null > /dev/null 2>&1 &
echo "loaded PID=$!"
# empty host (nat=72, insulator) on GPU1
setsid bash -c "CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp hse800_vc_empty.in > empty.out 2>&1" < /dev/null > /dev/null 2>&1 &
echo "empty PID=$!"
sleep 8
echo "=== launched, GPU status ==="
nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader 2>/dev/null
echo "=== loaded.out head (did it start?) ==="
head -5 loaded.out 2>/dev/null
