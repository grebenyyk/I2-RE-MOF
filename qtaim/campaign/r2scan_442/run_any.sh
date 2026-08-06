#!/bin/bash
# usage: run_any.sh <input_stem> <gpu>   (stem without .in)
stem=$1; gpu=$2
PW=/workspace/qe-7.5-libxc-gpu/PW/src/pw.x
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/r2scan_442
env -i PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/bin:/usr/bin:/bin \
    HOME=$HOME LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/lib64 \
    CUDA_VISIBLE_DEVICES=$gpu \
    mpirun --allow-run-as-root -np 1 $PW -nk 1 -nd 1 -nb 1 -nt 1 -inp ${stem}.in > ${stem}.out 2>&1
echo "${stem} (gpu${gpu}) exit=$?"
