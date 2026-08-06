#!/bin/bash
cd /workspace/qe-7.5-libxc-gpu
echo "=== GPU build start $(date) ==="
env -i PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/bin:/usr/bin:/bin \
    HOME=$HOME LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/lib64 \
    make -j 16 pw 2>&1
echo "=== GPU build exit=$? $(date) ==="
