#!/bin/bash
export PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:${LD_LIBRARY_PATH:-}
cd /workspace/qe-7.5-libxc
echo "=== build start $(date) ==="
make -j 16 pw 2>&1
echo "=== build exit=$? $(date) ==="
