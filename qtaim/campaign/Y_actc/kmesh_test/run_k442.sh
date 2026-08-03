#!/bin/bash
set -e
export PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_actc/kmesh_test
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_k442.in > scf_k442.out 2>&1
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/pp.x -inp pp_k442.in > pp_k442.out 2>&1
echo DONE > k442.done
