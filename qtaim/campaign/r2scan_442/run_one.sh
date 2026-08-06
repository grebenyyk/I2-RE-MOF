#!/bin/bash
# usage: run_one.sh <input_stem_rho>   e.g. 560  (looks for deriskNC_Yact_r560.in)
rho=$1
PW=/workspace/qe-7.5-libxc/PW/src/pw.x
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/r2scan_442
env -i PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:/usr/bin:/bin \
    HOME=$HOME LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib \
    mpirun --allow-run-as-root -np 32 $PW -inp deriskNC_Yact_r${rho}.in > deriskNC_Yact_r${rho}.out 2>&1
echo "r${rho} exit=$?"
