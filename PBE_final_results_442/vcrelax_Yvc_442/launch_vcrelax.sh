#!/bin/bash
# vc-relax (cell+ions) of I2@Y-MOFvc at 4x4x2 from the Gamma-converged geometry.
# Seed: rd_poly_Yvc.in structure (Gamma vc-relax result). Protocol: damp/damp,
# cell_dofree=all, degauss 0.04 relax footing (judge at 0.02 SP afterwards).
set -e
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442/vcrelax_Yvc"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$T"
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp vcrelax442_Yvc.in > vcrelax442_Yvc.out 2>&1
