#!/bin/bash
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
for f in pp_rho17_act_poly_lS2.in pp_val_act_poly_lS2.in pp_empty_val_act_poly_lS2.in pp_i2only_val_act_poly_lS2.in; do
  mpirun --allow-run-as-root -np 4 /workspace/qe-7.5/bin/pp.x -inp $f > ${f%.in}.out 2>&1 && echo "OK $f" || echo "FAIL $f"
done
echo "PPFIX DONE"
