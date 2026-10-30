#!/bin/bash
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
for f in pp_*.in; do
  o="${f%.in}.out"; c=$(grep fileout $f | awk '{print $3}' | tr -d "'")
  [ -f "$c" ] && { echo "SKIP $c (exists)"; continue; }
  mpirun --allow-run-as-root -np 4 /workspace/qe-7.5/bin/pp.x -inp $f > $o 2>&1 && echo "OK $c" || echo "FAIL $c"
done
echo "PP ALL DONE"
