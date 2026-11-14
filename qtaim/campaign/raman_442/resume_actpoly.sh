#!/bin/bash
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$BASE/act_poly_lS2"
echo "[$(date '+%H:%M:%S')] act_poly_lS2 resume (degauss 0.03, from save) start"
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp rerelax_act_poly_lS2.in > rerelax_act_poly_lS2_resume.out 2>&1
echo "[$(date '+%H:%M:%S')] act_poly_lS2 resume done rc=$? conv=$(grep -c 'bfgs converged' rerelax_act_poly_lS2_resume.out) stop=$(grep -c 'history already reset' rerelax_act_poly_lS2_resume.out) scffail=$(grep -c 'convergence NOT achieved' rerelax_act_poly_lS2_resume.out)"
echo "[$(date '+%H:%M:%S')] RESUME DONE"
