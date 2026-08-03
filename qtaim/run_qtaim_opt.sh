#!/usr/bin/env bash
# Tight SCF (60/800/1e-8/Gamma) + pp.x (all-electron cube) for the uniform_chain optimal
# structures, 2 dirs at a time (1 GPU each via CUDA_VISIBLE_DEVICES, NP=1).
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
run_one() {
  local d=$1 gpu=$2
  cd qtaim/$d
  CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf.in > scf.out 2>&1
  local rc=$?
  echo "$(date '+%H:%M:%S') $d pw rc=$rc $(grep -m1 'convergence has' scf.out || echo NOCONV)"
  CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE/bin/pp.x -inp pp.in > pp.out 2>&1
  echo "$(date '+%H:%M:%S') $d pp rc=$? cube=$(stat -c %s ae_rho.cube 2>/dev/null || echo 0)"
  cd ../..
}
export -f run_one; export QE MPI NVHPC CUDA
# GPU0 lane: 13Y_opt then 08Gd_opt ; GPU1 lane: 06Tb_opt
( run_one 13Y_opt 0; run_one 08Gd_opt 0 ) > qtaim/lane0.log 2>&1 &
( run_one 06Tb_opt 1 ) > qtaim/lane1.log 2>&1 &
wait
echo "=== ALL QTAIM SCF+PP DONE $(date '+%H:%M:%S') ==="
