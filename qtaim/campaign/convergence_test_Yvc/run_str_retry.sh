#!/bin/bash
# Standalone retry of the two hard strained-loaded calcs (Y_str, Tb_str) at
# mixing_beta=0.05 / ndim 24 / degauss 0.02, AFTER the main queue finishes.
# Runs them sequentially on GPU0. Launched by the babysitter cron.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
cd "$D"
for v in Y_str Tb_str; do
  [ -e "_run/redo442_loaded_${v}.done" ] && { echo "[$(date)] $v already done, skip" >> _run/redo442_status.log; continue; }
  echo "[$(date)] $v retry start (b0.05/dg0.02)" >> _run/redo442_status.log
  rm -rf "tmp_redo442_loaded_${v}"
  CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "redo442_loaded_${v}.in" > "redo442_loaded_${v}.out" 2>&1
  if grep -q "JOB DONE" "redo442_loaded_${v}.out"; then touch "_run/redo442_loaded_${v}.done"; echo "[$(date)] DONE $v retry" >> _run/redo442_status.log; else touch "_run/redo442_loaded_${v}.fail"; echo "[$(date)] FAIL $v retry" >> _run/redo442_status.log; fi
done
