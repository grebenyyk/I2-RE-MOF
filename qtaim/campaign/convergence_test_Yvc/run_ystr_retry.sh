#!/bin/bash
# Standalone retry of redo442_loaded_Y_str (the most-metallic system) after the
# main queue finishes, at mixing_beta=0.05 / ndim 24 / degauss 0.02. Launched by
# the babysitter cron when the queue is done and Y_str is not yet .done.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
cd "$D" && rm -rf tmp_redo442_loaded_Y_str
echo "[$(date)] Y_str retry start" >> _run/redo442_status.log
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp redo442_loaded_Y_str.in > redo442_loaded_Y_str.out 2>&1
if grep -q "JOB DONE" redo442_loaded_Y_str.out; then touch _run/redo442_loaded_Y_str.done; echo "[$(date)] DONE Y_str retry" >> _run/redo442_status.log; else touch _run/redo442_loaded_Y_str.fail; echo "[$(date)] FAIL Y_str retry" >> _run/redo442_status.log; fi
