#!/bin/bash
# wait for 4x4x2 SCF to finish, then run ph.x at Gamma
export PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/dfpt
cd "$T"
while true; do
  grep -aq "convergence has been achieved" scf_4x4x2.out 2>/dev/null && break
  sleep 30
done
echo "[$(date +%H:%M:%S)] SCF converged, launching ph.x at Gamma"
sleep 5
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/ph.x -inp ph4.in > ph4.out 2>&1
echo "[$(date +%H:%M:%S)] ph.x finished rc=$?"
echo PHDONE > ph.done
