#!/bin/bash
# SINGLE all-in-one DFPT driver for Y_deep 4x4x2: SCF -> ph.x(Gamma) -> done. Run ONCE.
export PATH=/workspace/openmpi/bin:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/dfpt
cd "$T"
echo "[$(date +%H:%M:%S)] SCF start"
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_4x4x2.in > scf_4x4x2.out 2>&1
if ! grep -aq "convergence has been achieved" scf_4x4x2.out; then echo "[$(date +%H:%M:%S)] SCF FAILED"; echo SCFFAIL>ph.done; exit 1; fi
echo "[$(date +%H:%M:%S)] SCF converged -> ph.x(Gamma) start"
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/ph.x -inp ph4.in > ph4.out 2>&1
echo "[$(date +%H:%M:%S)] ph.x done rc=$?"
echo PHDONE > ph.done
