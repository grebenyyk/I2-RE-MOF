#!/usr/bin/env bash
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
for d in 06Tb_qtaim 08Gd_qtaim; do
  echo "=== $(date '+%H:%M:%S') SCF $d ==="
  cd qtaim/$d
  mpirun --allow-run-as-root -np 2 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf.in > scf.out 2>&1
  echo "  pw.x rc=$? $(grep -m1 'convergence has' scf.out || echo NOCONV)"
  echo "=== $(date '+%H:%M:%S') pp.x $d ==="
  mpirun --allow-run-as-root -np 2 $QE/bin/pp.x -inp pp.in > pp.out 2>&1
  echo "  pp.x rc=$? cube=$(ls -la ae_rho.cube 2>/dev/null | awk '{print $5}')"
  cd ../..
done
echo "=== ALL QTAIM SCF+PP DONE $(date '+%H:%M:%S') ==="
