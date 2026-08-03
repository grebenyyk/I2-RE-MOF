#!/usr/bin/env bash
# Run SCF (conv_thr 1e-8) + pp.x (plot_num=21 all-electron density cube) for each conv dir.
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
PW=$QE/bin/pw.x; PP=$QE/bin/pp.x
DIRS="13Y_rho400_G 13Y_rho600_G 13Y_rho800_G 13Y_rho600_2x2x2"
for d in $DIRS; do
  echo "=== $(date '+%H:%M:%S') SCF $d ==="
  cd qtaim/$d
  NP=2 mpirun --allow-run-as-root -np 2 $PW -nk 1 -nd 1 -nb 1 -nt 1 -inp scf.in > scf.out 2>&1
  rc=$?
  echo "  pw.x rc=$rc  $(grep -m1 'convergence has' scf.out || echo '(no conv line)')"
  if [ $rc -eq 0 ]; then
    mpirun --allow-run-as-root -np 1 $PP -inp pp.in > pp.out 2>&1
    echo "  pp.x rc=$?  cube: $(ls -la ae_rho.cube 2>/dev/null | awk '{print $5}') bytes"
  fi
  cd ../..
done
echo "=== $(date '+%H:%M:%S') ALL CONV SCF+PP DONE ==="
