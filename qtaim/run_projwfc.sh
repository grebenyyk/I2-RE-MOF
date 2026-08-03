#!/usr/bin/env bash
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
QE=/workspace/qe-7.5; MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
for d in 13Y_rho800_G 06Tb_qtaim 08Gd_qtaim; do
  echo "=== $(date '+%H:%M:%S') SCF(wf_collect) $d ==="
  cd qtaim/$d
  mpirun --allow-run-as-root -np 2 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf.in > scf_wc.out 2>&1
  echo "  pw rc=$? $(grep -m1 'convergence has' scf_wc.out || echo NOCONV)"
  cat > proj.in <<E2
&projwfc
  prefix = '$d'
  outdir = './tmp'
  DeltaE = 0.01
  filpdos = 'pdos'
/
E2
  echo "=== projwfc $d ==="
  mpirun --allow-run-as-root -np 2 $QE/bin/projwfc.x -inp proj.in > proj.out 2>&1
  echo "  projwfc rc=$? pdosfiles=$(ls pdos_atm*.pdos_atm 2>/dev/null | wc -l)"
  cd ../..
done
echo "=== ALL PROJWFC DONE $(date '+%H:%M:%S') ==="
