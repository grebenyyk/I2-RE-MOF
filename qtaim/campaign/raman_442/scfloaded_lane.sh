#!/bin/bash
# Short SCFs (equilibrium geometry from re-relax save), shared on GPU1.
set -uo pipefail
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
for TAG in vc_poly_lS3 st_poly_lS5; do
  cd "$BASE/$TAG"
  echo "[$(date '+%H:%M:%S')] scf_loaded $TAG start"
  CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
    -inp scf_loaded.in > scf_loaded.out 2>&1
  echo "[$(date '+%H:%M:%S')] scf_loaded $TAG done conv=$(grep -c 'convergence has been achieved' scf_loaded.out)"
done
echo "[$(date '+%H:%M:%S')] SCFLOADED COMPLETE"
