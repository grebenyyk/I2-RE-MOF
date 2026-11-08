#!/bin/bash
# Resume an interrupted vc-relax from the .save (run on a new node or after a kill).
set -e
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442/vcrelax_Yvc"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd "$T"
sed -e "s/calculation = 'vc-relax'/calculation = 'vc-relax'/" -e '2i\  restart_mode = "restart",' vcrelax442_Yvc.in > vcrelax442_Yvc_resume.in
CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 \
  -inp vcrelax442_Yvc_resume.in > vcrelax442_Yvc.out 2>&1
