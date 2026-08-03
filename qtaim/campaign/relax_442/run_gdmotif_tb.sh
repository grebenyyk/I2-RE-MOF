#!/bin/bash
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
# wait for GPU1 to free (Tbst done)
while [ "$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i 1)" -gt 3000 ]; do sleep 30; done
cd "$T"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
rm -rf tmp_Tbact_gdmotif* 2>/dev/null || true
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE/bin/pw.x" -nk 1 -nd 1 -nb 1 -nt 1 -inp relax442_Tbact_gdmotif.in > relax442_Tbact_gdmotif.out 2>&1
