#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/st_vc
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 /workspace/qe-7.5/bin/pw.x \
  -nk 1 -nd 1 -nb 1 -nt 1 -inp rs_seed_st.in > rs_seed_st.out 2>&1
grep -q '^!' rs_seed_st.out && echo "SEED DONE E=$(grep '^!' rs_seed_st.out | tail -1 | awk '{print $5}')" || echo "SEED FAIL"
