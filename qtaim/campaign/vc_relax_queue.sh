#!/bin/bash
MPI=/workspace/openmpi; NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
QE=/workspace/qe-7.5/bin/pw.x
BASE=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign
QUEUE=("Gd_str" "Y_actc" "Tb_act")
while [ ${#QUEUE[@]} -gt 0 ]; do
  for gpu in 0 1; do
    used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | sed -n "$((gpu+1))p" | tr -d ' ')
    if [ "$used" -lt 1000 ] && [ ${#QUEUE[@]} -gt 0 ]; then
      sys="${QUEUE[0]}"; QUEUE=("${QUEUE[@]:1}")
      echo "[$(date +%H:%M:%S)] launching $sys on GPU$gpu"
      cd "$BASE/$sys" && rm -rf tmp_vc
      nohup env CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE -nk 1 -nd 1 -nb 1 -nt 1 -inp vc-relax.in > vc-relax.out 2>&1 &
      sleep 10  # let it allocate before checking next
    fi
  done
  sleep 60
done
echo "[$(date +%H:%M:%S)] all vc-relaxes launched"
