#!/bin/bash
# Queue the 9 empty-framework vc-relax runs on 2 L40S GPUs (1 calc/GPU, NP=1).
# Detached: launch with  nohup setsid bash run_empty_vc.sh > _empty_vc_run/queue.out 2>&1 &
# Polling/restartable: tracks PIDs + GPU mem; writes _empty_vc_run/status.log;
# touches <variant>.done / <variant>.fail on completion.
set -u
BASE=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign
QE=/workspace/qe-7.5/bin/pw.x
MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2

JOBS=(Y_str Y_actc Y_deep Gd_str Gd_act Gd_deep Tb_str Tb_act Tb_deep)
RUN=$BASE/_empty_vc_run
mkdir -p "$RUN"
STATUS=$RUN/status.log
echo "=== queue start $(date)  host=$(hostname) ===" >> "$STATUS"
echo $$ > "$RUN/queue.pid"

gpu_free(){ # $1=gpu -> echo 1 if memory.used < 500 MiB
  local m
  m=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$1" 2>/dev/null | tr -d ' \n')
  [ -z "$m" ] && { echo 1; return; }
  [ "$m" -lt 500 ] && echo 1 || echo 0
}

declare -A pid gpu_of
todo=("${JOBS[@]}")
while [ ${#todo[@]} -gt 0 ] || [ $((${#pid[@]})) -gt 0 ]; do
  # launch on free GPUs
  for g in 0 1; do
    if [ -z "${pid[$g]:-}" ] && [ ${#todo[@]} -gt 0 ] && [ "$(gpu_free $g)" = 1 ]; then
      v=${todo[0]}; todo=("${todo[@]:1}")
      echo "[$(date +%H:%M:%S)] START $v on GPU$g" >> "$STATUS"
      ( cd "$BASE/$v" && rm -rf tmp_empty_vc && \
        CUDA_VISIBLE_DEVICES=$g mpirun --allow-run-as-root -np 1 "$QE" \
          -nk 1 -nd 1 -nb 1 -nt 1 -inp vc_empty.in > vc_empty.out 2>&1 ) &
      pid[$g]=$!; gpu_of[$g]=$v
    fi
  done
  # reap finished
  for g in 0 1; do
    if [ -n "${pid[$g]:-}" ] && ! kill -0 "${pid[$g]}" 2>/dev/null; then
      v=${gpu_of[$g]}; wait "${pid[$g]}" 2>/dev/null
      if grep -q "JOB DONE" "$BASE/$v/vc_empty.out" 2>/dev/null; then
        echo "[$(date +%H:%M:%S)] DONE  $v (GPU$g)" >> "$STATUS"; touch "$RUN/$v.done"
      else
        echo "[$(date +%H:%M:%S)] FAIL  $v (GPU$g) -- see $v/vc_empty.out" >> "$STATUS"; touch "$RUN/$v.fail"
      fi
      unset 'pid[$g]'; unset 'gpu_of[$g]'
    fi
  done
  sleep 15
done
echo "=== queue finished $(date) ===" >> "$STATUS"
