#!/bin/bash
# Rerun the 5 fixed convergence inputs (gasI2 15 A box x3, loaded_6080 PDOS wf_collect-fix x2)
# after the main queue clears. 2-GPU detached, same pattern as run_convergence.sh.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
MPI=/workspace/openmpi
NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
RUN=$D/_run; mkdir -p "$RUN"; STATUS=$RUN/rerun5_status.log
echo "=== rerun5 start $(date) ===" >> "$STATUS"; echo $$ > "$RUN/rerun5.pid"
mapfile -t JOBS < <(awk 'NF>=1 && !/^#/{print $1}' "$D/rerun5.txt")
gpu_free(){ local m; m=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$1" 2>/dev/null | tr -d ' \n'); [ -z "$m" ] && { echo 1; return; }; [ "$m" -lt 500 ] && echo 1 || echo 0; }
declare -A pid tag_of; todo=("${JOBS[@]}")
while [ ${#todo[@]} -gt 0 ] || [ $((${#pid[@]})) -gt 0 ]; do
  for g in 0 1; do
    if [ -z "${pid[$g]:-}" ] && [ ${#todo[@]} -gt 0 ] && [ "$(gpu_free $g)" = 1 ]; then
      t=${todo[0]}; todo=("${todo[@]:1}")
      echo "[$(date +%H:%M:%S)] START $t GPU$g" >> "$STATUS"
      ( cd "$D" && rm -rf "tmp_${t}" && CUDA_VISIBLE_DEVICES=$g mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "${t}.in" > "${t}.out" 2>&1 ) &
      pid[$g]=$!; tag_of[$g]=$t
    fi
  done
  for g in 0 1; do
    if [ -n "${pid[$g]:-}" ] && ! kill -0 "${pid[$g]}" 2>/dev/null; then
      t=${tag_of[$g]}; wait "${pid[$g]}" 2>/dev/null
      if grep -q "JOB DONE" "$D/${t}.out" 2>/dev/null; then echo "[$(date +%H:%M:%S)] DONE $t" >> "$STATUS"; touch "$RUN/$t.done"
      else echo "[$(date +%H:%M:%S)] FAIL $t" >> "$STATUS"; touch "$RUN/$t.fail"; fi
      unset 'pid[$g]'; unset 'tag_of[$g]'
    fi
  done
  sleep 15
done
echo "=== rerun5 finished $(date) ===" >> "$STATUS"
