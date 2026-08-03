#!/bin/bash
# Resume the 4x4x2 redo, SKIPPING jobs already .done. Used after killing the
# diverging st-loaded calcs and regenerating them at mixing_beta=0.1.
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
RUN=$D/_run; mkdir -p "$RUN"; S=$RUN/redo442_status.log
echo "=== redo442 RESUME start $(date) ===" >> "$S"; echo $$ > "$RUN/redo442.pid"
mapfile -t ALL < <(awk 'NF{print $1}' "$D/redo442_manifest.txt")
JOBS=(); for t in "${ALL[@]}"; do [ -e "$RUN/${t}.done" ] || JOBS+=("$t"); done   # skip done
echo "resuming ${#JOBS[@]} jobs (skipping $((${#ALL[@]}-${#JOBS[@]})) done)" >> "$S"
gpu_free(){ local m; m=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$1" 2>/dev/null|tr -d ' \n'); [ -z "$m" ]&&{ echo 1;return;}; [ "$m" -lt 500 ]&&echo 1||echo 0; }
declare -A pid tag_of; todo=("${JOBS[@]}")
while [ ${#todo[@]} -gt 0 ]||[ $((${#pid[@]})) -gt 0 ]; do
  for g in 0 1; do
    if [ -z "${pid[$g]:-}" ]&&[ ${#todo[@]} -gt 0 ]&&[ "$(gpu_free $g)" = 1 ]; then
      t=${todo[0]}; todo=("${todo[@]:1}"); echo "[$(date +%H:%M:%S)] START $t GPU$g">>"$S"
      ( cd "$D"&&rm -rf "tmp_${t}"&&CUDA_VISIBLE_DEVICES=$g mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp "${t}.in">"${t}.out" 2>&1 )&
      pid[$g]=$!; tag_of[$g]=$t
    fi
  done
  for g in 0 1; do
    if [ -n "${pid[$g]:-}" ]&&! kill -0 "${pid[$g]}" 2>/dev/null; then
      t=${tag_of[$g]}; wait "${pid[$g]}" 2>/dev/null
      grep -q "JOB DONE" "$D/${t}.out" 2>/dev/null && { echo "[$(date +%H:%M:%S)] DONE $t">>"$S"; touch "$RUN/${t}.done"; } || { echo "[$(date +%H:%M:%S)] FAIL $t">>"$S"; touch "$RUN/${t}.fail"; }
      unset 'pid[$g]'; unset 'tag_of[$g]'
    fi
  done
  sleep 15
done
echo "=== redo442 RESUME finished $(date) ===" >> "$S"
