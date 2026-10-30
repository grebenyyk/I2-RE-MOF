#!/bin/bash
# density_442 act pilot: tight real-space SCFs (60/800, degauss 0.01, conv 1e-8, 442 k-mesh)
# on the l442 winner geometries: act_poly_lS2 (lowest polyiodide) + act_mol_rd (molecular,
# speciation contrast). 3 fields each: loaded / empty / i2only. 2 lanes, NP=1. Idempotent.
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
QE=/workspace/qe-7.5/bin/pw.x
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
LOG=$T/run_rs_act.log

exec 9>$T/.rs_act.lock; flock -n 9 || { echo "rs_act driver already running"; exit 0; }

lane() {
  local gpu=$1; shift
  for tag in "$@"; do
    if grep -q '^!' $T/${tag}.out 2>/dev/null; then
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu SKIP $tag" >> $LOG; continue
    fi
    echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu START $tag" >> $LOG
    cd $T && CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE \
      -nk 1 -nd 1 -nb 1 -nt 1 -inp ${tag}.in > ${tag}.out 2>&1
    if grep -q '^!' $T/${tag}.out; then
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu DONE $tag E=$(grep '^!' $T/${tag}.out | tail -1 | awk '{print $5}')" >> $LOG
    else
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu FAIL $tag" >> $LOG
    fi
  done
}

echo "[$(date '+%m-%d %H:%M:%S')] == rs_act batch start (6 jobs, 2 lanes) ==" >> $LOG
lane 0 rs_loaded_act_poly_lS2 rs_empty_act_poly_lS2 rs_i2only_act_poly_lS2 &
P0=$!
lane 1 rs_loaded_act_mol_rd rs_empty_act_mol_rd rs_i2only_act_mol_rd &
P1=$!
wait $P0 $P1
N=$(ls $T/rs_*.out 2>/dev/null | xargs grep -ls '^!' 2>/dev/null | wc -l)
echo "[$(date '+%m-%d %H:%M:%S')] == rs_act batch END: $N/6 with final energy ==" >> $LOG
