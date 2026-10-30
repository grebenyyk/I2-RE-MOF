#!/bin/bash
# decomp_442 batch: 6 frozen-fragment SCFs (isubl + hostld for the 3 variant winners,
# PBE-D3 442/0.02 — the footing of the clean SP02s) + 4 winner-confirmation SPs.
# Two lanes, one L40S each, NP=1. Idempotent (skips outs with final energy). flock-guarded.
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
QE=/workspace/qe-7.5/bin/pw.x
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/decomp_442
LOG=$T/run_dc442.log

exec 9>$T/.dc442.lock; flock -n 9 || { echo "dc442 driver already running"; exit 0; }

lane() {  # lane <gpu> <infile-basename>...
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

echo "[$(date '+%m-%d %H:%M:%S')] == dc442 batch start (10 jobs, 2 lanes) ==" >> $LOG
lane 0 sp_l442_act_poly_lS5 sp_l442_st_poly_S1 sp_l442_vc_poly_lS6 sp_l442_st_poly_lS4 &
P0=$!
lane 1 dc442_isubl_st_poly_lS5 dc442_hostld_st_poly_lS5 dc442_isubl_act_poly_lS2 dc442_hostld_act_poly_lS2 dc442_isubl_vc_poly_lS3 dc442_hostld_vc_poly_lS3 &
P1=$!
wait $P0 $P1
N=$(ls $T/*.out 2>/dev/null | xargs grep -ls '^!' 2>/dev/null | wc -l)
echo "[$(date '+%m-%d %H:%M:%S')] == dc442 batch END: $N/10 with final energy ==" >> $LOG
