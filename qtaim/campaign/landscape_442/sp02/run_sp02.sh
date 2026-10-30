#!/bin/bash
# Clean degauss=0.02 SCF single-points for the landscape_442 winners (+motif rivals).
# Two lanes, one L40S each, NP=1 (repo convention). Idempotent: skips jobs whose
# .out already holds a final "!  total energy". Logs to sp02/run_sp02.log.
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
QE=/workspace/qe-7.5/bin/pw.x
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/landscape_442
LOG=$T/sp02/run_sp02.log

exec 9>$T/sp02/.sp02.lock; flock -n 9 || { echo "another sp02 driver already running"; exit 0; }

lane() {  # lane <gpu> <tag> [<tag>...]
  local gpu=$1; shift
  for tag in "$@"; do
    if grep -q '^!' $T/sp02/sp_l442_${tag}.out 2>/dev/null; then
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu SKIP $tag (already has final energy)" >> $LOG; continue
    fi
    echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu START $tag" >> $LOG
    cd $T && CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE \
      -nk 1 -nd 1 -nb 1 -nt 1 -inp sp02/sp_l442_${tag}.in > sp02/sp_l442_${tag}.out 2>&1
    if grep -q '^!' $T/sp02/sp_l442_${tag}.out; then
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu DONE $tag E=$(grep '^!' $T/sp02/sp_l442_${tag}.out | tail -1 | awk '{print $5}')" >> $LOG
    else
      echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu FAIL $tag (no final energy)" >> $LOG
    fi
  done
}

echo "[$(date '+%m-%d %H:%M:%S')] == sp02 batch start (5 jobs, 2 lanes) ==" >> $LOG
lane 0 vc_poly_lS3 st_poly_lS5 st_mol_lS6 &
P0=$!
lane 1 act_poly_lS2 act_mol_rd &
P1=$!
wait $P0 $P1
echo "[$(date '+%m-%d %H:%M:%S')] == sp02 batch END: $(ls $T/sp02/sp_l442_*.out 2>/dev/null | xargs grep -ls '^!' 2>/dev/null | wc -l)/5 with final energy ==" >> $LOG
