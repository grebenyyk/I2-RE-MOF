#!/bin/bash
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
QE=/workspace/qe-7.5/bin/pw.x
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/st_vc
LOG=$T/run_rs_stvc.log
exec 9>$T/.rs_stvc.lock; flock -n 9 || exit 0
lane() { local gpu=$1; shift
  for tag in "$@"; do
    grep -q '^!' $T/${tag}.out 2>/dev/null && { echo "[$(date '+%H:%M:%S')] SKIP $tag" >> $LOG; continue; }
    echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu START $tag" >> $LOG
    cd $T && CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 $QE -nk 1 -nd 1 -nb 1 -nt 1 -inp ${tag}.in > ${tag}.out 2>&1
    grep -q '^!' $T/${tag}.out && echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu DONE $tag E=$(grep '^!' $T/${tag}.out | tail -1 | awk '{print $5}')" >> $LOG || echo "[$(date '+%m-%d %H:%M:%S')] GPU$gpu FAIL $tag" >> $LOG
  done; }
echo "[$(date '+%m-%d %H:%M:%S')] == rs_stvc batch start ==" >> $LOG
lane 0 rs_loaded_st_poly_lS5 rs_empty_st_poly_lS5 rs_i2only_st_poly_lS5 &
lane 1 rs_loaded_vc_poly_lS3 rs_empty_vc_poly_lS3 rs_i2only_vc_poly_lS3 &
wait
echo "[$(date '+%m-%d %H:%M:%S')] == rs_stvc batch END: $(ls $T/rs_*.out 2>/dev/null | xargs grep -ls '^!' 2>/dev/null | wc -l)/6 ==" >> $LOG
