#!/bin/bash
# Y_actc QTAIM redo driver: 3 SCFs -> pp/critic2/CT/ELF -> projwfc PDOS -> verify_bcp
# Runs in sequence; logs to run_qtaim.log. Correct geometry = final/Y-MOF_chain.POSCAR.
set -o pipefail
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_actc
QE=/workspace/qe-7.5
MPI=/workspace/openmpi; NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
PW="mpirun --allow-run-as-root -np 1 $QE/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1"
PP="$QE/bin/pp.x"; PROJ="$QE/bin/projwfc.x"
log(){ echo "[$(date +%H:%M:%S)] $*"; }

log "=== 1/5 SCFs (loaded GPU0, empty GPU1 parallel) ==="
CUDA_VISIBLE_DEVICES=0 $PW -inp scf_loaded.in > scf_loaded.out 2>&1 &
P0=$!
CUDA_VISIBLE_DEVICES=1 $PW -inp scf_empty.in > scf_empty.out 2>&1 &
P1=$!
wait $P0; L0=$?; wait $P1; L1=$?
log "loaded rc=$L0, empty rc=$L1"
log "=== i2only (GPU0) ==="
CUDA_VISIBLE_DEVICES=0 $PW -inp scf_i2only.in > scf_i2only.out 2>&1; L2=$?
log "i2only rc=$L2"
touch loaded.done empty.done i2only.done
# verify convergence
for c in loaded empty i2only; do
  E=$(grep -m1 '! *total energy' scf_${c}.out 2>/dev/null | awk '{print $NF}')
  log "  scf_${c}: E=${E:-FAIL}"
done

log "=== 2/5 analyze_campaign (pp cubes + deltarho + critic2 BCP/Bader CT + ELF) ==="
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign
python3 analyze_campaign.py Y_actc 2>&1 | tail -5

log "=== 3/5 projwfc PDOS ==="
cd Y_actc
$PROJ -inp proj.in > proj.out 2>&1; log "projwfc rc=$?"

log "=== 4/5 verify_bcp (BCP filter + corrected table) ==="
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign
python3 verify_bcp.py Y_actc 2>&1 | tail -3

log "=== 5/5 QTAIM DONE ==="
