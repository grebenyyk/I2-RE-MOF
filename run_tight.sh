#!/usr/bin/env bash
# Sequential tight-BFGS continuation for the 5 strongest binders.
# One calc at a time, NP=2 (both GPUs). nohup me: nohup bash run_tight.sh > tight.log 2>&1 &
set -uo pipefail
cd /workspace/dimitry/calc/qe/i2_mof
for m in 13_Y 06_Tb 12_Tb 08_Gd 17_Cd; do
  echo "=== $(date '+%H:%M:%S') START $m ==="
  NP=2 OMP_NUM_THREADS=1 bash run.sh "$m/loaded_tight.in" "$m/loaded_tight.out"
  echo "=== $(date '+%H:%M:%S') DONE $m rc=$? ==="
done
echo "=== $(date '+%H:%M:%S') ALL TIGHT DONE ==="
