#!/bin/bash
# One-command resume for the landscape_442 campaign. Safe to re-run (flock-guarded);
# the driver scans .out files and skips converged configs, reseeds partials.
# Run this on any node (after the current one dies) to continue the campaign.
set -e
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/landscape_442
# driver writes run.log internally (LOG=fopen append); send stdout to /dev/null to avoid dup
nohup python3 run_l442.py >/dev/null 2>&1 &
disown
echo "landscape_442 driver launched (PID $!). Tail: tail -f run.log"
