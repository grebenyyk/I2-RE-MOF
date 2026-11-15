#!/bin/bash
set -uo pipefail
MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/adduct"
export PATH="$MPI/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:${LD_LIBRARY_PATH:-}"
cd "$B"
for P in pp_rho17_adduct pp_val_adduct pp_elf_adduct pp_val_i2only_adduct; do
  echo "[$(date '+%H:%M:%S')] pp.x $P"
  mpirun --allow-run-as-root -np 2 "$QE/bin/pp.x" -inp $P.in > $P.out 2>&1
  echo "  rc=$? cube=$(ls -la --block-size=M ${P#pp_}*.cube 2>/dev/null | awk '{print $5}')"
done
echo PP DONE
