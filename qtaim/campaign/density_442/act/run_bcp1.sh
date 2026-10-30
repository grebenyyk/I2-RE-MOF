#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
export OMP_NUM_THREADS=1
ulimit -s unlimited 2>/dev/null || true
C2=/workspace/dimitry/calc/qe/i2_mof/tools/critic2_run.sh
for j in bcp_act_poly_lS2 bcp_act_mol_rd; do
  echo "RUN $j (OMP=1) $(date +%H:%M:%S)"
  $C2 $j.cri > $j.out 2>&1 && echo "OK $j" || echo "FAIL $j"
done
echo "BCP DONE"
