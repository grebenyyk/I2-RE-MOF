#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
C2=/workspace/dimitry/calc/qe/i2_mof/tools/critic2
for j in "bcp_act_poly_lS2" "ct_act_poly_lS2" "bcp_act_mol_rd" "ct_act_mol_rd" "ct_2x_act_poly_lS2"; do
  [ -s "$j.out" ] && { echo "SKIP $j"; continue; }
  echo "RUN $j $(date +%H:%M:%S)"
  $C2 $j.cri > $j.out 2>&1 && echo "OK $j" || echo "FAIL $j"
done
echo "CRITIC2 ALL DONE"
