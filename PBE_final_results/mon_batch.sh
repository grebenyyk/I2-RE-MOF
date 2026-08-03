#!/usr/bin/env bash
PBE=/workspace/dimitry/calc/qe/i2_mof/PBE_final_results
while true; do
  n=$(grep -c "DONE" $PBE/run_all_pbe.log 2>/dev/null)
  echo "[$(date +%H:%M)] $n/69 PBE SCFs done" 
  grep -q "ALL.*COMPLETE" $PBE/run_all_pbe.log 2>/dev/null && break
  sleep 120
done
echo "=== ALL 69 PBE+D3 SCFs COMPLETE ==="
echo "converged: $(grep -l 'convergence has been achieved' $PBE/*/*.out 2>/dev/null | wc -l)/69"
echo "sample output (rd_poly_Yvc):"
head -8 $PBE/adsorption_9systems/rd_poly_Yvc.out 2>/dev/null
