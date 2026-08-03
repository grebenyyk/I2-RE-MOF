#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 60); do
  grep -q "JOB DONE\|DONE\b" $D/gas_hse_76.out 2>/dev/null && break
  sleep 30
done
echo "=== QE 7.6 PAW gas HSE: DECISIVE RESULT ==="
echo "UNTESTED warning count: $(grep -ac 'compute_becpsi UNTESTED' $D/gas_hse_76.out 2>/dev/null)"
echo "PAW_xx_nrg not-running count: $(grep -ac 'PAW_xx_nrg not running' $D/gas_hse_76.out 2>/dev/null)"
echo "JOB DONE: $(grep -ac 'JOB DONE' $D/gas_hse_76.out 2>/dev/null)"
echo "convergence: $(grep -ac 'convergence has been achieved' $D/gas_hse_76.out 2>/dev/null)"
echo "=== !total energy lines ==="
grep -anE "!.*total energy" $D/gas_hse_76.out 2>/dev/null
echo "=== energy progression in EXX cycle (if any) ==="
grep -aE "total energy" $D/gas_hse_76.out 2>/dev/null | grep -vE "^!" | tail -10 | awk '{print NR": "$0}'
echo "=== EXX check ==="
python3 -c "
import re
t=open('$D/gas_hse_76.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  print(f'E_HSE(7.6 PAW I2 gas) = {H} Ry')
  print(f'diff vs PBE = {(H-P)*13.606*1000:.1f} meV')
  print('VERDICT: ' + ('GENUINE EXX (differs from PBE)' if abs((H-P)*13.606*1000) > 10 else 'still ~PBE (EXX missing/tiny)'))
else: print('no !total energy found')
"
