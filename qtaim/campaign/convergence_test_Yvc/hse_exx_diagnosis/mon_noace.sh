#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 80); do
  grep -q "JOB DONE\|EXIT=" $D/gas_hse_noace.out 2>/dev/null && break
  sleep 30
done
echo "=== final state ==="
grep -aE "convergence has been achieved|JOB DONE|EXIT=" $D/gas_hse_noace.out 2>/dev/null | tail -3
echo "=== ALL !total energy lines ==="
grep -anE "!.*total energy" $D/gas_hse_noace.out 2>/dev/null
echo "=== genuine HSE gas vs PBE ==="
python3 -c "
import re
t=open('$D/gas_hse_noace.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  print(f'genuine E_HSE(I2 gas, non-ACE) = {H} Ry')
  print(f'diff vs PBE = {(H-P)*13.606*1000:.1f} meV  (EXX contribution)')
else: print('no !total energy yet')
"
