#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 90); do
  grep -q "JOB DONE\|EXIT=" $D/gas_hse_np1.out 2>/dev/null && break
  sleep 30
done
echo "=== final state ==="
grep -aE "convergence has been achieved|JOB DONE|EXIT=" $D/gas_hse_np1.out 2>/dev/null | tail -4
echo "=== ALL !total energy lines ==="
grep -anE "!.*total energy" $D/gas_hse_np1.out 2>/dev/null
echo "=== EXX check ==="
python3 -c "
import re
t=open('$D/gas_hse_np1.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  d=(H-P)*13.606*1000
  print(f'final E_HSE(I2 gas, -np1) = {H} Ry')
  print(f'diff vs PBE gas = {d:.1f} meV  (expect +50 to +200 meV if EXX present; ~0 = EXX still missing)')
else: print('no !total energy found')
"
