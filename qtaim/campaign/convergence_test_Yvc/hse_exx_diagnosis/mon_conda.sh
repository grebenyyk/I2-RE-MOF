#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 50); do
  grep -q "JOB DONE\|EXIT=" $D/gas_hse_conda.out 2>/dev/null && break
  sleep 20
done
echo "=== final state ==="
grep -aE "convergence has been achieved|!.*total energy|Using ACE|ACE projected|JOB DONE|EXIT=|exxergy|EXX energy|highest occupied|lowest unoccupied" $D/gas_hse_conda.out 2>/dev/null | tail -12
echo "=== EXX check: E_HSE(conda CPU) vs E_PBE gas (-757.79993434) ==="
python3 -c "
import re
t=open('$D/gas_hse_conda.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  print(f'E_HSE(conda CPU I2 gas)={H}  diff vs PBE = {(H-P)*13.606*1000:.1f} meV')
else: print('no total energy found')
"
