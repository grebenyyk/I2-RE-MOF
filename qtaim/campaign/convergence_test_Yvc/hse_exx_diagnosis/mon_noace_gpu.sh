#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 40); do
  grep -q "JOB DONE\|EXIT=" $D/gas_hse_noace_gpu.out 2>/dev/null && break
  sleep 20
done
echo "=== non-ACE GPU gas HSE: final state ==="
grep -aE "convergence has been achieved|JOB DONE|EXIT=|Using ACE" $D/gas_hse_noace_gpu.out 2>/dev/null | tail -4
echo "=== !total energy lines ==="
grep -anE "!.*total energy" $D/gas_hse_noace_gpu.out 2>/dev/null
echo "=== total cpu time ==="
grep -aE "PWSCF.*CPU" $D/gas_hse_noace_gpu.out 2>/dev/null | tail -1
echo "=== EXX check (non-ACE GPU vs PBE -757.79993434) ==="
python3 -c "
import re
t=open('$D/gas_hse_noace_gpu.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  print(f'E_HSE(non-ACE GPU I2 gas) = {H} Ry ; diff vs PBE = {(H-P)*13.606*1000:.1f} meV')
else: print('no !total energy')
"
