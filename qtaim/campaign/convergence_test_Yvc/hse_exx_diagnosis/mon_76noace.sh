#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 80); do
  grep -q "JOB DONE" $D/gas_hse_76_noace.out 2>/dev/null && break
  sleep 30
done
echo "=== QE 7.6 + ace=.false. + PAW gas HSE: DECISIVE RESULT ==="
echo "UNTESTED count: $(grep -ac 'compute_becpsi UNTESTED' $D/gas_hse_76_noace.out 2>/dev/null)"
echo "PAW_xx_nrg not-running count: $(grep -ac 'PAW_xx_nrg not running' $D/gas_hse_76_noace.out 2>/dev/null)"
echo "JOB DONE: $(grep -ac 'JOB DONE' $D/gas_hse_76_noace.out 2>/dev/null)"
echo "convergence: $(grep -ac 'convergence has been achieved' $D/gas_hse_76_noace.out 2>/dev/null)"
echo "=== energy progression (looking for -584 jump OR stable ~-757.8) ==="
grep -aE "total energy" $D/gas_hse_76_noace.out 2>/dev/null | grep -vE "^!" | tail -12 | awk '{print NR": "$0}'
echo "=== !total energy ==="
grep -aE "!.*total energy" $D/gas_hse_76_noace.out 2>/dev/null
echo "=== EXX check ==="
python3 -c "
import re
t=open('$D/gas_hse_76_noace.out',errors='ignore').read()
m=re.findall(r'!\s*total energy\s*=\s*(-?\d+\.\d+)',t)
if m:
  H=float(m[-1]); P=-757.79993434
  d=(H-P)*13.606*1000
  print(f'E_HSE(7.6 ace=.false. PAW I2 gas) = {H} Ry ; diff vs PBE = {d:.1f} meV')
  if abs(d) > 10: print('VERDICT: *** GENUINE EXX COMPUTED *** (becpsi works in 7.6)')
  else: print('VERDICT: still ~PBE — becpsi still not computing EXX')
else: print('no !total energy')
"
