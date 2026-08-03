#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_vs_hse_442
cd "$D"
for i in $(seq 1 40); do
  done3=1
  for t in pbe_poly_Yvc pbe_empty_Yvc pbe_gas; do
    grep -q "convergence has been achieved\|JOB DONE" $t.out 2>/dev/null || done3=0
  done
  [ "$done3" = 1 ] && { echo "ALL 3 DONE"; break; }
  sleep 20
done
echo "=== PBE-D3 energies ==="
for t in pbe_poly_Yvc pbe_empty_Yvc pbe_gas; do echo "$t: $(grep -aE '!.*total energy' $t.out 2>/dev/null | tail -1)"; done
echo "=== PBE-D3 ΔEads(4x4x2) vs HSE (-111.98) ==="
python3 -c "
import re
def E(t): 
    m=re.search(r'total energy\s*=\s*(-?\d+\.\d+)',open(t+'.out',errors='ignore').read()); return float(m.group(1)) if m else None
L=E('pbe_poly_Yvc');E0=E('pbe_empty_Yvc');g=E('pbe_gas')
print(f'  PBE loaded={L} empty={E0} gas={g}')
dads=(L-E0-4*g)/4*1312.75
print(f'  PBE-D3 ΔEads(4x4x2) = {dads:.2f} kJ/I2   (HSE-D3 = -111.98)   diff = {dads-(-111.98):+.2f}')
"
