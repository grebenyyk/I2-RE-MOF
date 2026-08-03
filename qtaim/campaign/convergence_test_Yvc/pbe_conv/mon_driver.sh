#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_conv
for i in $(seq 1 80); do
  n=$(grep -c "DONE" $D/driver.log 2>/dev/null)
  [ "$n" -ge 9 ] && break
  sleep 30
done
echo "=== convergence driver: ALL 9 DONE ==="
cat $D/driver.log 2>/dev/null
echo "=== energies ==="
for f in pbe_kG pbe_k221 pbe_k663 pbe_poly_ec40 pbe_poly_ec60 pbe_empty_ec40 pbe_empty_ec60 pbe_gas_ec40 pbe_gas_ec60; do
  echo "$f: $(grep -aE '!.*total energy' $D/${f}.out 2>/dev/null | tail -1)"
done
