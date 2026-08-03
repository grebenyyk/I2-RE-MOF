#!/usr/bin/env bash
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
for i in $(seq 1 40); do
  grep -q "JOB DONE\|EXIT=" $D/gas_hse_cpu.out 2>/dev/null && break
  sleep 20
done
echo "=== final state ==="
grep -aE "convergence has been achieved|!.*total energy|Using ACE|ACE projected|JOB DONE|EXIT=|exxergy|EXX" $D/gas_hse_cpu.out 2>/dev/null | tail -10
