#!/usr/bin/env bash
# Poll HSE loaded+empty until both charge-density.dat are written (at SCF convergence,
# before any PAW-EXX post-step hang), then copy saves to captured/. Cap ~30 min.
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck
cd "$D"
captured=0
for i in $(seq 1 60); do   # 60 x 30s = 30 min
  L=$([ -f tmp_loaded/hse800_vc.save/charge-density.dat ] && echo y || echo n)
  E=$([ -f tmp_empty/hse800_vc_empty.save/charge-density.dat ] && echo y || echo n)
  LC=$(grep -ac "convergence has been achieved" loaded.out 2>/dev/null)
  EC=$(grep -ac "convergence has been achieved" empty.out 2>/dev/null)
  LI=$(grep -ac "estimated scf accuracy" loaded.out 2>/dev/null)
  EI=$(grep -ac "estimated scf accuracy" empty.out 2>/dev/null)
  LACC=$(grep -a "estimated scf accuracy" loaded.out 2>/dev/null | tail -1 | awk '{print $NF}')
  echo "[$(date +%H:%M:%S)] loaded: conv=$LC dens=$L iter=$LI acc=$LACC | empty: conv=$EC dens=$E iter=$EI"
  if [ "$L" = y ] && [ "$E" = y ]; then captured=1; echo "BOTH DENSITIES CAPTURED"; break; fi
  sleep 30
done
mkdir -p captured
[ "$L" = y ] && cp -r tmp_loaded/hse800_vc.save captured/hse800_vc.save 2>/dev/null
[ "$E" = y ] && cp -r tmp_empty/hse800_vc_empty.save captured/hse800_vc_empty.save 2>/dev/null
echo "=== convergence + energy ==="
grep -aiE "convergence has been achieved|!.*total energy" loaded.out 2>/dev/null | tail -2
grep -aiE "convergence has been achieved|!.*total energy" empty.out 2>/dev/null | tail -2
echo "=== captured densities ==="
ls -la captured/*/charge-density.dat 2>/dev/null
echo "=== captured save contents (loaded) ==="
ls captured/hse800_vc.save/ 2>/dev/null | head
[ "$captured" = 1 ] && echo "MONITOR_DONE_OK" || echo "MONITOR_DONE_TIMEOUT"
