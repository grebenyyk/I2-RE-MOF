#!/bin/bash
# Supervisor: keep run_442.py alive until all 48 done. Relaunches if it dies.
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/raman_442
while true; do
  ndone=$(python3 -c "
import re
n=0
for t in open('tags.txt').read().split('\n'):
    t=t.split()
    if len(t)<4: continue
    tag=t[3]; o=f'disp_{tag}.out'
    try:
        s=open(o).read()
        if 'convergence has been achieved' in s and 'Forces acting on atoms' in s: n+=1
    except: pass
print(n)")
  echo "[$(date +%H:%M:%S)] supervisor: $ndone/48 done" >> supervisor.log
  [ "$ndone" -ge 48 ] && { echo "ALL 48 DONE" >> supervisor.log; break; }
  if ! pgrep -f "python3 run_442.py" >/dev/null 2>&1; then
    echo "[$(date +%H:%M:%S)] supervisor: run_442.py dead, relaunching" >> supervisor.log
    python3 run_442.py >> run.live 2>&1 &
  fi
  sleep 90
done
