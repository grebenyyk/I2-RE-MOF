#!/usr/bin/env python3
"""Fix the activation-campaign killer: all 7 inputs shared prefix='activate_empty'
+ shared outdir='./tmp' (run.sh runs from BASE, no cd) -> concurrent jobs clobbered
each other's .save -> rc=127 / corrupted energies. Give each job a UNIQUE prefix and
a per-system absolute outdir so 2-at-a-time runs are fully isolated. Then clean the
polluted shared scratch and the failed outputs."""
import os, glob, shutil
BASE = "/workspace/dimitry/calc/qe/i2_mof"
SYSTEMS = ["13_Y", "06_Tb", "12_Tb", "05_Cd", "07_Cd", "10_Cd", "17_Cd"]
for s in SYSTEMS:
    f = os.path.join(BASE, s, "activate_empty", "relax.in")
    t = open(f).read()
    # unique prefix (was 'activate_empty' for ALL -> collision)
    t = t.replace("prefix = 'activate_empty'", f"prefix = '{s}_activate'")
    # per-system absolute outdir (was './tmp' = BASE/tmp shared)
    t = t.replace("outdir = './tmp'", f"outdir = '{BASE}/{s}/activate_empty/tmp'")
    open(f, "w").write(t)
    ok_pre = f"prefix = '{s}_activate'" in t
    ok_out = f"{s}/activate_empty/tmp'" in t
    print(f"{s}: unique_prefix={ok_pre} per_sys_outdir={ok_out}")
# remove polluted shared activate_empty.* scratch in BASE/tmp
removed = 0
for p in glob.glob(os.path.join(BASE, "tmp", "activate_empty*")):
    if os.path.isdir(p):
        shutil.rmtree(p)
    else:
        os.remove(p)
    removed += 1
print(f"removed {removed} polluted tmp/activate_empty* files")
# clean failed outputs so skip-logic + our eyes aren't confused
for s in SYSTEMS:
    o = os.path.join(BASE, s, "activate_empty", "relax.out")
    if os.path.exists(o):
        os.remove(o)
print("cleaned 7 failed relax.out")
