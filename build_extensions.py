#!/usr/bin/env python3
"""Build nstep=200 EXTENSION inputs for the 4 chain jobs that hit nstep=50.
Restart from their current (step-50) geometry extracted from loaded.out, so they
converge in a few more steps to the same tightness as the nstep=200 jobs."""
import re, os
BASE="/workspace/dimitry/calc/qe/i2_mof"
SYSTEMS=["06_Tb","17_Cd","12_Tb","05_Cd"]
for s in SYSTEMS:
    out=f"{BASE}/{s}/activated_i2/loaded.out"
    inp=f"{BASE}/{s}/activated_i2/loaded.in"
    txt=open(out,errors="ignore").read()
    # last ATOMIC_POSITIONS block: all atoms (framework frozen + I at step-50)
    blocks=txt.split("ATOMIC_POSITIONS")
    last=blocks[-1]
    pos=[]
    for l in last.splitlines():
        m=re.match(r"^\s*([A-Z][a-z]?)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",l)
        if m and m.group(1) in ("C","Gd","H","O","I","N","Cd","Tb","Y"):
            sym=m.group(1); ifp="1 1 1" if sym=="I" else "0 0 0"
            pos.append((sym,float(m.group(2)),float(m.group(3)),float(m.group(4)),ifp))
    if len(pos)<50: print(f"{s}: only {len(pos)} atoms parsed, skip"); continue
    # rebuild from loaded.in template, swapping positions + nstep + prefix/outdir
    it=open(inp).read()
    # bump nstep to 200 (already set, but ensure)
    it=re.sub(r"nstep\s*=\s*\d+","nstep = 200",it)
    it=it.replace(f"prefix = '{s}_act_i2_chain_uniform'",f"prefix = '{s}_act_i2_ext'").replace(f"prefix = '{s}_activated_i2_loaded'",f"prefix = '{s}_act_i2_ext'")
    it=re.sub(r"outdir = '.*?'",f"outdir = '{BASE}/{s}/activated_i2/tmp_ext'",it)
    # replace the ATOMIC_POSITIONS block positions
    body=it.split("ATOMIC_POSITIONS")[0]+"ATOMIC_POSITIONS crystal\n"
    for sym,x,y,z,ifp in pos:
        body+=f"{sym:<3s} {x: .12f} {y: .12f} {z: .12f} {ifp}\n"
    tail=it.split("K_POINTS")[1]
    body+="K_POINTS"+tail
    ext=f"{BASE}/{s}/activated_i2/loaded_ext.in"
    open(ext,"w").write(body)
    print(f"{s}: extension input -> {ext} ({len(pos)} atoms, nstep=200, restart from step-50)")
