#!/usr/bin/env python3
"""Compute PBE-D3(BJ) dispersion for all landscape configs via the dftd3 driver (func=pbe).
Same configs as batch_hse06d3.py but with func=pbe. Saves to pbe_d3_energies.json."""
import re,os,glob,json,subprocess,numpy as np
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RUN=f"{D}/run_d3.py"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
env=os.environ.copy()
env["PATH"]=MPI+"/bin:"+NVHPC+"/compilers/bin:"+env.get("PATH","")
env["LD_LIBRARY_PATH"]=MPI+"/lib:"+NVHPC+"/compilers/lib:"+NVHPC+"/math_libs/lib:"+env.get("LD_LIBRARY_PATH","")
def classify(tag):
    m=re.match(r'rd_poly_Y(act|vc|st)$',tag)
    if m: return m.group(1)
    if re.match(r'hse_lpoly|hse_lmol',tag): return 'x'
    if re.match(r'hmol_Y(act|vc|st)_S1$',tag): return 'x'
    if re.match(r'rd_mol_Yact$',tag): return 'x'
    if re.match(r'hmol_Yvc$',tag): return 'x'
    return None
tags=sorted(set(os.path.basename(o)[:-4] for o in glob.glob(f"{D}/*.out") if classify(os.path.basename(o)[:-4])))
tags+=["rd_empty_Yst","rd_empty_Yact","rd_empty_Yvc","rd_gas"]
res={}
for tag in tags:
    inp=f"{D}/{tag}.in"
    if not os.path.exists(inp): res[tag]=None; continue
    r=subprocess.run(["python3",RUN,inp,"pbe","4"],capture_output=True,text=True,env=env,cwd=D)
    out=(r.stdout+r.stderr).strip().splitlines()
    out=out[-1] if out else ""
    try: res[tag]=float(out)
    except ValueError: res[tag]=f"FAIL:{out[-300:]}"
json.dump(res,open(f"{D}/pbe_d3_energies.json","w"),indent=1)
ok=sum(1 for v in res.values() if isinstance(v,float))
print(f"PBE-D3 computed: {ok}/{len(tags)}")
for t in ["rd_poly_Yst","rd_poly_Yact","rd_poly_Yvc","rd_empty_Yst","rd_empty_Yact","rd_empty_Yvc","rd_gas"]:
    print(f"  {t}: {res.get(t)}")
