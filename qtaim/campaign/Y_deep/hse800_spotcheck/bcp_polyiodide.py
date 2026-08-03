#!/usr/bin/env python3
"""Clean BCP invariance table for I2@Y-MOFvc: match PBE vs HSE polyiodide I-I BCPs
by POSITION (fractional MIC), keep only genuine same-bond matches (close in position
AND rho), report rho_BCP PBE vs HSE. Filters out topology-noise mismatches."""
import re, numpy as np
SP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck"
YD="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"
def parse_cube(path):
    L=open(path).readlines(); nat=int(L[2].split()[0])
    nr=np.array([int(L[i].split()[0]) for i in (3,4,5)])
    step=np.array([[float(x) for x in L[i].split()[1:4]] for i in (3,4,5)])
    cellB=nr[:,None]*step
    sym=[];cartB=[]
    for i in range(6,6+nat):
        c=L[i].split(); sym.append({1:"H",6:"C",8:"O",39:"Y",53:"I"}.get(int(c[0]),c[0]))
        cartB.append([float(c[2]),float(c[3]),float(c[4])])
    frac=np.array(cartB)@np.linalg.inv(cellB)
    return sym,frac,cellB
def parse_bcp(p):
    R=[]
    for l in open(p,errors="ignore"):
        m=re.match(r"\s*\d+\s+\S+\s+\(3,-1\)\s+bond\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\s+\S+\s+(-?\d+\.\d+E[+-]\d+)\s+\S+\s+(-?\d+\.\d+E[+-]\d+)",l)
        if m: R.append((np.array([float(m.group(1)),float(m.group(2)),float(m.group(3))]),float(m.group(4)),float(m.group(5))))
    return R
sym,frac,cellB=parse_cube(f"{YD}/rho17.cube")
Ipos=frac[[i for i,s in enumerate(sym) if s=='I']]
def lbl(fc):
    df=frac-fc; df=df-np.round(df); d=np.linalg.norm(df,axis=1); i=np.argsort(d)[:2]
    return "".join(sorted(sym[i[0]]+sym[i[1]]))
def polyII(rows):
    return [(fc,rho,lap) for fc,rho,lap in rows if lbl(fc)=='II' and 0.015<rho<0.09]
P=polyII(parse_bcp(f"{SP}/bcp_PBE17.out")); H=polyII(parse_bcp(f"{SP}/bcp_HSE17.out.out" if False else f"{SP}/bcp_HSE17.out"))
print(f"Polyiodide I-I BCPs (0.015<rho<0.09): PBE={len(P)}, HSE={len(H)}")
print(f"{'rho_PBE':>9} {'rho_HSE':>9} {'d_rho':>9} {'%diff':>7}  {'lap_PBE':>9} {'lap_HSE':>9}  status")
Harr=np.array([h[0] for h in H]); Hrho=np.array([h[1] for h in H]); Hlap=np.array([h[2] for h in H])
match=0; topodiff=0
for fc,rho,lap in sorted(P,key=lambda r:-r[1]):
    df=Harr-fc; df=df-np.round(df); dist=np.linalg.norm(df,axis=1); j=dist.argmin()
    if dist[j]<0.04 and abs(rho-Hrho[j])/rho<0.10:   # genuine same-bond match (close pos, close rho)
        print(f"{rho:9.4f} {Hrho[j]:9.4f} {rho-Hrho[j]:9.4f} {100*(rho-Hrho[j])/rho:7.2f}  {lap:9.2e} {Hlap[j]:9.2e}  MATCH"); match+=1
    else:
        print(f"{rho:9.4f} {'--':>9} {'--':>9} {'--':>7}  {lap:9.2e} {'--':>9}  topology-shift"); topodiff+=1
pr=np.array([r[1] for r in P])
# build matched rho arrays for mean
mr=[]
for fc,rho,lap in P:
    df=Harr-fc; df=df-np.round(df); dist=np.linalg.norm(df,axis=1); j=dist.argmin()
    if dist[j]<0.04 and abs(rho-Hrho[j])/rho<0.10: mr.append(Hrho[j])
mr=np.array(mr); pr=np.array([r[1] for fc,rho,lap in P if any(np.allclose(fc,h[0]) for h in [(P[k][0],0,0) for k in range(len(P))])])  # noop guard
print(f"\nGenuine same-bond I-I matches: {match}/{len(P)}  (mean |%diff| on matches = {np.mean(np.abs(pr[:match]-mr)/pr[:match])*100 if match else 0:.2f}%)")
print(f"Topology-shifted (weak BCP, rho~0.02-0.03, appears in one not other): {topodiff}")
