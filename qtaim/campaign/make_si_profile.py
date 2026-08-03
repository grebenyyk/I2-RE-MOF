#!/usr/bin/env python3
"""Build an SI table: critic2 rho_BCP vs independent bond-path rho_profile,
for every accepted I-I and I-O bond across the 6 systems. Parses the existing
{sys}_bcp_verify.txt ([B]+[C] for I-I, [D] for I-O). NEW files."""
import re
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=[("Y_str","Y-MOF strained"),("Tb_str","Tb-MOF strained"),("Gd_str","Gd-MOF strained"),
     ("Y_act","Y-MOF activated"),("Tb_act","Tb-MOF activated"),("Gd_act","Gd-MOF activated")]

def parse(sys):
    """return ii=[(pair,d,rho_bcp,rho_prof)], io=[same]. rho in e/bohr3."""
    t=open(f"{BASE}/{sys}/{sys}_bcp_verify.txt").read()
    ii_b={}; ii_p={}; io=[]
    sec=None
    for l in t.splitlines():
        s=l.strip()
        if s.startswith("[B]"): sec="B"
        elif s.startswith("[C]"): sec="C"
        elif s.startswith("[D]"): sec="D"
        elif s.startswith("[") or s.startswith("="): sec=None
        if not s or s.startswith("pair") or s.startswith("[") or s.startswith("="): continue
        c=s.split()
        if sec=="B" and len(c)>=4 and c[0].startswith("I"):
            ii_b[c[0]]=(float(c[1]),float(c[2]))   # pair->(d,rho_bcp)
        elif sec=="C" and len(c)>=3 and c[0].startswith("I"):
            ii_p[c[0]]=float(c[2])                  # pair->rho_prof
        elif sec=="D" and len(c)>=4 and c[0].startswith("I"):
            io.append((c[0],float(c[1]),float(c[2]),float(c[3])))  # pair,d,rho_bcp,rho_prof
    ii=[(p,d,b,ii_p[p]) for p,(d,b) in ii_b.items() if p in ii_p]
    return ii,io

L=["SI: verification of QTAIM bond critical points (BCPs) on the PAW valence density","",
 "For each accepted I-I and I-O BCP (located with critic2 on plot_num=17), the density at the",
 "bond critical point (rho_BCP) is compared to an independent estimate (rho_prof): the minimum",
 "of rho(r) sampled along the straight internuclear line, t in [0.3,0.7], linear interpolation,",
 "on the PAW-augmentation-free valence field (plot_num=0). rho_prof uses none of critic2's",
 "CP-localization machinery. Agreement confirms the reported rho_BCP/grad^2 rho are intrinsic",
 "to the bonds and not artefacts of the PAW reconstruction. All rho in e/bohr^3.","",
 "I-I bonds:","",
 f"{'system':<20}{'bond':<9}{'d(A)':>6}{'rho_BCP':>10}{'rho_prof':>10}{'%diff':>8}"]
import statistics
diffs=[]
for sys,name in SYS:
    ii,io=parse(sys)
    for pair,d,b,p in sorted(ii,key=lambda x:x[1]):
        pd=100*abs(b-p)/b; diffs.append(pd)
        L.append(f"{name:<20}{pair:<9}{d:>6.2f}{b:>10.4f}{p:>10.4f}{pd:>7.1f}%")
L+=[f"\nI-O bonds:","",""]
L.append(f"{'system':<20}{'bond':<9}{'d(A)':>6}{'rho_BCP':>10}{'rho_prof':>10}{'%diff':>8}")
for sys,name in SYS:
    ii,io=parse(sys)
    for pair,d,b,p in sorted(io,key=lambda x:x[1]):
        pd=100*abs(b-p)/max(b,1e-9); diffs.append(pd)
        L.append(f"{name:<20}{pair:<9}{d:>6.2f}{b:>10.4f}{p:>10.4f}{pd:>7.1f}%")
L+=["",f"Mean |%diff| across all bonds: {statistics.mean(diffs):.2f}%, "
    f"max {max(diffs):.1f}% (n={len(diffs)})."]
open(f"{BASE}/SI_bcp_profile_check.txt","w").write("\n".join(L)+"\n")
print("\n".join(L[-3:]))
print(f"wrote SI_bcp_profile_check.txt ({len(diffs)} bonds)")
