#!/usr/bin/env python3
"""Consolidate per-system CT + BCP/EML tables into cross-system comparisons."""
import re, os
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=["Y_str","Tb_str","Gd_str","Y_act","Tb_act","Gd_act"]
METAL={"Y_str":"Y","Tb_str":"Tb","Gd_str":"Gd","Y_act":"Y","Tb_act":"Tb","Gd_act":"Gd"}
def read_ct(s):
    q={}; t=open(f"{BASE}/{s}/{s}_ct.txt").read().split("\n")[2:]
    for l in t:
        c=l.split()
        if len(c)>=4:
            try: q[c[0]]=float(c[2])
            except: pass
    return q
def read_bcp(s):
    bonds=[]; t=open(f"{BASE}/{s}/{s}_bcp_eml.txt").read().split("\n")[2:]
    for l in t:
        c=l.split()
        if len(c)>=8:
            try: bonds.append((c[0],float(c[1]),float(c[2]),float(c[3]),float(c[7])))
            except: pass
    return bonds
out=[]
out.append("="*72); out.append("CHARGE TRANSFER (Bader-basin δρ integral, e) — Σ over all atoms of species"); out.append("="*72)
out.append(f"{'system':<8}{'C(32)':>9}{'H(12)':>9}{'O(24)':>9}{'metal':>9}{'I(8)':>9}{'Σcheck':>9}")
for s in SYS:
    q=read_ct(s); m=METAL[s]
    out.append(f"{s:<8}{q.get('C',0):>+9.3f}{q.get('H',0):>+9.3f}{q.get('O',0):>+9.3f}{q.get(m,0):>+9.3f}{q.get('I',0):>+9.3f}{sum(q.values()):>+9.3f}")
out.append("")
out.append("="*72); out.append("BOND CPs + EML bond energies (Abramov; D_e=½|V|, kJ/mol) — I-related bonds"); out.append("="*72)
out.append(f"{'system':<8}{'nII':>5}{'ΣD_e II':>9}{'nIO':>5}{'ΣD_e IO':>9}{'nIM':>5}{'ΣD_e IM':>9}{'tot ΣD_e':>10}")
for s in SYS:
    b=read_bcp(s); m=METAL[s]; IM="".join(sorted(["I",m]))
    def grp(lab):
        bs=[x for x in b if x[0]==lab]; return len(bs),sum(x[4] for x in bs)
    nII,dII=grp("II"); nIO,dIO=grp("IO"); nIM,dIM=grp(IM)
    out.append(f"{s:<8}{nII:>5}{dII:>9.1f}{nIO:>5}{dIO:>9.1f}{nIM:>5}{dIM:>9.1f}{dII+dIO+dIM:>10.1f}")
out.append("")
out.append("="*72); out.append("I-I bond populations (per system): covalent (∇²ρ<0) vs closed-shell (∇²ρ>0)"); out.append("="*72)
out.append(f"{'system':<8}{'n cov':>6}{'ΣD_e cov':>10}{'n c-s':>6}{'ΣD_e c-s':>10}")
for s in SYS:
    b=read_bcp(s); ii=[x for x in b if x[0]=="II"]
    cov=[x for x in ii if x[3]<0]; cs=[x for x in ii if x[3]>=0]
    out.append(f"{s:<8}{len(cov):>6}{sum(x[4] for x in cov):>10.1f}{len(cs):>6}{sum(x[4] for x in cs):>10.1f}")
out.append("")
out.append("(ELF slices: <sys>_elf.png per system; per-bond detail in <sys>_bcp_eml.txt)")
open(f"{BASE}/consolidated_tables.txt","w").write("\n".join(out)+"\n")
print("\n".join(out))
