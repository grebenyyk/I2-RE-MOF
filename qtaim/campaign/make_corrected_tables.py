#!/usr/bin/env python3
"""Consolidate the 6 per-system corrected BCP tables into one cross-system table
(text + docx) and an explicit changelog of the corrected rows. NEW files only."""
import os, re
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=[("Y_str","Y-MOF strained"),("Tb_str","Tb-MOF strained"),("Gd_str","Gd-MOF strained"),
     ("Y_act","Y-MOF activated"),("Tb_act","Tb-MOF activated"),("Gd_act","Gd-MOF activated")]

# ---------- old table (the one with the monotonicity violations) ----------
OLD={}
for blk in open(f"{BASE}/QTAIM_bond_table_all6.txt").read().split("Y-MOF")[1:]:
    pass
cur=None
for l in open(f"{BASE}/QTAIM_bond_table_all6.txt"):
    m=re.match(r"\s*(Y|Tb|Gd)-MOF\s+(strained|activated)",l)
    if m: cur=(m.group(1)+"_"+{"strained":"str","activated":"act"}[m.group(2)]); continue
    m=re.match(r"\s*(I\d+-I\d+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)",l)
    if m and cur: OLD.setdefault(cur,{})[m.group(1)]=(float(m.group(2)),float(m.group(3)),float(m.group(4)))

def parse_corr(sys):
    ii=[];io=[];sec=None
    for l in open(f"{BASE}/{sys}/{sys}_bcp_corrected.txt"):
        if "CORRECTED I-I" in l: sec="ii"; continue
        if "CORRECTED I-O" in l: sec="io"; continue
        c=l.split()
        if l.startswith("#") or not c or not c[0].startswith("I"): continue
        # pair d rho_b lap_b rho_A3 lap_A5 type [fld]
        if len(c)>=7:
            pair=c[0]; d=float(c[1]); rho=float(c[2]); lap=float(c[3]); rhoA=float(c[4]); lapA=float(c[5])
            typ=c[6]; fld=c[7] if len(c)>7 else ""
            (ii if sec=="ii" else io).append((pair,d,rho,lap,rhoA,lapA,typ,fld))
    return ii,io

# ---------- consolidated text ----------
L=["="*80,"CORRECTED QTAIM I-I bond table (all 6 systems; PAW/shell & near-nuclear artifacts removed)","="*80,
   "Re-located on clean valence field (plot_num=0) + rho17, pooled. Keep rho in [0.02,0.085] e/bohr3,",
   "lap<+0.3, not same-atom termini. Verified monotone (rho falls with d) + bond-path profile agrees.",
   "rho in e/A3, lap in e/A5. type: cov=covalent(lap<0), cs=closed-shell(lap>0).",
   "-"*80]
allcorr={}
for sys,name in SYS:
    ii,io=parse_corr(sys); allcorr[sys]=(ii,io)
    L.append(f"\n{name}")
    L.append(f"  {'bond':<8}{'d,A':>6}{'rho':>8}{'lap':>9}{'type':>5}  (vs old table)")
    for pair,d,rho,lap,rhoA,lapA,typ,fld in sorted(ii,key=lambda x:x[1]):
        old=OLD.get(sys,{}).get(pair)
        tag=""
        if old:
            od,orho,olap=old
            if abs(orho-rhoA)>0.05: tag=f"  OLD rho={orho:.3f} {'(SPIKE)' if orho>rhoA*1.6 else ''}"
        else:
            tag="  (new/not in old table)"
        t="cov" if lap<-0.02 else ("cs" if lap>0.02 else "bord")
        L.append(f"  {pair:<8}{d:>6.2f}{rhoA:>8.3f}{lapA:>9.3f}{t:>5}{tag}")
    # removed bonds (in old, not in new corrected)
    newpairs={p for p,*_ in ii}
    for pair,(od,orho,olap) in OLD.get(sys,{}).items():
        if pair not in newpairs and orho>0.4:
            L.append(f"  {pair:<8}{od:>6.2f}{orho:>8.3f}{olap:>9.3f}{'':>5}  REMOVED (was SPIKE, no real BCP)")
open(f"{BASE}/QTAIM_bond_table_corrected.txt","w").write("\n".join(L)+"\n")

# ---------- changelog of the 4(+more) corrected rows ----------
CH=["="*70,"CHANGELOG: rows corrected (old monotonicity-violating -> new)","="*70,
    f"{'system':<9}{'bond':<9}{'d,A':>6}{'old rho':>10}{'old lap':>10}  ->  {'new rho':>9}{'new lap':>9}  status"]
nfix=0
for sys,name in SYS:
    ii,io=allcorr[sys]; newmap={p:(d,rhoA,lapA) for p,d,rho,lap,rhoA,lapA,typ,fld in ii}
    newpairs=set(newmap)
    for pair,(od,orho,olap) in OLD.get(sys,{}).items():
        if pair in newmap:
            nd,nrho,nlap=newmap[pair]
            if abs(orho-nrho)>0.05:
                CH.append(f"{sys:<9}{pair:<9}{od:>6.2f}{orho:>10.3f}{olap:>10.3f}  ->  {nrho:>9.3f}{nlap:>9.3f}  corrected"); nfix+=1
        else:
            if orho>0.4:
                CH.append(f"{sys:<9}{pair:<9}{od:>6.2f}{orho:>10.3f}{olap:>10.3f}  ->  {'---':>9}{'---':>9}  removed (spike)"); nfix+=1
CH.append("-"*70); CH.append(f"{nfix} rows corrected across 6 systems.")
open(f"{BASE}/BCP_changelog.txt","w").write("\n".join(CH)+"\n")
print("\n".join(CH))
print("\nwrote QTAIM_bond_table_corrected.txt, BCP_changelog.txt")
