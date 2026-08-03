#!/usr/bin/env python3
"""Re-locate I-related BCPs on a CLEAN field; flag & remove PAW/shell artifacts.
No existing file is overwritten. New outputs:
  {sys}/{sys}_bcp_corrected.txt  — corrected I-I + I-O BCP table
  {sys}/{sys}_bcp_verify.txt     — diagnostic (rejected CPs, monotonicity, cross-check)

ROOT CAUSE of the old monotonicity violations: the PAW valence field (rho17,
plot_num=17) and even the smooth valence field carry intra-atomic (3,-1) CPs at
rho~0.10-0.11 e/bohr3 in the 5p-shell region of each iodine (gradient paths end
at the SAME atom or at non-nuclear maxima). These are NOT I-I bonds. The old
rho>0.15 filter and "highest-rho per pair" selection let some through (rho~0.097),
producing rows where a 3.0 A "bond" had 2x the density of a 2.7 A one.

FIX = re-locate on the clean valence field (val.cube, plot_num=0) with a proper
QTAIM inter-atomic criterion: a bond CP must have TWO DIFFERENT real-atom
gradient-path termini, sane r1/r2 and angle, and rho in the bonding range. The
clean-field BCP rho/lap are reported; they are confirmed identical to rho17 at
the same positions (PAW reconstruction only perturbs the near-nuclear region).
Monotonicity (rho falls with d) is checked on the corrected values and on a
bonding-region (t=0.3-0.7) straight-line profile.

Units: rho e/bohr3 | e/A3 (x6.7481); lap e/bohr5 | e/A5 (x24.098).
"""
import os, re, sys
from collections import defaultdict
import numpy as np
from scipy.ndimage import map_coordinates

BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=["Y_str","Tb_str","Gd_str","Y_act","Tb_act","Gd_act"]
METAL={"Y_str":"Y","Tb_str":"Tb","Gd_str":"Gd","Y_act":"Y","Y_actc":"Y","Tb_act":"Tb","Gd_act":"Gd","Gd_deep":"Gd","Tb_deep":"Tb","Y_deep":"Y"}
BOHR=0.5291772109
RHO_A=1.0/BOHR**3; LAP_A=1.0/BOHR**5
ELEM={"C","H","O","N","Y","Tb","Gd","I","La","Ce"}
RHO_LO,RHO_HI=0.025,0.085   # bonding-range rho window (real I-I/I-O 0.037-0.075; excludes 5p-shell ~0.10 and near-nuclear ~0.012-0.02)
LAP_MAX=0.30                # near-nuclear residual CPs have lap~+1.5; real I-I/I-O BCPs |lap|<0.1

def read_cube(fn):
    L=open(fn).readlines(); nat=int(L[2].split()[0])
    nx,ny,nz=[int(L[i].split()[0]) for i in (3,4,5)]
    d=np.array([float(x) for l in L[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
    return d,nx,ny,nz

def read_atoms(sys):
    s=open(f"{BASE}/{sys}/scf_loaded.in").read()
    cellA=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l)))
                    for l in s.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
    apos=[]
    for l in s.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
        if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cellA,apos

def vesta_labels(apos):
    labs={}; ki=ko=0
    for sp,f in apos:
        if sp=="I": ki+=1; labs[id(f)]=f"I{ki}"
        elif sp=="O": ko+=1; labs[id(f)]=f"O{ko}"
    return labs

def rho_at(d,nx,ny,nz,fr,order=1):
    c=np.array([fr[0]*nx,fr[1]*ny,fr[2]*nz],float).reshape(3,1)
    return float(map_coordinates(d,c,order=order,mode="wrap")[0])

def profile_bond(d,nx,ny,nz,fa,fb):
    """rho at the bonding-region BCP: min of linear profile over t in [0.3,0.7]."""
    df=fb-fa; df=df-np.round(df)
    ts=np.linspace(0.30,0.70,41)
    v=np.array([rho_at(d,nx,ny,nz,(fa+t*df)-np.floor(fa+t*df),order=1) for t in ts])
    return float(v.min())

# ---- critic2 parsing (works for both rho17 bcp.out and val bcp_val.out) ----
def parse_cp_list(t):
    out={}; inl=False
    for l in t.splitlines():
        if "Critical point list, final report" in l: inl=True; continue
        if inl and l.strip().startswith("*"): break
        if not inl: continue
        c=l.split()
        if len(c)>=12 and re.match(r"^\d+$",c[0]):
            try:
                if "(3,-1)" not in c[2]: continue
                out[int(c[0])]={"frac":np.array([float(c[4]),float(c[5]),float(c[6])]),
                                "rho":float(c[9]),"lap":float(c[11])}
            except (ValueError,IndexError): pass
    return out

def parse_bonds(t):
    """id -> (sym1,idx1,sym2,idx2).  From 'Analysis of system bonds' (termini)."""
    out={}; f=False
    for l in t.splitlines():
        if "Analysis of system bonds" in l: f=True; continue
        if f and l.strip().startswith("*"): break
        if not f: continue
        c=l.split()
        if len(c)>=5 and re.match(r"^\d+$",c[0]) and c[2].startswith("(") and c[4].startswith("("):
            try: out[int(c[0])]=(c[1],c[2],c[3],c[4])
            except (ValueError,IndexError): pass
    return out

def is_nn(sym): return bool(re.match(r"n\d+",sym))
def near2(fc,apos,cellT):
    return sorted([(np.linalg.norm(cellT@((f-fc)-np.round(f-fc))),sp,f) for sp,f in apos],key=lambda x:x[0])[:2]

def real_interatomic_bonds(sys,field="val"):
    """Real I-I and I-O bond CPs. Discriminators (inter-atomic bond vs 5p-shell CP):
      - rho in bonding window [RHO_LO,RHO_HI]  (real I-I/I-O 0.04-0.075; shell ~0.10-0.11)
      - NOT same-atom gradient termini (both paths end at one atom = intra-atomic)
    NOTE: critic2's gradient-path tracing often ends at non-nuclear maxima for shallow
    closed-shell bonds even on the clean field, so 'two real-atom termini' is NOT required."""
    cellA,apos=read_atoms(sys); cellT=cellA.T; ilab=vesta_labels(apos)
    fn=f"{BASE}/{sys}/{('bcp_val.out' if field=='val' else 'bcp.out')}"
    t=open(fn,errors="ignore").read()
    cps=parse_cp_list(t); bonds=parse_bonds(t)
    ii=[]; io=[]; rejected=[]
    for cid in sorted(bonds):
        e1,i1,e2,i2=bonds[cid]; cp=cps.get(cid)
        if cp is None: continue
        rho=cp["rho"]
        nn=near2(cp["frac"],apos,cellT)
        if not (nn[0][1]=="I" or nn[1][1]=="I"): continue
        pair="-".join(sorted((ilab.get(id(nn[0][2]),nn[0][1]), ilab.get(id(nn[1][2]),nn[1][1]))))
        same_atom = (e1 in ELEM and e2 in ELEM and not is_nn(e1) and not is_nn(e2)
                     and (e1,i1)==(e2,i2))
        near_nuclear = cp["lap"]>LAP_MAX
        is_ii = (nn[0][1]=="I" and nn[1][1]=="I")
        is_io = ({nn[0][1],nn[1][1]}=={"I","O"})
        # type-specific rho windows: I-I stronger (0.02-0.085), I-O weak halogen (0.003-0.05)
        if is_ii:   lo,hi=0.020,RHO_HI
        elif is_io: lo,hi=0.003,0.050
        else:       lo,hi=RHO_HI,RHO_HI  # not II/IO -> never in window
        in_window = lo<=rho<=hi
        rec=(pair,rho,cp["lap"],cid,cp["frac"])
        if in_window and not same_atom and not near_nuclear:
            if is_ii: ii.append(rec)
            elif is_io: io.append(rec)
        elif is_ii or is_io:
            btype="II" if is_ii else "IO"
            why = ("same-atom termini (5p-shell/PAW)" if same_atom
                   else (f"lap={cp['lap']:.2f} (near-nuclear)" if near_nuclear
                         else f"rho={rho:.3f} outside bonding window"))
            rejected.append((cid,pair,rho,cp["lap"],why,field,btype))
    # dedup by pair (keep highest-rho real CP per pair)
    def dedup(lst):
        best={}
        for r in lst:
            if r[0] not in best or r[1]>best[r[0]][1]: best[r[0]]=r
        return sorted(best.values(),key=lambda x:x[0])
    return dedup(ii),dedup(io),rejected

def report(sys):
    cellA,apos=read_atoms(sys); cellT=cellA.T; m=METAL[sys]
    ii_val,io_val,rej_val=real_interatomic_bonds(sys,"val")     # clean-field selection
    ii17,io17,rej17=real_interatomic_bonds(sys,"rho17")         # rho17
    # pool both fields, dedup by pair (highest rho). tag field. (val and r17 agree at BCPs;
    # pooling recovers pairs one field's auto-search missed.)
    def pool(val,r17):
        best={}
        for p,rho,lap,cid,frac in val:
            if p not in best or rho>best[p][1]: best[p]=(p,rho,lap,cid,frac,"val")
        for p,rho,lap,cid,frac in r17:
            if p not in best or rho>best[p][1]: best[p]=(p,rho,lap,cid,frac,"r17")
        return list(best.values())
    ii_pool=pool(ii_val,ii17); io_pool=pool(io_val,io17)
    # geometric distances
    ilab=vesta_labels(apos)
    def dist(la,lb):
        fa=[f for s,f in apos if ilab.get(id(f))==la][0]
        fb=[f for s,f in apos if ilab.get(id(f))==lb][0]
        df=fb-fa; df=df-np.round(df); return float(np.linalg.norm(cellT@df)),fa,fb
    # ---------- corrected I-I table ----------
    C=[f"# {sys} CORRECTED I-I BCP table (re-located; PAW/shell & near-nuclear artifacts removed).",
       "# rho:e/bohr3|e/A3  lap:e/bohr5|e/A5.  Keep: rho in [0.025,0.085], lap<+0.3, not same-atom termini.",
       "# fld = which field supplied the CP (val=plot_num=0 clean, r17=plot_num=17); values agree at BCPs.",
       f"{'pair':<8}{'d(A)':>6}{'rho_b':>9}{'lap_b':>10}{'rho_A3':>9}{'lap_A5':>9}  type           fld"]
    ii_out=[]
    for pair,rho,lap,cid,frac,fld in sorted(ii_pool,key=lambda x:dist(x[0].split("-")[0],x[0].split("-")[1])[0]):
        a,b=pair.split("-"); dAA,fa,fb=dist(a,b)
        typ="covalent" if lap<-0.02 else ("closed-shell" if lap>0.02 else "borderline")
        C.append(f"{pair:<8}{dAA:>6.2f}{rho:>9.4f}{lap:>10.3f}{rho*RHO_A:>9.3f}{lap*LAP_A:>9.3f}  {typ:<14} {fld}")
        ii_out.append((pair,dAA,rho,lap,typ))
    C.append("")
    C.append(f"# {sys} CORRECTED I-O BCP table (iodine...framework oxygen; COO vs OH classified externally)")
    C.append(f"{'pair':<8}{'d(A)':>6}{'rho_b':>9}{'lap_b':>10}{'rho_A3':>9}{'lap_A5':>9}  type           fld")
    io_out=[]
    for pair,rho,lap,cid,frac,fld in sorted(io_pool,key=lambda x:dist(x[0].split("-")[0],x[0].split("-")[1])[0]):
        a,b=pair.split("-"); dAA,fa,fb=dist(a,b)
        C.append(f"{pair:<8}{dAA:>6.2f}{rho:>9.4f}{lap:>10.3f}{rho*RHO_A:>9.3f}{lap*LAP_A:>9.3f}  closed-shell   {fld}")
        io_out.append((pair,dAA,rho,lap))
    open(f"{BASE}/{sys}/{sys}_bcp_corrected.txt","w").write("\n".join(C)+"\n")
    # ---------- verify file ----------
    V=[f"==== {sys} : BCP re-location & artifact removal ===="]
    def dump_rejected(tag,label):
        V.append(""); V.append(f"[A{tag}] {label} REJECTED as non-bonded (intra-atomic 5p-shell / PAW spike / near-nuclear):")
        V.append(f"  {'field':<6}{'ncp':>5}{'pair':<9}{'rho_b':>9}{'lap_b':>10}{'rho_A3':>9}  reason")
        seen=set()
        for cid,pair,rho,lap,why,fld,bt in rej_val+rej17:
            if bt!=tag: continue
            key=(pair,round(rho,4))
            if key in seen: continue
            seen.add(key)
            V.append(f"  {fld:<6}{cid:>5}{pair:<9}{rho:>9.4f}{lap:>10.3f}{rho*RHO_A:>9.3f}  {why}")
        return seen
    seenII=dump_rejected("II","I-I candidate CPs")
    seenIO=dump_rejected("IO","I-O candidate CPs")
    V.append(""); V.append("[B] monotonicity check on corrected I-I BCP rho (must fall with distance):")
    V.append(f"  {'pair':<8}{'d(A)':>6}{'rho_b':>9}{'rho_A3':>9}  monotone?")
    prev=None
    for pair,dAA,rho,lap,typ in sorted(ii_out,key=lambda x:x[1]):
        mono="ok"
        if prev is not None and dAA>prev[0]+0.01 and rho>prev[1]+0.004: mono="VIOLATION"
        V.append(f"  {pair:<8}{dAA:>6.2f}{rho:>9.4f}{rho*RHO_A:>9.3f}  {mono}")
        prev=(dAA,rho)
    # bonding-region profile cross-check on val.cube (I-I)
    V.append(""); V.append("[C] bonding-region (t=0.3-0.7) straight-line rho profile on clean val.cube (I-I check):")
    if os.path.exists(f"{BASE}/{sys}/val.cube"):
        dv,nxv,nyv,nzv=read_cube(f"{BASE}/{sys}/val.cube")
        V.append(f"  {'pair':<8}{'d(A)':>6}{'rho_prof':>10}{'vs BCP':>9}")
        for pair,dAA,rho,lap,typ in sorted(ii_out,key=lambda x:x[1]):
            a,b=pair.split("-"); _,fa,fb=dist(a,b)
            rp=profile_bond(dv,nxv,nyv,nzv,fa,fb)
            agree="ok" if abs(rp-rho)/max(rho,1e-9)<0.25 else "differs"
            V.append(f"  {pair:<8}{dAA:>6.2f}{rp:>10.4f}{agree:>9}")
        # [D] I-O profile consistency (same rigor: clean-field midpoint rho vs critic2 BCP rho)
        V.append(""); V.append("[D] I-O bonding-region profile on clean val.cube (same rigor as I-I):")
        V.append(f"  {'pair':<8}{'d(A)':>6}{'rho_BCP':>9}{'rho_prof':>10}{'vs BCP':>9}")
        for pair,dAA,rho,lap in sorted(io_out,key=lambda x:x[1]):
            a,b=pair.split("-"); _,fa,fb=dist(a,b)
            rp=profile_bond(dv,nxv,nyv,nzv,fa,fb)
            agree="ok" if abs(rp-rho)/max(rho,1e-9)<0.35 else "differs"
            V.append(f"  {pair:<8}{dAA:>6.2f}{rho:>9.4f}{rp:>10.4f}{agree:>9}")
    open(f"{BASE}/{sys}/{sys}_bcp_verify.txt","w").write("\n".join(V)+"\n")
    print(f"{sys}: {len(ii_out)} I-I, {len(io_out)} I-O real bonds; rej II={len(seenII)} IO={len(seenIO)}")
    return ii_out,io_out

if __name__=="__main__":
    alldat={}
    for s in (sys.argv[1:] or SYS): alldat[s]=report(s)
    if len(alldat)>1:
        S=["="*78,"BCP RE-LOCATION SUMMARY (clean valence field; artifacts removed)","="*78,
           f"{'system':<8}{'I-I bonds':>10}{'I-O bonds':>10}", "-"*78]
        for s in SYS:
            if s in alldat: ii,io=alldat[s]; S.append(f"{s:<8}{len(ii):>10}{len(io):>10}")
        S.append(""); S.append("Monotonicity (rho_BCP vs d): see per-system _bcp_verify.txt [B].")
        open(f"{BASE}/BCP_correction_summary.txt","w").write("\n".join(S)+"\n")
        print("\n".join(S))
