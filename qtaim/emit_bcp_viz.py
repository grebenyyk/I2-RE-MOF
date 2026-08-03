import re, numpy as np
SCF="/workspace/dimitry/calc/qe/i2_mof/qtaim/13Y_opt/scf.in"
BCP="/workspace/dimitry/calc/qe/i2_mof/qtaim/13Y_opt/bcp17.out"
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/13Y_opt/13Y_with_BCPs.vasp"
BOHR=0.5291772109
# --- trusted cell + atoms from scf.in ---
s=open(SCF).read()
cell_Å=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l))) for l in s.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
apos=[]
for l in s.split("ATOMIC_POSITIONS")[1].splitlines():
    m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",l)
    if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
# --- CPs from critic2 crystallographic coords (true fractional) ---
t=open(BCP,errors="ignore").read()
parts=re.split(r"\+ Critical point no\. (\d+)",t)
cps=[]
for i in range(1,len(parts),2):
    b=parts[i+1]
    typ=re.search(r"Type :\s*(\S+(?:,\S+)?)",b)
    crys=re.search(r"Crystallographic coordinates:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",b)
    f=re.search(r"Field value \(f\):\s*(-?\d+\.\d+E[+-]\d+)",b)
    lap=re.search(r"Laplacian \(del2 f\):\s*(-?\d+\.\d+E[+-]\d+)",b)
    if typ and crys and f and lap:
        cps.append((int(parts[i]),typ.group(1),np.array([float(crys.group(1)),float(crys.group(2)),float(crys.group(3))]),float(f.group(1)),float(lap.group(1))))
def wrap(f): return f-np.floor(f)
def mic_d(fc,fa):  # fractional mic distance in Å
    d=wrap(fc-fa+0.5)-0.5
    return np.linalg.norm(cell_Å.T@d)  # cell_Å rows are a,b,c
bcp=[c for c in cps if c[1]=="(3,-1)"]
# label each BCP by 2 nearest atoms (true frame)
I_bcp=[]
for n,typ,fc,rho,lap in bcp:
    ds=sorted((mic_d(fc,a[1]),k) for k,a in enumerate(apos))
    lab="".join(sorted([apos[ds[0][1]][0],apos[ds[1][1]][0]]))
    if "I" not in lab: continue
    minr=min(ds[0][0],ds[1][0]); blen=ds[0][0]+ds[1][0]
    if minr<0.55*BOHR*1.0: minrA=minr; 
    # 0.55 Å threshold for artifact
    art = minr < 0.55
    I_bcp.append((n,lab,blen,minr,rho,lap,fc,art))
DUMMY={"II":"X","IO":"F","IY":"Cl","CI":"S","HI":"P"}
phys=[r for r in I_bcp if not r[7]]
bykind={}
for n,lab,blen,minr,rho,lap,fc,art in phys:
    bykind.setdefault(lab,[]).append(fc)
# --- write POSCAR ---
spec=["C","H","O","Y","I"]; order=[o for o in ["II","IO","IY","CI","HI"] if o in bykind]
cnt={s:sum(1 for a in apos if a[0]==s) for s in spec}
L=["13Y uniform_chain (strained) + I-related BCPs (valence dens plot_num=17) [FIXED coords]","1.0"]
for r in cell_Å: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
L.append("   "+"   ".join(f"{s:>4s}" for s in spec)+"   "+"   ".join(f"{DUMMY[o]:>4s}" for o in order))
L.append("   "+"   ".join(f"{cnt[s]:4d}" for s in spec)+"   "+"   ".join(f"{len(bykind[o]):4d}" for o in order))
L.append("Direct")
for s in spec:
    for a in apos:
        if a[0]==s:
            w=wrap(a[1]); L.append(f"   {w[0]: .12f}   {w[1]: .12f}   {w[2]: .12f}")
for o in order:
    for fc in bykind[o]:
        w=wrap(fc); L.append(f"   {w[0]: .12f}   {w[1]: .12f}   {w[2]: .12f}")
open(OUT,"w").write("\n".join(L)+"\n")
# sanity: verify a couple of framework atoms match scf.in after wrap, and BCP placement
print("cell Å:\n",cell_Å)
print(f"atoms: C{cnt['C']} H{cnt['H']} O{cnt['O']} Y{cnt['Y']} I{cnt['I']} = {sum(cnt.values())}")
print(f"I-related physical BCPs kept: "+" ".join(f"{DUMMY[o]}={o}({len(bykind[o])})" for o in order)+f" = {sum(len(v) for v in bykind.values())}")
print(f"\nfirst 3 wrapped C atoms (should match scf.in 0.386/0.157/0.068...):")
for a in [x for x in apos if x[0]=='C'][:3]: print("  ",np.round(wrap(a[1]),4))
# verify I-I BCP ncp145 sits between 2 I at ~2.95A
print(f"\nflux check — I-I BCP ncp145: nearest 2 atoms distance sum:")
n=145;c=[c for c in cps if c[0]==n][0]; ds=sorted((mic_d(c[2],a[1]),k,apos[k][0]) for k,a in enumerate(apos))[:4]
for d,k,s in ds: print(f"   {s}({k}): {d:.3f} Å")
# --- corrected table ---
print(f"\n{'='*72}\nCORRECTED I-related BCP table (true fractional frame):\n{'ncp':>4} {'bond':<5} {'d(Å)':>6} {'minR(Å)':>7} {'ρ(e/bohr³)':>11} {'∇²ρ':>10}  status")
print("-"*60)
from collections import defaultdict
g=defaultdict(list)
for n,lab,blen,minr,rho,lap,fc,art in sorted(I_bcp):
    st="artifact" if art else "physical"
    print(f"{n:>4} {lab:<5} {blen:>6.2f} {minr:>7.2f} {rho:>11.3e} {lap:>10.2e}  {st}")
    if not art: g[lab].append((blen,rho,lap))
print("\nphysical BCPs by type:")
for lab,rows in sorted(g.items()):
    print(f"  {lab} ({len(rows)}): "+" | ".join(f"d={r[0]:.2f} ρ={r[1]:.2e} ∇²ρ={r[2]:.2e}" for r in rows))
