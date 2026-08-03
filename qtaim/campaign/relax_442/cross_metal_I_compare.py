import re,numpy as np
A="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"

def parse(fn):
    raw=open(fn).read()
    cp=re.search(r"CELL_PARAMETERS angstrom\n(.*?)\nATOMIC",raw,re.S).group(1).split("\n")
    cell=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)])
    pos=re.search(r"ATOMIC_POSITIONS crystal\n(.*?)\nK_POINTS",raw,re.S).group(1).strip().split("\n")
    I=np.array([[float(x) for x in l.split()[1:4]] for l in pos if l.startswith("I ")])
    return cell,I

# Γ-relaxed geometries (production rd_poly)
metals=["Y","Gd","Tb"]
cells={}; I_gamma={}
for m in metals:
    cell,I=parse(f"{A}/rd_poly_{m}act.in")
    cells[m]=cell; I_gamma[m]=I

print("="*70)
print("Q1: Lattice comparison (act variant)")
print("="*70)
for m in metals:
    a,b,c=np.linalg.norm(cells[m][0]),np.linalg.norm(cells[m][1]),np.linalg.norm(cells[m][2])
    print(f"  {m}act: a={a:.3f} b={b:.3f} c={c:.3f} Å")
print(f"  Max cell-vector diff Y↔Gd: {np.abs(cells['Y']-cells['Gd']).max():.3f} Å")
print(f"  Max cell-vector diff Y↔Tb: {np.abs(cells['Y']-cells['Tb']).max():.3f} Å")

print(f"\n{'='*70}")
print("Q2: Γ-relaxed I positions across metals (how different are they?)")
print("="*70)
# Compare I positions in Cartesian (each in its own cell)
for m in metals:
    Ic=I_gamma[m]@cells[m]
    print(f"\n  {m}act I atoms (Cartesian, own cell):")
    for i in range(8):
        print(f"    I{i+1}: [{Ic[i,0]:.3f}, {Ic[i,1]:.3f}, {Ic[i,2]:.3f}]")

# Cross-metal I position comparison: transfer Y frac coords to Gd/Tb cells
print(f"\n{'='*70}")
print("Q1: Y act frac coords → Gd/Tb cells: how far from native Γ-relaxed I?")
print("="*70)
for target in ["Gd","Tb"]:
    # Y fractional I coords placed in target cell (Cartesian)
    Y_in_target = I_gamma["Y"] @ cells[target]
    # Target native I coords (Cartesian)
    native = I_gamma[target] @ cells[target]
    # displacement (min-image in target cell)
    inv=np.linalg.inv(cells[target])
    df=(I_gamma["Y"]-I_gamma[target])@inv; df=df-np.round(df)
    dc=df@cells[target]
    dists=np.linalg.norm(dc,axis=1)
    print(f"\n  Y→{target} act:")
    print(f"    RMSD: {np.sqrt((dists**2).mean()):.3f} Å | max: {dists.max():.3f} Å (I{np.argmax(dists)+1})")
    for i in range(8):
        print(f"    I{i+1}: Δ={dists[i]:.3f} Å")

# Also: Y act 4×4×2-relaxed (damp) I coords → Gd/Tb cells
print(f"\n{'='*70}")
print("Q1b: Y act 4×4×2-RELAXED frac coords → Gd/Tb cells")
print("="*70)
# extract Y act damp final I positions
raw_damp=open(f"{T}/relax442_Yact_damp.out",errors="ignore").read()
parts=raw_damp.split("ATOMIC_POSITIONS")
last_I=[]
for part in reversed(parts[1:]):
    lines=part.strip().split("\n")[1:]
    test=[]
    for l in lines:
        p=l.split()
        if len(p)>=4 and p[0]=="I": test.append([float(p[1]),float(p[2]),float(p[3])])
    if len(test)==8: last_I=test; break
I_yact_relaxed=np.array(last_I)
print(f"  Y act 4×4×2-relaxed I frac coords extracted ({len(I_yact_relaxed)} atoms)")
for target in ["Gd","Tb"]:
    df=(I_yact_relaxed-I_gamma[target])@np.linalg.inv(cells[target]); df=df-np.round(df)
    dc=df@cells[target]
    dists=np.linalg.norm(dc,axis=1)
    print(f"\n  Y-relaxed → {target} act:")
    print(f"    RMSD from {target} Γ-relaxed: {np.sqrt((dists**2).mean()):.3f} Å | max: {dists.max():.3f} Å")

# Q2: detailed comparison of Γ-relaxed I positions (Cartesian) Y vs Gd vs Tb
print(f"\n{'='*70}")
print("Q2: Pairwise Γ-relaxed I comparison (Cartesian, each in own cell)")
print("="*70)
for a,b in [("Y","Gd"),("Y","Tb"),("Gd","Tb")]:
    Ia=I_gamma[a]@cells[a]; Ib=I_gamma[b]@cells[b]
    # just compare position-by-position (same atom ordering)
    dists=np.linalg.norm(Ia-Ib,axis=1)
    print(f"  {a}↔{b}: RMSD={np.sqrt((dists**2).mean()):.3f} Å | max={dists.max():.3f} Å (I{np.argmax(dists)+1})")
    for i in range(8):
        print(f"    I{i+1}: {dists[i]:.3f} Å")
