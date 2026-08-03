import re,numpy as np
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
A="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"

for s in ["Yvc","Gdvc"]:
    inp=open(f"{A}/rd_poly_{s}.in").read()
    pos=re.search(r"ATOMIC_POSITIONS crystal\n(.*?)\nK_POINTS",inp,re.S).group(1).strip().split("\n")
    seed_I=np.array([[float(x) for x in l.split()[1:4]] for l in pos if l.startswith("I ")])
    cp=re.search(r"CELL_PARAMETERS angstrom\n(.*?)\nATOMIC",inp,re.S).group(1).split("\n")
    cell=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)])

    # read all lines, find I-lines, take last 8
    lines=open(f"{T}/relax442_{s}.out",errors="ignore").readlines()
    i_lines=[l for l in lines if l.lstrip().startswith("I ") and len(l.split())>=4]
    final_I=np.array([[float(x) for x in l.split()[1:4]] for l in i_lines[-8:]])

    df=final_I-seed_I; df=df-np.round(df); dc=df@cell
    dists=np.linalg.norm(dc,axis=1)
    print(f"=== {s}: Γ-relaxed → 4×4×2-relaxed ===")
    print(f"RMSD={np.sqrt((dists**2).mean()):.3f} Å | max={dists.max():.3f} Å (I{np.argmax(dists)+1})")
    for i in range(8): print(f"  I{i+1}: {dists[i]:.3f} Å")
    print(f"\n  I-I distances (Å): Γ-seed → 4×4×2-relaxed")
    for i,j in [(0,4),(0,1),(1,2),(1,3),(2,3),(4,5),(5,6),(5,7),(6,7)]:
        def dij(fr): d=fr[i]-fr[j]; d=d-np.round(d); return np.linalg.norm(d@cell)
        ds,df_=dij(seed_I),dij(final_I)
        print(f"  I{i+1}-I{j+1}: {ds:.3f} → {df_:.3f} (Δ={df_-ds:+.3f} Å)")
    print()
