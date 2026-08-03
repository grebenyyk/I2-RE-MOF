import math
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test"
def bcps(path):
    out=[]
    for line in open(path,errors="ignore"):
        if "(3,-1)" not in line or "bond" not in line: continue
        f=line.split()
        try: x,y,z,rho=float(f[4]),float(f[5]),float(f[6]),float(f[9])
        except: continue
        out.append((x,y,z,rho))
    return out
def mind(p,q):
    s=0.0
    for i in range(3):
        d=abs(p[i]-q[i])%1.0
        if d>0.5: d=1.0-d
        s+=d*d
    return math.sqrt(s)
sets={t:bcps(f"{T}/bcp_{t}_full.out") for t in ["gamma","k221","k442"]}
for t,s in sets.items(): print(f"{t}: {len(s)} bond CPs found")
def matchfrac(A,B,tol):
    n=sum(1 for a in A if min(mind(a,b) for b in B)<tol)
    return n,len(A)
for tol in [0.01,0.005]:
    print(f"\n--- tol={tol} fractional (~{tol*9:.2f} Å) ---")
    for a,b in [("gamma","k221"),("gamma","k442"),("k221","k442")]:
        nA,_=matchfrac(sets[a],sets[b],tol); nB,_=matchfrac(sets[b],sets[a],tol)
        print(f"  {a}->{b}: {nA}/{len(sets[a])} {a} BCPs have a {b} counterpart;  {b}->{a}: {nB}/{len(sets[b])} {b} BCPs have a {a} counterpart")
# I-I BCPs specifically: do the curated chain bonds appear on all meshes?
print("\n=== I-I chain BCPs (rho 0.045-0.063) count per mesh ===")
for t,s in sets.items():
    ii=[p for p in s if 0.045<=p[3]<=0.063]
    print(f"  {t}: {len(ii)} BCPs in I-I chain range")
