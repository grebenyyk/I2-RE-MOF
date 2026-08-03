import numpy as np, os
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"; T=C+"/kmesh_test"
def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline(); nat=int(f.readline().split()[0])
        gx=f.readline().split(); gy=f.readline().split(); gz=f.readline().split()
        nx,ny,nz=int(gx[0]),int(gy[0]),int(gz[0])
        for _ in range(nat): f.readline()
        return np.fromstring(" ".join(f.readlines()),sep=" ").reshape(nx,ny,nz),nx,ny,nz
def rms(a,b): return np.sqrt(((a-b)**2).mean())
print("================ PAIRWISE k-mesh convergence (Y_deep) ================")
for field,globref in [("DENSITY rho17","rho17"),("ELF","elf")]:
    cubes={}
    refpath=C+"/elf.cube" if field.startswith("ELF") else C+"/rho17.cube"
    g,_,_,_=read_cube(refpath); cubes["G"]=g
    for t in ["k221","k442","k663"]:
        p=f"{T}/{globref}_{t}.cube" if field.startswith("DENSITY") else f"{T}/elf_{t}.cube"
        if os.path.exists(p) and os.path.getsize(p)>1e7:
            cubes[t],_,_,_=read_cube(p)
    keys=["G","k221","k442","k663"]; keys=[k for k in keys if k in cubes]
    print(f"\n--- {field}: RMS difference, all pairs (ELF in [0,1], rho in e/bohr3) ---")
    rr=np.sqrt((cubes["G"]**2).mean())
    print(f"    {'vs G':>16}{'G-k221':>10}{'G-k442':>10}{'G-k663':>10}{'k221-k442':>11}{'k442-k663':>11}")
    vals={}
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            vals[f"{a}-{b}"]=rms(cubes[a],cubes[b])
    def gpair(x): return vals.get(x,float('nan'))
    if field.startswith("ELF"):
        print(f"    {'':>16}{gpair('G-k221'):>10.4f}{gpair('G-k442'):>10.4f}{gpair('G-k663'):>10.4f}{gpair('k221-k442'):>11.4f}{gpair('k442-k663'):>11.4f}")
        print(f"    (ELF: k442-k663 = {gpair('k442-k663'):.4f} => {'CONVERGED' if gpair('k442-k663')<0.005 else 'still moving'})")
    else:
        print(f"    {'':>16}{gpair('G-k221')/rr*100:>9.3f}%{gpair('G-k442')/rr*100:>9.3f}%{gpair('G-k663')/rr*100:>9.3f}%{gpair('k221-k442')/rr*100:>10.3f}%{gpair('k442-k663')/rr*100:>10.3f}%")
        print(f"    (rho: k442-k663 = {gpair('k442-k663')/rr*100:.3f}% of rho_rms => {'CONVERGED' if gpair('k442-k663')/rr*100<0.1 else 'still moving'})")
