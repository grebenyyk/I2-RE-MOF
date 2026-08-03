import numpy as np
def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        gx=f.readline().split(); gy=f.readline().split(); gz=f.readline().split()
        nx,ny,nz=int(gx[0]),int(gy[0]),int(gz[0])
        for _ in range(nat): f.readline()
        data=np.fromstring(" ".join(f.readlines()),sep=" ").reshape(nx,ny,nz)
    return data
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
eD=read_cube(f"{C}/Y_deep/elf.cube")
eA=read_cube(f"{C}/Y_actc/elf.cube")
def stats(e,lbl):
    flat=e.ravel()
    h,_=np.histogram(flat,bins=np.linspace(0,1,21))
    return flat,h
fD,hD=stats(eD,"Y_deep")
fA,hA=stats(eA,"Y_actc")
print("=== ELF value distribution: Y_deep (narrow band) vs Y_actc (wider band) ===")
print(f"{'':32}{'Y_deep':>12}{'Y_actc':>12}")
print(f"{'mean ELF':32}{fD.mean():>12.4f}{fA.mean():>12.4f}")
print(f"{'median ELF':32}{np.median(fD):>12.4f}{np.median(fA):>12.4f}")
# metallic delocalization signature: ELF in [0.4,0.6] (the "metallic" channel, neither localized nor avoided)
metD=np.mean((fD>=0.4)&(fD<=0.6)); metA=np.mean((fA>=0.4)&(fA<=0.6))
locD=np.mean(fD>0.7); locA=np.mean(fA>0.7)
print(f"{'frac ELF in [0.4,0.6] (metallic)':32}{metD:>12.4f}{metA:>12.4f}")
print(f"{'frac ELF >0.7 (localized)':32}{locD:>12.4f}{locA:>12.4f}")
print("\n=== histogram (fraction of voxels per 0.05 ELF bin) ===")
print(f"{'ELF bin':>8}{'Y_deep':>10}{'Y_actc':>10}{'deep-actc':>10}")
c=(hD[0]+hA[0]);
for i in range(20):
    d=hD[i]/hD.sum(); a=hA[i]/hA.sum()
    print(f"{i/20:>8.2f}{d:>10.4f}{a:>10.4f}{d-a:>+10.4f}")
print("\n(Y_deep has MORE voxels in the metallic [0.4-0.6] channel => more delocalized/metallic ELF,")
print(" consistent with its narrower band. Direct cube subtraction impossible: different cells/grids.)")
