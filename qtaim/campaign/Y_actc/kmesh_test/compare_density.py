import numpy as np, os
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_actc/kmesh_test"
def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        gx=f.readline().split(); gy=f.readline().split(); gz=f.readline().split()
        nx,ny,nz=int(gx[0]),int(gy[0]),int(gz[0])
        for _ in range(nat): f.readline()
        data=np.fromstring(" ".join(f.readlines()),sep=" ")
    return data.reshape((nx,ny,nz))
print("reading cubes (160x180x300)...",flush=True)
g=read_cube(f"{T}/rho17_gamma.cube")
rrms=np.sqrt((g**2).mean()); rmax=g.max()
print(f"  gamma: shape {g.shape} rho_rms={rrms:.5e} rho_max={rmax:.4f}")
for tag in ["k221","k442","k663"]:
    p=f"{T}/rho17_{tag}.cube"
    if not (os.path.exists(p) and os.path.getsize(p)>1e7):
        print(f"  gamma vs {tag}: (cube not ready)"); continue
    d=read_cube(p); diff=g-d
    rms=np.sqrt((diff**2).mean()); mx=np.abs(diff).max()
    print(f"  gamma vs {tag}: RMS_drho={rms:.4e} = {rms/rrms*100:.4f}% of rho_rms | max|drho|={mx:.4e} = {mx/rmax*100:.4f}% of rho_max | bias={diff.mean():+.3e}")
print("\n<0.1% relative RMS => density k-converged at Gamma.")
