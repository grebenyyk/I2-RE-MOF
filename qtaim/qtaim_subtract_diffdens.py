import numpy as np, sys
LD="13Y_opt/ae_rho.cube"; EM="13Y_empty/ae_rho.cube"; IO="13Y_i2only/ae_rho.cube"; OUT="13Y_diffdens.cube"
def read_cube(fn):
    with open(fn) as f:
        lines=f.readlines()
    nat=int(lines[2].split()[0])
    nx,ny,nz=[int(lines[i].split()[0]) for i in (3,4,5)]
    hdr=lines[:6+nat]   # 2 title + origin + 3 grid + nat atom lines
    data=np.array([float(x) for l in lines[6+nat:] for x in l.split()])
    data=data.reshape((nx,ny,nz))
    return hdr,data,nx,ny,nz
hL,dL,nx,ny,nz=read_cube(LD)
_,dE,_,_,_=read_cube(EM)
_,dI,_,_,_=read_cube(IO)
assert dE.shape==dL.shape==dI.shape, "grid mismatch!"
dd=dL-dE-dI
# write diff cube using loaded header (atoms shown)
with open(OUT,"w") as o:
    o.writelines(hL)
    flat=dd.ravel(order='F')  # cube stores z-fastest in Fortran order... actually QE cube is C-contiguous with z fastest: reshape default C; ravel F may be wrong. Match by using same order as parsed (parsed flat->reshape default C). Re-emit in C order same as input.
    # Re-derive: input parsed as reshape((nx,ny,nz)) default C; data lines fill [i][j][k] with k fastest. Re-emit identically:
    flat=dd.ravel(order='C')
    for i in range(0,len(flat),6):
        o.write("".join(f"{v:13.5E}" for v in flat[i:i+6])+"\n")
# stats
print(f"δρ cube written: {OUT}  grid {nx}x{ny}x{nz}")
print(f"δρ min={dd.min():.4e}  max={dd.max():.4e}  e/bohr³")
vol=np.dot(np.array([float(x) for x in hL[3].split()[1:4]]),np.cross( np.array([float(x) for x in hL[4].split()[1:4]]), np.array([float(x) for x in hL[5].split()[1:4]]) ))
# cell vol from vectors; but careful: hL[3] is nx + 3 vec comps. vol = |a_vec . (b x c)| where a_vec = (line4 vec)*? Actually grid vec line = (a_cell/nx). So cell vol = (nx*ny*nz)*|det of grid-step matrix|
ga=np.array([float(x) for x in hL[3].split()[1:4]]); gb=np.array([float(x) for x in hL[4].split()[1:4]]); gc=np.array([float(x) for x in hL[5].split()[1:4]])
cellvol=abs(np.dot(ga*nx,np.cross(gb*ny,gc*nz)))
integral=dd.sum()*abs(np.dot(ga,np.cross(gb,gc)))  # voxel vol × sum
print(f"voxel vol={abs(np.dot(ga,np.cross(gb,gc))):.4e} bohr³, cell vol={cellvol:.1f} bohr³")
print(f"∫δρ dV = {integral:+.3f} e  (should be ~0: charge conserved)")
# net charge transferred framework->I: integrate δρ over the I region vs framework region
print(f"\nload into VESTA with isosurface ±0.005 e/bohr³ (blue=accumulation on I, red=depletion on framework)")
