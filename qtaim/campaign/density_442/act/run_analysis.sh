#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
echo "=== 2x cubes via Fourier interpolation ==="
python3 <<'PYEOF'
import numpy as np
A='./'
def read_cube(p):
    with open(p) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        dims=[int(f.readline().split()[0]) for _ in range(3)]
        for _ in range(nat): f.readline()
        return dims, np.loadtxt(f).reshape(dims)
def interp2x(x):
    F=np.fft.fftn(x); nx,ny,nz=x.shape
    P=np.zeros((2*nx,2*ny,2*nz),dtype=complex)
    h=[slice(0,(n+1)//2) for n in (nx,ny,nz)]; t=[slice(-(n//2),None) for n in (nx,ny,nz)]
    for a in (h,t):
        for b in (h,t):
            for c in (h,t):
                P[a[0],b[1],c[2]]=F[a[0],b[1],c[2]]
    return np.real(np.fft.ifftn(P))*8
def write_cube_like(src, data, out):
    with open(src) as f: n=int(f.readlines()[2].split()[0])
    with open(src) as f, open(out,'w') as g:
        for i in range(6+n): g.write(f.readline())
        g.write("\n".join(f"{v:.5e}" for v in np.ravel(data))+"\n")
for name in ('rho17_act_poly_lS2','val_act_poly_lS2','empty_val_act_poly_lS2','i2only_val_act_poly_lS2'):
    d,x=read_cube(A+name+'.cube'); assert d==[160,180,300],(name,d)
    x2=interp2x(x); write_cube_like(A+name+'.cube',x2,A+name+'_2x.cube')
    print(f'{name}_2x: {x2.shape} sum*8-consistent: {abs(x2.sum()-8*x.sum())<1e-2*x.sum()}')
_,vl=read_cube(A+'val_act_poly_lS2_2x.cube'); _,ve=read_cube(A+'empty_val_act_poly_lS2_2x.cube'); _,vi=read_cube(A+'i2only_val_act_poly_lS2_2x.cube')
write_cube_like(A+'val_act_poly_lS2_2x.cube',vl-ve-vi,A+'deltarho_2x_act_poly_lS2.cube')
print('deltarho_2x written')
PYEOF
echo "=== critic2 runs ==="
C2=/workspace/dimitry/calc/qe/i2_mof/tools/critic2
for j in bcp_act_poly_lS2 ct_act_poly_lS2 bcp_act_mol_rd ct_act_mol_rd ct_2x_act_poly_lS2; do
  echo "RUN $j $(date +%H:%M:%S)"
  $C2 $j.cri > $j.out 2>&1 && echo "OK $j" || echo "FAIL $j"
done
echo "ANALYSIS ALL DONE"
