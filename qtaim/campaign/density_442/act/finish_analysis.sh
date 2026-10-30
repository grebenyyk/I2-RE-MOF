#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/act
python3 <<'PYEOF'
import numpy as np
A='./'
# 1. fix 2x headers (dims + voxel vectors) from pristine 1x grid lines
def header_1x(p):
    L=open(p).readlines()
    nat=int(L[2].split()[0]); grid=L[3:6]
    return L,nat,grid
srcL,nat,grid1=header_1x('rho17_act_poly_lS2.cube')
grid2=[]
for g in grid1:
    t=g.split(); n=int(t[0]); v=[float(x) for x in t[1:4]]
    grid2.append(f"{2*n} "+" ".join(f"{x/2:.9f}" for x in v)+"\n")
for name,nat in (('rho17_act_poly_lS2',80),('val_act_poly_lS2',80),('empty_val_act_poly_lS2',72),('i2only_val_act_poly_lS2',8)):
    L=open(A+name+'_2x.cube').readlines()
    data_idx=6+nat
    out=L[:3]+grid2+L[6:data_idx]+L[data_idx:]
    open(A+name+'_2x.cube','w').writelines(out)
    print(f'fixed header: {name}_2x (nat={nat})')
# 2. deltarho_2x
def read_cube(p):
    with open(p) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        dims=[int(f.readline().split()[0]) for _ in range(3)]
        for _ in range(3): f.readline()
        for _ in range(nat): f.readline()
        return dims, np.loadtxt(f).reshape(dims)
_,vl=read_cube(A+'val_act_poly_lS2_2x.cube'); _,ve=read_cube(A+'empty_val_act_poly_lS2_2x.cube'); _,vi=read_cube(A+'i2only_val_act_poly_lS2_2x.cube')
out=srcL[:3]+grid2+srcL[6:6+80]
open(A+'deltarho_2x_act_poly_lS2.cube','w').writelines(out)
with open(A+'deltarho_2x_act_poly_lS2.cube','a') as g:
    g.write("\n".join(f"{v:.5e}" for v in np.ravel(vl-ve-vi))+"\n")
d,x=read_cube(A+'deltarho_2x_act_poly_lS2.cube')
print('deltarho_2x grid:',d)
# 3. critic2 inputs
for tag in ('act_poly_lS2','act_mol_rd'):
    open(A+f'bcp_{tag}.cri','w').write(f"crystal rho17_{tag}.cube\nload rho17_{tag}.cube id rho\nauto\n")
    open(A+f'ct_{tag}.cri','w').write(f"crystal rho17_{tag}.cube\nload rho17_{tag}.cube id rho\nload deltarho_{tag}.cube id dr\nintegrable dr\nbader\n")
open(A+'ct_2x_act_poly_lS2.cri','w').write("crystal rho17_2x_act_poly_lS2.cube\nload rho17_2x_act_poly_lS2.cube id rho\nload deltarho_2x_act_poly_lS2.cube id dr\nintegrable dr\nbader\n")
print('cri files written')
PYEOF
C2=/workspace/dimitry/calc/qe/i2_mof/tools/critic2_run.sh
for j in bcp_act_poly_lS2 ct_act_poly_lS2 bcp_act_mol_rd ct_act_mol_rd ct_2x_act_poly_lS2; do
  echo "RUN $j $(date +%H:%M:%S)"
  $C2 $j.cri > $j.out 2>&1 && echo "OK $j" || echo "FAIL $j"
done
echo "ANALYSIS ALL DONE"
