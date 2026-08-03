#!/usr/bin/env python3
"""Build iodine-only relax inputs for the Y_act (activated) pore with the polyiodide
chain rotated +120 deg and +135 deg about the pore (c) axis to test whether it
relaxes to the st-like Y-anchored motif. Framework frozen (if_pos 0 0 0), I free
(1 1 1). Settings: production (50/400, 4x4x2, PBE-D3, conv 1e-6, etot 1e-5, forc 5e-3, bfgs)."""
import re,numpy as np,os
A="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_act_rot"; os.makedirs(T,exist_ok=True)
raw=open(f"{A}/rd_poly_Yact.in").read()
cp=re.search(r"CELL_PARAMETERS.*?\n(.*?\n.*?\n.*?)\nATOMIC",raw,re.S).group(1).split("\n")
cell=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)]); inv=np.linalg.inv(cell)
pos=re.search(r"(ATOMIC_POSITIONS crystal\n.*?\n)K_POINTS",raw,re.S).group(1)
atoms=[(l.split()[0],np.array([float(x) for x in l.split()[1:4]])) for l in pos.strip().split("\n")[1:]]
def cart(f): return f[0]*cell[0]+f[1]*cell[1]+f[2]*cell[2]
def rot_about(v,ax,t):
    c=np.cos(t/2); s=np.sin(t/2); qv=s*ax
    return v+2*np.cross(qv,np.cross(qv,v)+c*v)
I_idx=[i for i,(s,_) in enumerate(atoms) if s=="I"]
Icart=np.array([cart(atoms[i][1]) for i in I_idx]); cen=Icart.mean(0)
caxis=cell[2]/np.linalg.norm(cell[2])
template="""&CONTROL
  calculation = 'relax'
  prefix = 'Y_act_rot{tag}'
  outdir = './tmp_rot{tag}'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-5
  forc_conv_thr = 5.0d-3
/
&SYSTEM
  ibrav = 0
  nat = 80
  ntyp = 5
  ecutwfc = 50
  ecutrho = 400
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.01
  nspin = 1
  vdw_corr = 'grimme-d3'
  dftd3_version = 4
/
&ELECTRONS
  conv_thr = 1.0d-6
  diagonalization = 'david'
  mixing_mode = 'local-TF'
  mixing_beta = 0.2
  mixing_ndim = 16
  electron_maxstep = 300
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
C       12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF
H        1.008  H.pbe-kjpaw.UPF
O       15.999  O.pbe-kjpaw.UPF
Y       88.906  Y.pbe-spn-kjpaw_psl.1.0.0.UPF
I      126.904  I.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
{celllines}
ATOMIC_POSITIONS crystal
{poslines}
K_POINTS automatic
  4 4 2 0 0 0
"""
celllines="\n".join("  ".join(f"{x:.10f}" for x in cell[i]) for i in range(3))
for deg in [120,135]:
    th=np.radians(deg)
    Irot=np.array([rot_about(p-cen,caxis,th)+cen for p in Icart])
    Irot_frac=(Irot@inv)
    lines=[]
    for i,(s,f) in enumerate(atoms):
        if i in I_idx:
            j=I_idx.index(i); ff=Irot_frac[j]; flag="1 1 1"
        else: ff=f; flag="0 0 0"
        lines.append(f"{s} {ff[0]:.12f} {ff[1]:.12f} {ff[2]:.12f} {flag}")
    inp=template.format(tag=deg,celllines=celllines,poslines="\n".join(lines))
    open(f"{T}/relax_{deg}.in","w").write(inp)
    print(f"wrote relax_{deg}.in (chain rotated {deg}deg about c-axis; framework frozen, I free)")
