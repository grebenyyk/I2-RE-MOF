#!/usr/bin/env python3
"""Build SCF + pp + critic2 inputs for the 3 DEEP-MINIMUM systems (sandwiched I2 +
uniform chain). Same settings as the fixed-cell campaign: 60/800/Gamma, degauss 0.01,
mv, local-TF mixing 0.2 ndim 16, david, conv 1e-8. Loaded(80) + empty(72, frozen fw)
+ i2only(8), all at the SAME deep cell/positions (isolates electronic reorganization)."""
import os, re
import numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"
OUT=f"{BASE}/qtaim/campaign"
# deep geometry source: vc-relax output (Gd/Tb seeded, Y chain-start that found deep)
SYS=[("Gd_deep",f"{OUT}/Gd_act/vc-relax_seedY.out","Gd"),
     ("Tb_deep",f"{OUT}/Tb_act/vc-relax_seedY.out","Tb"),
     ("Y_deep", f"{OUT}/Y_actc/vc-relax_t1.out","Y")]
PPMAP={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
       "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
       "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Tb":158.925,"Gd":157.25,"I":126.904}
def read_deep(out):
    t=open(out,errors="ignore").read()
    cell=None
    for blk in re.findall(r"CELL_PARAMETERS\s*\(?[^\n]*\n((?:[^\n]+\n){3})",t)[::-1]:
        n=re.findall(r"-?\d+\.\d+",blk)
        if len(n)>=9: cell=[[float(n[i*3+j]) for j in range(3)] for i in range(3)]; break
    b=re.split(r"ATOMIC_POSITIONS[^\n]*\n",t); atoms=[]
    for ln in b[-1].split("\n"):
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",ln)
        if m and m.group(1) in MASS and len(atoms)<80:
            atoms.append((m.group(1),np.array([float(m.group(i)) for i in (2,3,4)])))
    return cell,atoms
def write_scf(path,prefix,cell,atoms,species):
    ntyp=len(set(species)); outdir='./tmp_'+prefix.split('_')[-1]
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n  tstress = .true.\n/\n",
       "&SYSTEM\n  ibrav = 0\n",f"  nat = {len(atoms)}\n",f"  ntyp = {ntyp}\n",
       "  ecutwfc = 60\n  ecutrho = 800\n  occupations = 'smearing'\n  smearing = 'mv'\n",
       "  degauss = 0.01\n  nspin = 1\n  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/\n",
       "&ELECTRONS\n  conv_thr = 1.0d-8\n  diagonalization = 'david'\n  diago_david_ndim = 2\n",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.2\n  mixing_ndim = 16\n  electron_maxstep = 300\n/\n",
       "ATOMIC_SPECIES\n"]
    for s in species: L.append(f"{s:<3s} {MASS[s]:8.3f}  {PPMAP[s]}\n")
    L.append("CELL_PARAMETERS angstrom\n")
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}\n")
    L.append("ATOMIC_POSITIONS crystal\n")
    for s,fr in atoms: L.append(f"{s:<3s} {fr[0]: .12f} {fr[1]: .12f} {fr[2]: .12f}\n")
    L.append("K_POINTS automatic\n  1 1 1 0 0 0\n")
    open(path,"w").write("".join(L))
def write_pp(path,prefix,plotnum,outfile):
    open(path,"w").write(f"&INPUTPP\n  prefix = '{prefix}'\n  outdir = './tmp'\n  plot_num = {plotnum}\n  filplot = 'fld'\n/\n&PLOT\n  nfile = 1\n  filepp(1) = 'fld'\n  weight(1) = 1.0\n  iflag = 3\n  output_format = 6\n  fileout = '{outfile}'\n/\n")
for name,src,metal in SYS:
    d=f"{OUT}/{name}"; os.makedirs(d,exist_ok=True)
    cell,atoms=read_deep(src)
    species=["C","H","O",metal,"I"]
    write_scf(f"{d}/scf_loaded.in",f"{name}_loaded",cell,atoms,species)
    fw=[(s,fr) for s,fr in atoms if s!="I"]; io=[(s,fr) for s,fr in atoms if s=="I"]
    write_scf(f"{d}/scf_empty.in",f"{name}_empty",cell,fw,["C","H","O",metal])
    write_scf(f"{d}/scf_i2only.in",f"{name}_i2only",cell,io,["I"])
    write_pp(f"{d}/pp_loaded_17.in",f"{name}_loaded",17,"rho17.cube")
    write_pp(f"{d}/pp_loaded_8.in",f"{name}_loaded",8,"elf.cube")
    write_pp(f"{d}/pp_loaded_0.in",f"{name}_loaded",0,"val.cube")
    write_pp(f"{d}/pp_empty_0.in",f"{name}_empty",0,"val.cube")
    write_pp(f"{d}/pp_i2only_0.in",f"{name}_i2only",0,"val.cube")
    # critic2 inputs
    open(f"{d}/bcp.cri","w").write("crystal rho17.cube\nload rho17.cube id rho\nauto\n")
    open(f"{d}/bcp_val.cri","w").write("crystal val.cube\nload val.cube id rho\nauto\n")
    open(f"{d}/ct.cri","w").write("crystal rho17.cube\nload rho17.cube id rho\nload deltarho.cube id dr\nintegrable dr\nbader\n")
    print(f"{name}: {len(atoms)} atoms ({len(fw)} fw + {len(io)} I), metal={metal}")
print("done")
