#!/usr/bin/env python3
"""Build SCF + pp inputs for the 6-system real-space campaign.
Unified settings: 60/800/Gamma, degauss 0.01, mv, mixing 0.2, david_ndim 4, conv 1e-8.
Per system: loaded SCF (full), empty SCF (framework), I-only SCF (8 I)."""
import os, re, sys
import numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"
OUT=f"{BASE}/qtaim/campaign"
SYS=[("Y_str","13_Y/uniform_chain/POSCAR_final",None),
     ("Tb_str","06_Tb/uniform_chain/POSCAR_final",None),
     ("Gd_str","08_Gd/uniform_chain/POSCAR_final",None),
     ("Y_act","13_Y/activated_i2/POSCAR_final_optimized",None),
     ("Tb_act",("06_Tb/activated_i2/loaded_ext.out","06_Tb/activated_i2/loaded.in"),"qe"),
     ("Gd_act","08_Gd/activated_chain/POSCAR_final_optimized",None)]
PPMAP={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
       "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
       "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Tb":158.925,"Gd":157.25,"I":126.904}
def read_poscar(p):
    import ase.io; a=ase.io.read(p)
    cell=a.cell.array; sym=a.get_chemical_symbols(); fr=a.get_scaled_positions(wrap=False)
    return cell,[(s,f) for s,f in zip(sym,fr)]
def read_qe_out(src):
    pout,pin = (src if isinstance(src,tuple) else (src,src))
    tin=open(pin,errors="ignore").read()
    cellp=re.findall(r"CELL_PARAMETERS[^\n]*\n[^\n]*\n[^\n]*\n[^\n]*\n",tin)
    lines=cellp[-1].strip().split("\n")[1:]
    cell=np.array([[float(x) for x in l.split()[:3]] for l in lines])
    t=open(pout,errors="ignore").read()
    blk=re.split(r"ATOMIC_POSITIONS",t)[-1]  # last block = final coords
    rows=[]
    for l in blk.splitlines():
        m=re.match(r"^\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",l)
        if m and m.group(1) in MASS: rows.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cell,rows
def write_scf(path,prefix,cell,atoms,species):
    fw=[s for s in atoms if s!="I"]; io=[a for a in atoms if a=="I"]
    ntyp=len({s for s,_ in atoms})
    outdir='./tmp_'+prefix.split('_')[-1]
    L=[f"""&CONTROL
  calculation = 'scf'
  prefix = '{prefix}'
  outdir = '{outdir}'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  tprnfor = .true.
  tstress = .true.
/
&SYSTEM
  ibrav = 0
  nat = {len(atoms)}
  ntyp = {ntyp}
  ecutwfc = 60
  ecutrho = 800
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.01
  nspin = 1
  vdw_corr = 'grimme-d3'
  dftd3_version = 4
/
&ELECTRONS
  conv_thr = 1.0d-8
  diagonalization = 'david'
  diago_david_ndim = 2
  mixing_mode = 'local-TF'
  mixing_beta = 0.2
  mixing_ndim = 16
  electron_maxstep = 300
/
ATOMIC_SPECIES
"""]
    for s in species: L.append(f"{s:<3s} {MASS[s]:8.3f}  {PPMAP[s]}\n")
    L.append("CELL_PARAMETERS angstrom\n")
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}\n")
    L.append("ATOMIC_POSITIONS crystal\n")
    for s,fr in atoms: L.append(f"{s:<3s} {fr[0]: .12f} {fr[1]: .12f} {fr[2]: .12f}\n")
    L.append("K_POINTS automatic\n  1 1 1 0 0 0\n")
    open(path,"w").write("".join(L))
def write_pp(path,prefix,plotnum,outfile):
    open(path,"w").write(f"""&INPUTPP
  prefix = '{prefix}'
  outdir = './tmp'
  plot_num = {plotnum}
  filplot = 'fld'
/
&PLOT
  nfile = 1
  filepp(1) = 'fld'
  weight(1) = 1.0
  iflag = 3
  output_format = 6
  fileout = '{outfile}'
/
""")
n=0
for name,src,kind in SYS:
    d=f"{OUT}/{name}"; os.makedirs(d,exist_ok=True)
    cell,atoms=(read_qe_out((f"{BASE}/{src[0]}",f"{BASE}/{src[1]}")) if kind=="qe" else read_poscar(f"{BASE}/{src}"))
    species=sorted({s for s,_ in atoms},key=lambda s:["C","H","O","Y","Tb","Gd","I"].index(s))
    write_scf(f"{d}/scf_loaded.in",f"{name}_loaded",cell,atoms,species)
    fw=[(s,fr) for s,fr in atoms if s!="I"]; ios=[(s,fr) for s,fr in atoms if s=="I"]
    write_scf(f"{d}/scf_empty.in",f"{name}_empty",cell,fw,[s for s in species if s!="I"])
    write_scf(f"{d}/scf_i2only.in",f"{name}_i2only",cell,ios,["I"])
    write_pp(f"{d}/pp_loaded_17.in",f"{name}_loaded",17,"rho17.cube")
    write_pp(f"{d}/pp_loaded_9.in",f"{name}_loaded",9,"gtau.cube")
    write_pp(f"{d}/pp_loaded_8.in",f"{name}_loaded",8,"elf.cube")
    write_pp(f"{d}/pp_loaded_0.in",f"{name}_loaded",0,"val.cube")
    write_pp(f"{d}/pp_empty_0.in",f"{name}_empty",0,"val.cube")
    write_pp(f"{d}/pp_i2only_0.in",f"{name}_i2only",0,"val.cube")
    nI=len(ios); nFw=len(fw)
    print(f"{name}: {len(atoms)} atoms ({nFw} fw + {nI} I), {len(species)} species {species}")
    n+=1
print(f"\nbuilt {n} systems in {OUT}")
