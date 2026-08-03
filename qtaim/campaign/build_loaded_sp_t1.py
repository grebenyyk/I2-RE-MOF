#!/usr/bin/env python3
"""Single-point SCF on the vc-relaxed loaded geometry at Table-1 level
(degauss 0.01, conv_thr 1e-5 = loaded.in). Gives loaded_vc energy at the
exact Table-1 smearing, on the vc-relaxed geometry found at degauss 0.03.
"""
import re
MOF="/workspace/dimitry/calc/qe/i2_mof"; CAMP=f"{MOF}/qtaim/campaign"
SYS={"Gd":(f"{MOF}/08_Gd/loaded.in", f"{CAMP}/Gd_act/vc-relax_t1.out", f"{CAMP}/Gd_act"),
     "Tb":(f"{MOF}/06_Tb/loaded.in", f"{CAMP}/Tb_act/vc-relax_t1.out", f"{CAMP}/Tb_act"),
     "Y" :(f"{MOF}/13_Y/loaded.in", f"{CAMP}/Y_actc/vc-relax_t1.out", f"{CAMP}/Y_actc")}
def species_block(inp):
    m=re.search(r"ATOMIC_SPECIES\n(.*?)\n(?:CELL_PARAMETERS|ATOMIC_POSITIONS)", inp, re.S)
    return m.group(1).strip().split("\n")
def last_cell(out):
    t=open(out,errors="ignore").read()
    c=re.findall(r"CELL_PARAMETERS[^\n]*\n((?:\s*-?[\d.eE+-]+[^\n]*\n){3})",t)
    return [x.strip() for x in c[-1].strip().split("\n")] if c else None
def last_pos(out):
    t=open(out,errors="ignore").read(); b=re.split(r"ATOMIC_POSITIONS[^\n]*\n",t); a=[]
    for ln in b[-1].split("\n"):
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",ln)
        if m and len(a)<80: a.append((m.group(1),m.group(2),m.group(3),m.group(4)))
    return a
for m,(inp,out,tgt) in SYS.items():
    if "JOB DONE" not in open(out,errors="ignore").read():
        print(f"{m}: vc-relax_t1 not done yet"); continue
    spec=species_block(open(inp).read()); cell=last_cell(out); pos=last_pos(out)
    if not cell or len(pos)!=80: print(f"{m}: geometry incomplete ({len(pos)})"); continue
    L=[f"&CONTROL","  calculation = 'scf'",f"  prefix = '{m}_loadsp01'",
       f"  outdir = '{tgt}/tmp_sp01'","  pseudo_dir = '/workspace/qe-7.5/pseudo'",
       "  tstress = .false.","  tprnfor = .false.","/",
       "&SYSTEM","  ibrav = 0","  nat = 80","  ntyp = 5",
       "  ecutwfc = 50","  ecutrho = 400","  occupations = 'smearing'","  smearing = 'mv'",
       "  degauss = 0.01","  nspin = 1","  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
       "&ELECTRONS","  conv_thr = 1.0d-5","  diagonalization = 'david'","  mixing_mode = 'local-TF'",
       "  mixing_beta = 0.2","  mixing_ndim = 16","  electron_maxstep = 300","/",
       "ATOMIC_SPECIES"]+spec+["CELL_PARAMETERS angstrom"]+[f"  {c}" for c in cell]+\
      ["ATOMIC_POSITIONS crystal"]+[f"{s:<3} {x} {y} {z}" for (s,x,y,z) in pos]+\
      ["K_POINTS automatic","  1 1 1 0 0 0",""]
    open(f"{tgt}/scf_loaded_sp01.in","w").write("\n".join(L))
    print(f"{m}: loaded_sp01 single-point built (80 atoms, 50/400/dg0.01/conv1e-5)")
