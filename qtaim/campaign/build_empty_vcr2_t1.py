#!/usr/bin/env python3
"""Build empty-framework RELAX inputs at the vc-relaxed cell, Table-1 level
(50/400, degauss 0.01). Starts from the vc-relaxed loaded framework (iodine
stripped) so the empty framework finds its own equilibrium at the loaded cell.
PPs/species copied verbatim from campaign loaded.in.
"""
import re
MOF="/workspace/dimitry/calc/qe/i2_mof"
CAMP=f"{MOF}/qtaim/campaign"
SYS={"Gd":(f"{MOF}/08_Gd/loaded.in", f"{CAMP}/Gd_act/vc-relax_t1.out", f"{CAMP}/Gd_act"),
     "Tb":(f"{MOF}/06_Tb/loaded.in", f"{CAMP}/Tb_act/vc-relax_t1.out", f"{CAMP}/Tb_act"),
     "Y" :(f"{MOF}/13_Y/loaded.in", f"{CAMP}/Y_actc/vc-relax_t1.out", f"{CAMP}/Y_actc")}

def species_block(inp):
    m=re.search(r"ATOMIC_SPECIES\n(.*?)\n(?:CELL_PARAMETERS|ATOMIC_POSITIONS)", inp, re.S)
    return m.group(1).strip().split("\n") if m else None

def last_cell(out):
    txt=open(out,errors="ignore").read()
    cells=re.findall(r"CELL_PARAMETERS[^\n]*\n((?:\s*-?[\d.eE+-]+[^\n]*\n){3})", txt)
    return [c.strip() for c in cells[-1].strip().split("\n")] if cells else None

def last_pos(out, dropI=True):
    txt=open(out,errors="ignore").read()
    blocks=re.split(r"ATOMIC_POSITIONS[^\n]*\n", txt)
    atoms=[]
    for ln in blocks[-1].split("\n"):
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", ln)
        if m:
            s=m.group(1)
            if not(dropI and s=="I"):
                atoms.append((s,m.group(2),m.group(3),m.group(4)))
        if len(atoms)>=72 and dropI: break
    return atoms

for m,(inp,out,tgt) in SYS.items():
    spec=species_block(open(inp).read())
    cell=last_cell(out); pos=last_pos(out, dropI=True)
    if cell is None or not pos:
        print(f"{m}: vc-relax_t1.out not ready (no final cell/pos)"); continue
    uniq=[]
    for s,*_ in pos:
        if s not in uniq: uniq.append(s)
    L=[f"&CONTROL","  calculation = 'relax'",f"  prefix = '{m}_empty_vct2'",
       f"  outdir = '{tgt}/tmp_evct2'","  pseudo_dir = '/workspace/qe-7.5/pseudo'",
       "  tstress = .true.","  tprnfor = .true.","  nstep = 200",
       "  etot_conv_thr = 1.0d-5","  forc_conv_thr = 1.0d-3","/",
       "&SYSTEM","  ibrav = 0",f"  nat = {len(pos)}",f"  ntyp = {len(uniq)}",
       "  ecutwfc = 50","  ecutrho = 400","  occupations = 'smearing'","  smearing = 'mv'",
       "  degauss = 0.01","  nspin = 1","  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
       "&ELECTRONS","  conv_thr = 1.0d-6","  diagonalization = 'david'","  mixing_mode = 'local-TF'",
       "  mixing_beta = 0.2",
       "  mixing_ndim = 16",
       "  electron_maxstep = 300","/","&IONS","  ion_dynamics = 'bfgs'","/",
       "ATOMIC_SPECIES"]+spec+["CELL_PARAMETERS angstrom"]+[f"  {c}" for c in cell]+\
      ["ATOMIC_POSITIONS crystal"]+[f"{s:<3} {x} {y} {z}" for (s,x,y,z) in pos]+\
      ["K_POINTS automatic","  1 1 1 0 0 0",""]
    open(f"{tgt}/scf_empty_vct2.in","w").write("\n".join(L))
    print(f"{m}: empty_vct2 built ({len(pos)} fw atoms, 50/400/dg0.01, at vc cell)")
