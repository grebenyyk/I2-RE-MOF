#!/usr/bin/env python3
"""Build vc-relax inputs at Table-1 level (ecutwfc=50, ecutrho=400, degauss=0.01)
starting from the EXACT Table-1 chain geometry. Copies SYSTEM/ELECTRONS/species/cell
verbatim from the campaign loaded.in; only swaps in final positions + vc-relax flags.
"""
import re
MOF="/workspace/dimitry/calc/qe/i2_mof"
CAMP=f"{MOF}/qtaim/campaign"
# metal: (loaded.in, chain-loaded.out, target subdir, prefix)
SYS={"Gd":(f"{MOF}/08_Gd/loaded.in", f"{MOF}/08_Gd/activated_chain/loaded.out", f"{CAMP}/Gd_act","Gd_vct1"),
     "Tb":(f"{MOF}/06_Tb/loaded.in", f"{MOF}/06_Tb/activated_i2/loaded_ext.out", f"{CAMP}/Tb_act","Tb_vct1"),
     "Y" :(f"{MOF}/13_Y/loaded.in", f"{MOF}/13_Y/activated_i2_kickB/loaded.out", f"{CAMP}/Y_actc","Y_vct1")}

def section(text, name):
    m=re.search(rf"&{name}\b(.*?)(/)\s*\n", text, re.S)
    return m.group(1) if m else None

def extract_cell(inp):
    m=re.search(r"CELL_PARAMETERS[^\n]*\n((?:\s*-?[\d.\s]+\n){3})", inp, re.S)
    return m.group(1).strip().split("\n") if m else None

def extract_species(inp):
    m=re.search(r"ATOMIC_SPECIES\n(.*?)\n(?:CELL_PARAMETERS|ATOMIC_POSITIONS)", inp, re.S)
    return m.group(1).strip().split("\n") if m else None

def extract_last_pos(out):
    blocks=re.split(r"ATOMIC_POSITIONS[^\n]*\n", open(out,errors="ignore").read())
    last=blocks[-1]
    atoms=[]
    for ln in last.split("\n"):
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", ln)
        if m and len(atoms)<80:
            atoms.append((m.group(1),m.group(2),m.group(3),m.group(4)))
    return atoms

for m,(inp,out,tgt,pfx) in SYS.items():
    itext=open(inp).read()
    sys_s=section(itext,"SYSTEM"); el_s=section(itext,"ELECTRONS")
    cell=extract_cell(itext); spec=extract_species(itext); pos=extract_last_pos(out)
    nat=len(pos)
    metal=[l.split()[0] for l in spec if l.split()[0] in ("Gd","Tb","Y")][0]
    outdir="./tmp_vct1"
    L=[f"&CONTROL",
       f"  calculation = 'vc-relax'",
       f"  prefix = '{pfx}'",
       f"  outdir = '{outdir}'",
       f"  pseudo_dir = '/workspace/qe-7.5/pseudo'",
       f"  tstress = .true.",
       f"  tprnfor = .true.",
       f"  nstep = 200",
       f"  etot_conv_thr = 1.0d-5",
       f"  forc_conv_thr = 1.0d-3",
       "/",
       f"&SYSTEM", sys_s.strip(),
       "/",
       f"&ELECTRONS",
       f"  conv_thr = 1.0d-8",
       f"  diagonalization = 'david'",
       f"  mixing_beta = 0.2",
       f"  electron_maxstep = 300",
       "/",
       f"&IONS",
       f"  ion_dynamics = 'bfgs'",
       "/",
       f"&CELL",
       f"  cell_dynamics = 'bfgs'",
       f"  press_conv_thr = 0.5d0",
       f"  cell_dofree = 'all'",
       "/",
       "ATOMIC_SPECIES"]
    L+=spec
    L+=["CELL_PARAMETERS angstrom"]+[f"  {c.strip()}" for c in cell]
    L+=["ATOMIC_POSITIONS crystal"]
    L+=[f"{s:<3} {x} {y} {z}" for (s,x,y,z) in pos]
    L+=["K_POINTS automatic\n  1 1 1 0 0 0",""]
    open(f"{tgt}/vc-relax_t1.in","w").write("\n".join(L))
    print(f"{m}: vc-relax_t1.in written ({nat} atoms, metal={metal}, 50/400/dg0.01)")
