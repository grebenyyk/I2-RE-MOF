#!/usr/bin/env python3
"""Package the activated-pore campaign results into one self-contained dir:
  final_activated_results/
    table.txt                      (the decomposed table, activated + strained brackets)
    final/      <System>_<motif>.POSCAR      (16 converged geometries = the table energies)
    initial/    <System>_<motif>_seed.POSCAR (16 starting geometries / seeds)
    frameworks/ <System>_activated.POSCAR    (8 activated empty frameworks)
Names: Gd-MOF, Tb-MOF, Y-MOF, 12-Tb, 17-Cd, 07-Cd, 10-Cd, 05-Cd ; motifs: chain, mol."""
import re, os, shutil
import numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"; DST=f"{BASE}/final_activated_results"
for d in ["final","initial","frameworks"]: os.makedirs(f"{DST}/{d}",exist_ok=True)
SYS=[("08_Gd","Gd-MOF"),("06_Tb","Tb-MOF"),("13_Y","Y-MOF"),("12_Tb","12-Tb"),
     ("17_Cd","17-Cd"),("07_Cd","07-Cd"),("10_Cd","10-Cd"),("05_Cd","05-Cd")]
SPECIES=("C","Gd","H","O","I","N","Cd","Tb","Y")
PSE=["C","Gd","H","O","I","N","Cd","Tb","Y"]
def parse_qe_out_atoms(fn):
    """last ATOMIC_POSITIONS block -> [(sym,x,y,z)] in fractional crystal."""
    t=open(fn,errors="ignore").read(); best=None
    for b in re.split(r"ATOMIC_POSITIONS",t)[1:]:
        rows=[]
        for l in b.splitlines():
            m=re.match(r"^\s*([A-Z][a-z]?)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",l)
            if m and m.group(1) in SPECIES: rows.append((m.group(1),float(m.group(2)),float(m.group(3)),float(m.group(4))))
        if len(rows)>=50: best=rows
    return best
def parse_qe_in_atoms(fn):
    """ATOMIC_POSITIONS block from a QE input -> [(sym,x,y,z)] fractional."""
    t=open(fn,errors="ignore").read(); body=t.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0]
    rows=[]
    for l in body.splitlines():
        m=re.match(r"^\s*([A-Z][a-z]?)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",l)
        if m and m.group(1) in SPECIES: rows.append((m.group(1),float(m.group(2)),float(m.group(3)),float(m.group(4))))
    return rows
def cell_of(syscode):
    # read cell from the activated framework POSCAR
    import ase.io
    return ase.io.read(f"{BASE}/{syscode}/activate_empty/POSCAR_activated_empty").cell.array
def write_poscar(path, atoms, cell, title):
    sym=[a[0] for a in atoms]
    order=[];[order.append(s) for s in sym if s not in order]
    cnt=[sym.count(s) for s in order]
    L=[title,"1.0"]
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L.append("   "+"   ".join(f"{s:>4s}" for s in order)); L.append("   "+"   ".join(f"{c:4d}" for c in cnt)); L.append("Direct")
    for s in order:
        for (sp,x,y,z) in atoms:
            if sp==s: L.append(f"   {x: .12f}   {y: .12f}   {z: .12f}")
    open(path,"w").write("\n".join(L)+"\n")
# sources: (syscode, motif) -> (final_out_or_poscar, initial_in)
CHAIN_SRC={"08_Gd":("08_Gd/activated_chain/POSCAR_final_optimized","08_Gd/activated_chain/loaded.in"),
 "06_Tb":("06_Tb/activated_i2/loaded_ext.out","06_Tb/activated_i2/loaded.in"),
 "13_Y":("13_Y/activated_i2/loaded_ext.out","13_Y/activated_i2/loaded.in"),
 "12_Tb":("12_Tb/activated_i2/loaded_ext.out","12_Tb/activated_i2/loaded.in"),
 "17_Cd":("17_Cd/activated_i2/loaded_ext.out","17_Cd/activated_i2/loaded.in"),
 "07_Cd":("07_Cd/activated_i2/loaded.out","07_Cd/activated_i2/loaded.in"),
 "10_Cd":("10_Cd/activated_i2/loaded.out","10_Cd/activated_i2/loaded.in"),
 "05_Cd":("05_Cd/activated_i2/loaded_ext.out","05_Cd/activated_i2/loaded.in")}
MOL_SRC={"08_Gd":("08_Gd/activated_i2/POSCAR_final_optimized","08_Gd/activated_i2/loaded.in")}
for sc,_ in SYS:
    if sc!="08_Gd": MOL_SRC[sc]=(f"{sc}/activated_i2_mol/loaded.out",f"{sc}/activated_i2_mol/loaded.in")
def get_atoms(src):
    p=f"{BASE}/{src}"
    if "POSCAR" in src:  # any VASP POSCAR file -> ASE read (fractional)
        import ase.io; a=ase.io.read(p); sym=a.get_chemical_symbols(); fr=a.get_scaled_positions(wrap=False)
        return [(s,f[0],f[1],f[2]) for s,f in zip(sym,fr)]
    return parse_qe_out_atoms(p)
n_written=0
for sc,name in SYS:
    cell=cell_of(sc)
    # framework
    fw=get_atoms(f"{sc}/activate_empty/POSCAR_activated_empty")
    write_poscar(f"{DST}/frameworks/{name}_activated.POSCAR", fw, cell, f"{name} activated empty framework")
    for motif,SRC in [("chain",CHAIN_SRC),("mol",MOL_SRC)]:
        fin_src, ini_in = SRC[sc]
        fa=get_atoms(fin_src)
        write_poscar(f"{DST}/final/{name}_{motif}.POSCAR", fa, cell, f"{name} {motif} FINAL (converged, activated pore)")
        ia=parse_qe_in_atoms(f"{BASE}/{ini_in}")
        write_poscar(f"{DST}/initial/{name}_{motif}_seed.POSCAR", ia, cell, f"{name} {motif} INITIAL seed (activated fw + I2 seed)")
        n_written+=2
print(f"wrote {n_written} POSCARs (16 final + 16 initial) + 8 frameworks to {DST}")
