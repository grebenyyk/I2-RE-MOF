#!/usr/bin/env python3
"""Export the new molecular/clustered geometries as human-readable CIFs.
Sources: the HSE input .in files (= the relaxed geometries that were HSE'd),
or the relax .out (last ATOMIC_POSITIONS) for in-progress mmol_Yst."""
import re,os,numpy as np
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
DIR=f"{OUT}/new_geometries_cif"; os.makedirs(DIR,exist_ok=True)
def read_in(fn):
    L=open(fn).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and re.match(r'^[A-Z]',s[0]):
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    return cell,atoms
def read_out_last(fn):
    L=open(fn).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][-1]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    ai=[i for i,l in enumerate(L) if "ATOMIC_POSITIONS" in l][-1]
    atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES","End"): break
        if len(s)>=4 and re.match(r'^[A-Z]',s[0]):
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    return cell,atoms
def cell_geo(cell):
    a,b,c=[np.linalg.norm(cell[i]) for i in range(3)]
    al=np.degrees(np.arccos(np.dot(cell[1],cell[2])/(np.linalg.norm(cell[1])*np.linalg.norm(cell[2]))))
    be=np.degrees(np.arccos(np.dot(cell[0],cell[2])/(np.linalg.norm(cell[0])*np.linalg.norm(cell[2]))))
    ga=np.degrees(np.arccos(np.dot(cell[0],cell[1])/(np.linalg.norm(cell[0])*np.linalg.norm(cell[1]))))
    return a,b,c,al,be,ga
def write_cif(cell,atoms,name,note):
    a,b,c,al,be,ga=cell_geo(cell)
    s=[f"data_{name}",f"# {note}","# generated from the molecular/polyiodide HSE landscape scan (PBE/4x4x2 relaxed -> HSE/4x4x2)",
       "_symmetry_space_group_name_H-M   'P 1'","_symmetry_Int_Tables_number   1",
       f"_cell_length_a   {a:.6f}",f"_cell_length_b   {b:.6f}",f"_cell_length_c   {c:.6f}",
       f"_cell_angle_alpha   {al:.6f}",f"_cell_angle_beta   {be:.6f}",f"_cell_angle_gamma   {ga:.6f}",
       "loop_","_atom_site_type_symbol","_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z"]
    for e,x,y,z in atoms: s.append(f"{e} {x:.6f} {y:.6f} {z:.6f}")
    open(f"{DIR}/{name}.cif","w").write("\n".join(s)+"\n")
    return f"{name}.cif ({len(atoms)} atoms)"
# (source, output name, note). source = .in (HSE input = relaxed geom) or .out (in-progress relax)
JOBS=[
 (f"{OUT}/rd_mol_Yact.in","Y-MOFact_molecular-4I2","Y-MOF act: 4 molecular I2 (relax442_mol, RASPA-seeded). HSE: molecular favored by -16.7 kJ/I2 vs polyiodide."),
 (f"{OUT}/hmol_Yvc.in","Y-MOFvc_molecular-4I2","Y-MOF vc: 4 molecular I2 (void-placed). HSE: molecular +57.8 kJ/I2 above polyiodide (polyiodide favored)."),
 (f"{OUT}/hmol_Yvc_S1.in","Y-MOFvc_clustered-alt-polyiodide","Y-MOF vc: clustered alternative polyiodide (perturbed-polyiodide relax). HSE: -6.4 kJ/I2 below the original CIF polyiodide."),
 (f"{OUT}/hmol_Yst_S1.in","Y-MOFst_clustered-alt-polyiodide","Y-MOF st: clustered alternative polyiodide. HSE: -18.1 kJ/I2 below the original CIF polyiodide."),
]
for src,name,note in JOBS:
    src_read = read_out_last(src) if src.endswith(".out") else read_in(src)
    cell,atoms=src_read
    print("wrote",write_cif(cell,atoms,name,note))
# st molecular (mmol_Yst, in-progress relax) — from the current .out
if os.path.exists(f"{OUT}/mmol_Yst.out"):
    try:
        cell,atoms=read_out_last(f"{OUT}/mmol_Yst.out")
        if len(atoms)>=70: print("wrote",write_cif(cell,atoms,"Y-MOFst_molecular-4I2_INPROGRESS","Y-MOF st: 4 molecular I2 (void-placed), relax IN PROGRESS (energy not final)."))
    except Exception as e: print(f"mmol_Yst: {e}")
print(f"\nall CIFs in {DIR}/")
