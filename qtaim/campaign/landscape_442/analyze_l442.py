#!/usr/bin/env python3
"""End-state analysis for the landscape_442 campaign. Run when configs converge.
- Ranks converged configs per variant by relaxed energy (lowest = true minimum).
- Writes CIFs (P1, fractional [0,1)) for ALL converged geometries -> cifs/.
- Generates sp02/sp_l442_<winner>.in (clean 0.02 SP) for the 3 variant winners.
- If those SP .out exist, computes ΔEads per I2 vs production Table 1."""
import re,math,numpy as np,os,glob
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/landscape_442"
L="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/landscape_3variants"
RY2KJ=13.6057*96.485
TAGS=["vc_poly_S1","vc_poly_lS3","vc_poly_lS6","vc_mol_lS4","vc_mol_ref","vc_mol_lS5","vc_mol_lS3",
      "act_mol_rd","act_poly_lS2","act_poly_lS5","act_mol_lS2","act_mol_lS6",
      "st_poly_S1","st_poly_lS5","st_poly_lS6","st_poly_lS4","st_mol_lS6","st_mol_lS3"]
def var(tag): return tag.split("_")[0]
def readf(p):
    try: return open(p,errors="ignore").read()
    except: return ""
def Eloaded(tag):
    m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(f"{T}/l442_{tag}.out"))
    return float(m[-1]) if m else None
def conv(tag):
    s=readf(f"{T}/l442_{tag}.out"); return ("bfgs converged" in s or "Damped Dynamics: convergence achieved" in s)
def cell_in(tag):
    raw=readf(f"{T}/l442_{tag}.in")
    cm=re.search(r"CELL_PARAMETERS[^\n]*\n(.*?)\n\s*ATOMIC",raw,re.S).group(1).strip().split("\n")[:3]
    return np.array([[float(x) for x in l.split()[:3]] for l in cm])
def atoms_out(tag):
    raw=readf(f"{T}/l442_{tag}.out")
    blk=raw.split("ATOMIC_POSITIONS")[-1].split("K_POINTS")[0].strip().split("\n")[1:]
    return [(l.split()[0],float(l.split()[1]),float(l.split()[2]),float(l.split()[3])) for l in blk if l.split() and l.split()[0] in ("C","H","O","Y","I")]
def write_cif(tag,cell,atoms):
    a1,a2,a3=cell; la,lb,lc=np.linalg.norm(a1),np.linalg.norm(a2),np.linalg.norm(a3)
    al=math.degrees(math.acos(np.dot(a2,a3)/(lb*lc))); be=math.degrees(math.acos(np.dot(a1,a3)/(la*lc))); ga=math.degrees(math.acos(np.dot(a1,a2)/(la*lb)))
    cnt={}; lines=[f"data_l442_{tag}",f"_cell_length_a {la:.6f}",f"_cell_length_b {lb:.6f}",f"_cell_length_c {lc:.6f}",
        f"_cell_angle_alpha {al:.4f}",f"_cell_angle_beta {be:.4f}",f"_cell_angle_gamma {ga:.4f}",
        "_symmetry_space_group_name_H-M 'P 1'","_symmetry_Int_Tables_number 1","loop_","_atom_site_label",
        "_atom_site_type_symbol","_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z"]
    for sp,x,y,z in atoms:
        cnt[sp]=cnt.get(sp,0)+1; x,y,z=x%1.0,y%1.0,z%1.0
        lines.append(f"{sp}{cnt[sp]} {sp} {x:.6f} {y:.6f} {z:.6f}")
    open(f"{T}/cifs/l442_{tag}_relaxed.cif","w").write("\n".join(lines)+"\n")
os.makedirs(f"{T}/cifs",exist_ok=True); os.makedirs(f"{T}/sp02",exist_ok=True)
# rank + CIFs
print("=== converged configs, ranked per variant by relaxed energy (degauss=0.04) ===")
winners={}
for v in ["vc","act","st"]:
    vt=[(t,Eloaded(t)) for t in TAGS if var(t)==v and conv(t)]
    vt=[t for t in vt if t[1] is not None]
    if not vt: print(f"\n{v}: none converged yet"); continue
    vt.sort(key=lambda x:x[1])
    print(f"\n{v}:")
    for t,e in vt:
        write_cif(t,cell_in(t),atoms_out(t))
        print(f"   {t:16} E={e:.6f} Ry")
    winners[v]=vt[0][0]
# generate 0.02 SP inputs for winners
for v,t in winners.items():
    raw=readf(f"{T}/l442_{t}.in")
    pos=atoms_out(t)
    hdr=raw.split("ATOMIC_POSITIONS")[0].replace("calculation = 'relax'","calculation = 'scf'")
    hdr=re.sub(r"degauss = [0-9.]+","degauss = 0.02",hdr)
    hdr=re.sub(r"&IONS.*?/","",hdr,flags=re.S)
    hdr=re.sub(r"prefix = '[^']*'",f"prefix = 'sp_l442_{t}'",hdr)
    hdr=re.sub(r"outdir = '[^']*'",f"outdir = './tmp_sp_l442_{t}'",hdr)
    poslines=[f"{sp} {x:.12f} {y:.12f} {z:.12f}" for sp,x,y,z in pos]
    open(f"{T}/sp02/sp_l442_{t}.in","w").write(hdr+"ATOMIC_POSITIONS crystal\n"+"\n".join(poslines)+"\nK_POINTS"+raw.split("K_POINTS")[1])
    print(f"winner {v}: {t} -> sp02/sp_l442_{t}.in (0.02 SP)")
# ΔEads if SPs done
def Esp(tag):
    m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(f"{T}/sp02/sp_l442_{tag}.out"))
    return float(m[-1]) if m else None
try:
    gas=float(re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(f"{L}/rd_gas.out"))[-1])
    empty={re.search(r'rd_empty_(\w+)',f).group(1):float(re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(f))[-1]) for f in glob.glob(f"{L}/rd_empty_*.out")}
except Exception: gas=None; empty={}
print("\n=== ΔEads (kJ/I2, clean 0.02 SP) ===")
for v,t in winners.items():
    el=Esp(t)
    if el is None or gas is None or v not in empty: print(f"  {v} ({t}): SP not done yet"); continue
    da=(el-empty[v]-4*gas)*RY2KJ/4.0
    print(f"  {v} ({t}): ΔEads = {da:.2f} kJ/I2")
