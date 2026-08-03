#!/usr/bin/env python3
"""Generate 18 I-only 4x4x2 relax inputs from the landscape_3variant Γ-relaxed seeds.
Recipe: framework frozen (if_pos 0 0 0), I free (1 1 1), damp, degauss=0.04, conv_thr=1e-6,
4x4x2, PBE-D3. Template namelists/ATOMIC_SPECIES from relax442_Yvc.in."""
import re,os
T_relax="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442/relax442_Yvc.in"
L="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/landscape_3variants"
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/landscape_442"
# (tag, seed_file) — 18 distinct basins, lowest-Γ-energy rep each
SEEDS=[
 ("vc_poly_S1",  f"{L}/hmol_Yvc_S1.in"),
 ("vc_poly_lS3", f"{L}/hse_lpoly_Yvc_S3.in"),
 ("vc_poly_lS6", f"{L}/hse_lpoly_Yvc_S6.in"),
 ("vc_mol_lS4",  f"{L}/hse_lmol_Yvc_S4.in"),
 ("vc_mol_ref",  f"{L}/hmol_Yvc.in"),
 ("vc_mol_lS5",  f"{L}/hse_lmol_Yvc_S5.in"),
 ("vc_mol_lS3",  f"{L}/hse_lmol_Yvc_S3.in"),
 ("act_mol_rd",  f"{L}/rd_mol_Yact.in"),
 ("act_poly_lS2",f"{L}/hse_lpoly_Yact_S2.in"),
 ("act_poly_lS5",f"{L}/hse_lpoly_Yact_S5.in"),
 ("act_mol_lS2", f"{L}/hse_lmol_Yact_S2.in"),
 ("act_mol_lS6", f"{L}/hse_lmol_Yact_S6.in"),
 ("st_poly_S1",  f"{L}/hmol_Yst_S1.in"),
 ("st_poly_lS5", f"{L}/hse_lpoly_Yst_S5.in"),
 ("st_poly_lS6", f"{L}/hse_lpoly_Yst_S6.in"),
 ("st_poly_lS4", f"{L}/hse_lpoly_Yst_S4.in"),
 ("st_mol_lS6",  f"{L}/hse_lmol_Yst_S6.in"),
 ("st_mol_lS3",  f"{L}/hse_lmol_Yst_S3.in"),
]
# template header (namelists+ATOMIC_SPECIES, before CELL_PARAMETERS) and footer (K_POINTS)
tmpl=open(T_relax).read()
HEADER=tmpl.split("CELL_PARAMETERS")[0]
FOOTER="K_POINTS"+tmpl.split("K_POINTS")[1]
if "nstep" not in HEADER:
    HEADER=HEADER.replace("  forc_conv_thr = 5.0d-3\n","  forc_conv_thr = 5.0d-3\n  nstep = 200\n",1)

# force recipe in header
HEADER=re.sub(r"degauss\s*=\s*[0-9.]+","degauss = 0.04",HEADER)
HEADER=re.sub(r"ion_dynamics\s*=\s*'[a-z]+'","ion_dynamics = 'damp'",HEADER)
HEADER=re.sub(r"nstep\s*=\s*[0-9]+","nstep = 200",HEADER)
def seed_atoms(fn):
    raw=open(fn).read()
    cell=re.search(r"CELL_PARAMETERS[^\n]*\n(.*?)\n\s*ATOMIC",raw,re.S).group(1).strip().split("\n")[:3]
    pos=re.search(r"ATOMIC_POSITIONS crystal\n(.*?)\nK_POINTS",raw,re.S).group(1).strip().split("\n")
    atoms=[(l.split()[0],l.split()[1:4]) for l in pos if l.split() and l.split()[0] in ("C","H","O","Y","I")]
    return cell,atoms
for tag,seed in SEEDS:
    cell,atoms=seed_atoms(seed)
    assert len(atoms)==80, f"{tag}: {len(atoms)} atoms (expect 80)"
    nI=sum(1 for sp,_ in atoms if sp=="I")
    assert nI==8, f"{tag}: {nI} I atoms (expect 8)"
    hdr=HEADER
    hdr=re.sub(r"prefix\s*=\s* '[^']*'",f"prefix = 'l442_{tag}'",hdr)
    hdr=re.sub(r"outdir\s*=\s*'[^']*'",f"outdir = './tmp_l442_{tag}'",hdr)
    celllines="\n".join("  "+"  ".join(f"{x:.10f}" for x in [float(v) for v in l.split()[:3]]) for l in cell)
    poslines=[]
    for sp,frac in atoms:
        flag="1 1 1" if sp=="I" else "0 0 0"
        poslines.append(f"{sp} {float(frac[0]):.12f} {float(frac[1]):.12f} {float(frac[2]):.12f} {flag}")
    out=hdr+"CELL_PARAMETERS angstrom\n"+celllines+"\nATOMIC_POSITIONS crystal\n"+"\n".join(poslines)+"\n"+FOOTER
    open(f"{OUT}/l442_{tag}.in","w").write(out)
print(f"wrote 18 l442_*.in to {OUT}")
# sanity: verify one
t="act_mol_rd"; s=open(f"{OUT}/l442_{t}.in").read()
dg=re.search(r"degauss = ([0-9.]+)",s).group(1)
io=re.search(r"ion_dynamics = '([a-z]+)'",s).group(1)
ns=(re.search(r"nstep = (\d+)",s) or [None,"(default)"])[1]
nl=len([l for l in s.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].split("\n") if l.strip() and l.split() and l.split()[0] in ("C","H","O","Y","I")])
print(f"sample {t}: degauss={dg} ion={io} nstep={ns} nat_lines={nl}")
