#!/usr/bin/env python3
"""Stage-2 MOLECULAR: place DISCRETE I2 (RASPA/GCMC molecular seed from
poscars_unitcell/POSCAR_XX_..._I2loaded) on the ACTIVATED framework, relax I only.
This is the molecular-I2 analogue of 08_Gd/activated_i2 (-64.5 kJ/mol), to compare
against the chain-seeded result (build_activated_i2.py, -132 for Gd). Cell is the
same fixed cell, so RASPA fractional I2 coords transfer directly."""
import os, re, numpy as np, ase.io
BASE="/workspace/dimitry/calc/qe/i2_mof"
# system -> RASPA molecular I2-loaded POSCAR (the GCMC adsorption sites, discrete I2)
RASPA={"13_Y":"POSCAR_13_2020_Y_pts_3_ASR_1_I2loaded_1000Pa_UC",
       "06_Tb":"POSCAR_06_2019_Tb_pts_3_ASR_1_I2loaded_1000Pa_UC",
       "08_Gd":"POSCAR_08_2020_Gd_pts_3_ASR_1_I2loaded_1000Pa_UC",
       "12_Tb":"POSCAR_12_2020_Tb_nan_3_ASR_7_I2loaded_1000Pa_UC",
       "05_Cd":"POSCAR_05_oc8b00378_si_004_ASR_pacman_I2loaded_1000Pa_UC",
       "07_Cd":"POSCAR_07_oc8b00378_si_006_ASR_pacman_I2loaded_1000Pa_UC",
       "10_Cd":"POSCAR_10_oc8b00378_si_003_ASR_pacman_I2loaded_1000Pa_UC",
       "17_Cd":"POSCAR_17_2016_Cd_nan_3_ASR_1_I2loaded_1000Pa_UC"}
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF",
     "O":"O.pbe-kjpaw.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF","N":"N.pbe-n-kjpaw_psl.1.0.0.UPF","Cd":"Cd.pbe-n-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"Gd":157.25,"H":1.008,"O":15.999,"I":126.904,"N":14.007,"Cd":112.414,"Tb":158.925,"Y":88.906}
DEGAUSS={s:0.01 for s in RASPA};
for s in ["05_Cd","07_Cd","10_Cd","17_Cd"]: DEGAUSS[s]=0.005  # small-gap Cd
BIG={"12_Tb","05_Cd","07_Cd","10_Cd"}  # need david_ndim=2 to fit VRAM
def read_poscar_frac(p):
    """return (species_list, frac_coords, cell) from a VASP POSCAR (ase)."""
    a=ase.io.read(p); return a.get_chemical_symbols(), a.get_scaled_positions(wrap=False), a.cell.array
for s in RASPA:
    fpos=f"{BASE}/{s}/activate_empty/POSCAR_activated_empty"
    if not os.path.exists(fpos): print(f"{s}: no activated POSCAR, skip"); continue
    fsym,ffrac,cell=read_poscar_frac(fpos)   # framework (no I)
    rsym,rfrac,_=read_poscar_frac(f"{BASE}/poscars_unitcell/{RASPA[s]}")
    # RASPA I atoms = the molecular I2 seed
    Ifrac=[rfrac[i] for i in range(len(rsym)) if rsym[i]=="I"]
    nI=len(Ifrac)
    allsym=list(fsym)+["I"]*nI; allfrac=list(ffrac)+Ifrac
    uniq=[]; [uniq.append(x) for x in allsym if x not in uniq]
    outdir=f"{BASE}/{s}/activated_i2_mol"; os.makedirs(outdir,exist_ok=True)
    david=2 if s in BIG else 4
    L=[f"&CONTROL","  calculation = 'relax'",f"  prefix = '{s}_act_i2_mol'",
       f"  outdir = '{outdir}/tmp'","  pseudo_dir = '/workspace/qe-7.5/pseudo'","  tprnfor = .true.",
       "  etot_conv_thr = 1.0d-5","  forc_conv_thr = 1.0d-3","/",
       "&SYSTEM","  ibrav = 0",f"  nat = {len(allsym)}",f"  ntyp = {len(uniq)}","  ecutwfc = 50","  ecutrho = 400",
       "  occupations = 'smearing'","  smearing = 'mv'",f"  degauss = {DEGAUSS[s]}","  nspin = 1",
       "  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
       "&ELECTRONS","  conv_thr = 1.0d-5","  diagonalization = 'david'",f"  diago_david_ndim = {david}",
       "  mixing_mode = 'local-TF'","  mixing_beta = 0.2","  mixing_ndim = 8","  electron_maxstep = 300","/",
       "&IONS","  ion_dynamics = 'bfgs'","/",
       "ATOMIC_SPECIES"]
    for sp in uniq: L.append(f"{sp:<3s} {MASS[sp]:>10.3f}  {PSE[sp]}")
    L+=["CELL_PARAMETERS angstrom"]
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L+=["ATOMIC_POSITIONS crystal"]
    for sp,fr in zip(allsym,allfrac):
        ifp="0 0 0" if sp!="I" else "1 1 1"
        L.append(f"{sp:<3s} {fr[0]: .12f} {fr[1]: .12f} {fr[2]: .12f} {ifp}")
    L+=["K_POINTS gamma",""]
    open(f"{outdir}/loaded.in","w").write("\n".join(L))
    print(f"{s}: fw {len(fsym)} frozen + {nI} I movable (RASPA molecular seed), david={david} -> {outdir}/loaded.in")
