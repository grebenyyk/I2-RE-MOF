#!/usr/bin/env python3
"""Chain-on-activated seeded from the UNIFORM chain (not herringbone).
For 13_Y, 06_Tb (and optionally 08_Gd): read the uniform-chain strained final
(<sys>/uniform_chain/POSCAR_final), take its I atoms, seed the ACTIVATED framework,
relax I only. Overwrites the herringbone-seeded <sys>/activated_i2/loaded.in."""
import os, numpy as np, ase.io
BASE="/workspace/dimitry/calc/qe/i2_mof"
SYSTEMS=["13_Y","06_Tb"]   # 08_Gd already done
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF",
     "O":"O.pbe-kjpaw.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF","N":"N.pbe-n-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"Gd":157.25,"H":1.008,"O":15.999,"I":126.904,"N":14.007,"Tb":158.925,"Y":88.906}
for s in SYSTEMS:
    fpos=f"{BASE}/{s}/activate_empty/POSCAR_activated_empty"
    cpos=f"{BASE}/{s}/uniform_chain/POSCAR_final"
    if not os.path.exists(cpos): print(f"{s}: no uniform_chain/POSCAR_final, skip"); continue
    fa=ase.io.read(fpos); ca=ase.io.read(cpos)
    fsym=fa.get_chemical_symbols(); ffrac=fa.get_scaled_positions(wrap=False); cell=fa.cell.array
    csym=ca.get_chemical_symbols(); cfrac=ca.get_scaled_positions(wrap=False)
    Ifrac=[cfrac[i] for i in range(len(csym)) if csym[i]=="I"]
    allsym=list(fsym)+["I"]*len(Ifrac); allfrac=list(ffrac)+Ifrac
    uniq=[]; [uniq.append(x) for x in allsym if x not in uniq]
    outdir=f"{BASE}/{s}/activated_i2"; os.makedirs(outdir,exist_ok=True)
    L=[f"&CONTROL","  calculation = 'relax'",f"  prefix = '{s}_act_i2_chain_uniform'",
       f"  outdir = '{outdir}/tmp'","  pseudo_dir = '/workspace/qe-7.5/pseudo'","  tprnfor = .true.",
       "  etot_conv_thr = 1.0d-5","  forc_conv_thr = 1.0d-3","/",
       "&SYSTEM","  ibrav = 0",f"  nat = {len(allsym)}",f"  ntyp = {len(uniq)}","  ecutwfc = 50","  ecutrho = 400",
       "  occupations = 'smearing'","  smearing = 'mv'","  degauss = 0.01","  nspin = 1",
       "  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
       "&ELECTRONS","  conv_thr = 1.0d-5","  diagonalization = 'david'","  diago_david_ndim = 4",
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
    # verify motif: uniform chain (≥5 short contacts in 2.7-3.2 range)
    P=np.array(Ifrac)@cell; inv=np.linalg.inv(cell)
    def mic(a,b): f=(a-b)@inv; f-=np.round(f); return np.linalg.norm(f@cell)
    short=sorted(round(mic(P[i],P[j]),2) for a,i in enumerate(range(len(P))) for j in range(a+1,len(P)) if mic(P[i],P[j])<3.4)
    print(f"{s}: uniform-chain seed -> {outdir}/loaded.in; {len(Ifrac)} I, short contacts {short}")
