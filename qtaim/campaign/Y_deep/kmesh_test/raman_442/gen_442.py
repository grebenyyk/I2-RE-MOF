#!/usr/bin/env python3
"""Generate the 48 displaced 4x4x2 SCF inputs for the I-sublattice Hessian
(same procedure as the Gamma Raman: 8 I atoms x 3 axes x +/-0.01 A, framework
frozen, degauss 0.03, 60/800). Reads Y_deep/scf_loaded.in; writes to raman_442/.
Seeds each outdir with the Gamma charge density for consistent convergence."""
import re, os, numpy as np, shutil
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"
T=C+"/kmesh_test/raman_442"; os.makedirs(T, exist_ok=True)
DELTA=0.01
raw=open(C+"/scf_loaded.in").read()
# parse cell + fractional positions
cp=re.search(r"CELL_PARAMETERS.*?\n(.*?\n.*?\n.*?)\nATOMIC", raw, re.S).group(1).split("\n")
A=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)]); Ainv=np.linalg.inv(A)
pos_block=re.search(r"(ATOMIC_POSITIONS crystal\n.*?\n)K_POINTS", raw, re.S).group(1)
atoms=[(ln.split()[0], np.array([float(x) for x in ln.split()[1:4]])) for ln in pos_block.strip().split("\n")[1:]]
I_idx=[i for i,(s,_) in enumerate(atoms) if s=="I"]   # last 8
assert len(I_idx)==8
GSRC=C+"/tmp_loaded/Y_deep_loaded.save"   # Gamma save (charge density for seeding)
tags=[]
for ii,ia in enumerate(I_idx):       # I atom 0..7 (QE index 73..80)
    for ax in range(3):
        for sign,sl in [(1,"p"),(-1,"n")]:
            tag=f"I{ii+1}{chr(97+3+ax)}{sl}" if False else f"I{73+ii}{['x','y','z'][ax]}{sl}"
            dcart=np.zeros(3); dcart[ax]=sign*DELTA
            fshift=Ainv@dcart
            lines=[]
            for k,(s,f) in enumerate(atoms):
                ff=f+fshift if k==ia else f
                lines.append(f"{s} {ff[0]:.12f} {ff[1]:.12f} {ff[2]:.12f}")
            inp=raw.replace("prefix = 'Y_deep_loaded'", f"prefix = 'Y_deep_442_{tag}'")
            inp=re.sub(r"outdir = '[^']*'", f"outdir = './tmp_{tag}'", inp)
            inp=re.sub(r"ATOMIC_POSITIONS crystal\n.*?\nK_POINTS", "ATOMIC_POSITIONS crystal\n"+"\n".join(lines)+"\nK_POINTS", inp, flags=re.S)
            inp=re.sub(r"degauss\s*=\s*0\.\d+", "degauss = 0.03", inp)
            inp=re.sub(r"K_POINTS.*?\n.*?\n", "K_POINTS automatic\n  4 4 2 0 0 0\n", inp, count=1, flags=re.S)
            open(f"{T}/disp_{tag}.in","w").write(inp)
            # seed outdir with Gamma charge density (state-consistent starting guess)
            sd=f"{T}/tmp_{tag}/Y_deep_442_{tag}.save"; os.makedirs(sd, exist_ok=True)
            if os.path.exists(f"{GSRC}/charge-density.dat"):
                shutil.copy(f"{GSRC}/charge-density.dat", f"{sd}/charge-density.dat")
            tags.append((ii,ax,sign,tag))
with open(f"{T}/tags.txt","w") as f:
    for ii,ax,sign,tag in tags: f.write(f"{ii} {ax} {sign} {tag}\n")
print(f"generated {len(tags)} displaced 4x4x2 inputs + seeded Gamma density -> {T}")
