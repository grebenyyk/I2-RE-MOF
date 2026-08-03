import re,os,numpy as np
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"; T=C+"/kmesh_test/raman_fd"
os.makedirs(T,exist_ok=True)
raw=open(C+"/scf_loaded.in").read()
cp=re.search(r"CELL_PARAMETERS.*?\n(.*?\n.*?\n.*?)\nATOMIC", raw, re.S).group(1).split("\n")
A=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)])  # rows a1,a2,a3
a1,a2,a3=A
Ainv=np.linalg.inv(A)
pos=re.search(r"ATOMIC_POSITIONS.*?crystal\n(.*?)\nK_POINTS", raw, re.S).group(1).strip().split("\n")
atoms=[(ln.split()[0], np.array([float(x) for x in ln.split()[1:4]])) for ln in pos]
I_idx=[i for i,(s,_) in enumerate(atoms) if s=="I"]
def cart(f): return f[0]*a1+f[1]*a2+f[2]*a3
r2=cart(atoms[I_idx[1]][1]); r3=cart(atoms[I_idx[2]][1])
d=r3-r2; bl=np.linalg.norm(d); u=d/bl
print(f"I2-I3 bond = {bl:.3f} A, axis=[{u[0]:.3f},{u[1]:.3f},{u[2]:.3f}]")
DELTA=0.01
def make_inp(kmesh,sign,tag):
    fshift=Ainv@(sign*DELTA*u); f2new=atoms[I_idx[1]][1]+fshift
    lines=[]
    for i,(s,f) in enumerate(atoms):
        ff=f2new if i==I_idx[1] else f
        lines.append(f"{s} {ff[0]:.12f} {ff[1]:.12f} {ff[2]:.12f}")
    inp=re.sub(r"prefix = 'Y_deep_loaded'",f"prefix = 'raman_fd_{tag}'",raw)
    inp=re.sub(r"outdir = './tmp_loaded'",f"outdir = './tmp_{tag}'",inp)
    inp=re.sub(r"ATOMIC_POSITIONS.*?crystal\n(.*?)\nK_POINTS", f"ATOMIC_POSITIONS crystal\n"+chr(10).join(lines)+f"\nK_POINTS", inp, flags=re.S)
    inp=re.sub(r"degauss\s*=\s*0\.\d+", "degauss = 0.03", inp)
    inp=re.sub(r"K_POINTS.*?\n.*?\n", f"K_POINTS automatic\n  {kmesh}\n", inp, count=1, flags=re.S)
    open(f"{T}/scf_{tag}.in","w").write(inp)
for km,klbl in [("1 1 1 0 0 0","g"),("4 4 2 0 0 0","k442")]:
    for sg,slbl in [(1,"p"),(-1,"m")]:
        make_inp(km,sg,f"{klbl}{slbl}")
print("made 4 inputs (gp/gm = Gamma +/-, k442p/k442m = 4x4x2 +/-)")
