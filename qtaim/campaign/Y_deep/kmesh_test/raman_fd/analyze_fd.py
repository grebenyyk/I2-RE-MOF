import re,numpy as np
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"; T=C+"/kmesh_test/raman_fd"
DELTA=0.01
# recompute bond axis u (cartesian) from equilibrium
raw=open(C+"/scf_loaded.in").read()
cp=re.search(r"CELL_PARAMETERS.*?\n(.*?\n.*?\n.*?)\nATOMIC",raw,re.S).group(1).split("\n")
A=np.array([[float(x) for x in cp[i].split()[:3]] for i in range(3)]); a1,a2,a3=A
pos=re.search(r"ATOMIC_POSITIONS.*?crystal\n(.*?)\nK_POINTS",raw,re.S).group(1).strip().split("\n")
I=[np.array([float(x) for x in l.split()[1:4]]) for l in pos if l.startswith("I ")]
def cart(f): return f[0]*a1+f[1]*a2+f[2]*a3
r2,r3=cart(I[1]),cart(I[2]); d=r3-r2; u=d/np.linalg.norm(d)
I2_atom=74; I3_atom=75   # 1-indexed (last 8 atoms are I: 73..80)
def forces(path,atom):
    t=open(path,errors="ignore").read()
    m=re.search(rf"atom\s+{atom}\s+type\s+\d+\s+force\s*=\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",t)
    return np.array([float(m.group(i)) for i in (1,2,3)]) if m else None
res={}
for tag in ["gp","gm","k442p","k442m"]:
    f2=forces(f"{T}/scf_{tag}.out",I2_atom); f3=forces(f"{T}/scf_{tag}.out",I3_atom)
    res[tag]=(f2,f3)
# force constant (diagonal, I2 along bond axis): k = [F(-d)-F(+d)]/(2d)
def kdiag(pm_p,pm_m):  # force on I2 along u
    return (np.dot(res[pm_m][0],u)-np.dot(res[pm_p][0],u))/(2*DELTA)
kg=kdiag("gp","gm"); k4=kdiag("k442p","k442m")
# also off-diagonal: force on I3 along u
def koff(pm_p,pm_m): return (np.dot(res[pm_m][1],u)-np.dot(res[pm_p][1],u))/(2*DELTA)
offg=koff("gp","gm"); off4=koff("k442p","k442m")
# symmetric-stretch effective k ~ kdiag - koff (Ry/au per A); ratio is unit-free
print(f"I2-I3 bond axis u=[{u[0]:.3f},{u[1]:.3f},{u[2]:.3f}], d={DELTA} A")
print(f"\nForce on I2 (displaced) along bond, +d/-d:")
for tag in ["gp","gm","k442p","k442m"]:
    print(f"  {tag}: F_I2.u = {np.dot(res[tag][0],u):+.6f} Ry/au")
print(f"\nDiagonal force constant k(I2,I2)_u:  Gamma={kg:.4f}  4x4x2={k4:.4f} Ry/au/A   ratio k442/kG={k4/kg:.4f}")
print(f"Off-diagonal k(I3,I2)_u:             Gamma={offg:.4f}  4x4x2={off4:.4f}   ratio={off4/offg:.4f}")
estretch_g=kg-offg; estretch_4=k4-off4
print(f"Effective stretch constant (diag-off): Gamma={estretch_g:.4f} 4x4x2={estretch_4:.4f}  ratio={estretch_4/estretch_g:.4f}")
ratio=estretch_4/estretch_g
print(f"\n>>> Implied I2-I3 stretch frequency shift: dW/W = 0.5*(k442/kG - 1) = {0.5*(ratio-1)*100:+.2f}%")
print(f"    (e.g. a 150 cm-1 mode would shift by {0.5*(ratio-1)*150:+.2f} cm-1)")
