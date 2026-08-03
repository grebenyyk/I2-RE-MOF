import numpy as np
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"
T=C+"/kmesh_test"
def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        gx=f.readline().split(); gy=f.readline().split(); gz=f.readline().split()
        nx,ny,nz=int(gx[0]),int(gy[0]),int(gz[0])
        for _ in range(nat): f.readline()
        return np.fromstring(" ".join(f.readlines()),sep=" ").reshape(nx,ny,nz),nx,ny,nz
g,nx,ny,nz=read_cube(C+"/elf.cube")          # Γ ELF
a,_,_,_=read_cube(T+"/elf_k221.cube")
b,_,_,_=read_cube(T+"/elf_k442.cube")
print(f"grids: gamma {g.shape} k221 {a.shape} k442 {b.shape} (must match)")
print("\n=== ELF: Γ vs k-mesh (RMS, direct subtraction, same grid) ===")
for tag,d in [("k221",a),("k442",b)]:
    diff=g-d
    rms=np.sqrt((diff**2).mean()); mx=np.abs(diff).max()
    print(f"  gamma vs {tag}: RMS_dELF={rms:.4e} ({rms*100:.3f}% abs, ELF in [0,1]) | max|dELF|={mx:.4f}")
print("\n=== WHERE does ELF change? (by ELF level) ===")
print(f"{'ELF region':30} {'RMS_dELF(k442)':>16}")
for lbl,mask in [("core/lone-pair (ELF>0.8)", g>0.8),
                 ("localized (0.6<ELF<0.8)",(g>0.6)&(g<=0.8)),
                 ("metallic channel (0.3<ELF<0.6)",(g>0.3)&(g<=0.6)),
                 ("low/delocalized (ELF<0.3)", g<=0.3)]:
    if mask.sum()>0:
        rms=np.sqrt(((g-b)[mask]**2).mean())
        print(f"{lbl:30} {rms:16.4e}")
# ELF at the I-I BCP coords (does the bond ELF shift with k?)
print("\n=== ELF at gamma I-I BCP coordinates ===")
import re
coords=[]
for line in open(T+"/bcp_gamma_full.out",errors="ignore"):
    if "(3,-1)" not in line or "bond" not in line: continue
    f=line.split()
    try: rho=float(f[9])
    except: continue
    if 0.045<=rho<=0.063: coords.append((float(f[4]),float(f[5]),float(f[6])))
def interp(e,fc):
    x=fc[0]%1.0*nx;y=fc[1]%1.0*ny;z=fc[2]%1.0*nz
    i=int(x);j=int(y);k=int(z);fx=x-i;fy=y-j;fz=z-k
    i0=i%nx;j0=j%ny;k0=k%nz;i1=(i+1)%nx;j1=(j+1)%ny;k1=(k+1)%nz
    def c(ii,jj,kk): return e[ii,jj,kk]
    return ((1-fx)*(1-fy)*(1-fz)*c(i0,j0,k0)+fx*(1-fy)*(1-fz)*c(i1,j0,k0)+(1-fx)*fy*(1-fz)*c(i0,j1,k0)+fx*fy*(1-fz)*c(i1,j1,k0)+
            (1-fx)*(1-fy)*fz*c(i0,j0,k1)+fx*(1-fy)*fz*c(i1,j0,k1)+(1-fx)*fy*fz*c(i0,j1,k1)+fx*fy*fz*c(i1,j1,k1))
print("%9s %9s %9s %7s %7s"%("ELF_G","ELF_k221","ELF_k442","k221%","k442%"))
r2=[];r4=[]
for fc in coords:
    eg=interp(g,fc);ea=interp(a,fc);eb=interp(b,fc)
    print("%9.4f %9.4f %9.4f %+7.2f %+7.2f"%(eg,ea,eb,(ea-eg),(eb-eg)))
    r2.append(abs(ea-eg));r4.append(abs(eb-eg))
print("\nELF at I-I BCPs: k221 mean|dELF|=%.4f max=%.4f | k442 mean|dELF|=%.4f max=%.4f"%(sum(r2)/len(r2),max(r2),sum(r4)/len(r4),max(r4)))
