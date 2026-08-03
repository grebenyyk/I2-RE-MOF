import numpy as np, math
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_actc/kmesh_test"
def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline()
        nat=int(f.readline().split()[0])
        gx=f.readline().split(); gy=f.readline().split(); gz=f.readline().split()
        nx,ny,nz=int(gx[0]),int(gy[0]),int(gz[0])
        for _ in range(nat): f.readline()
        data=np.fromstring(" ".join(f.readlines()),sep=" ").reshape(nx,ny,nz)
    return data,nx,ny,nz
# gamma I-I BCP coords from critic2 gamma output (rho in curated I-I range)
coords=[]
for line in open(T+"/bcp_gamma_full.out",errors="ignore"):
    if "(3,-1)" not in line or "bond" not in line: continue
    f=line.split()
    try: rho=float(f[9])
    except: continue
    if 0.045<=rho<=0.063:   # the 7 curated polyiodide I-I chain BCPs (0.0487-0.0616)
        coords.append((float(f[4]),float(f[5]),float(f[6]),rho))
print("gamma I-I BCP coords found:",len(coords))
def interp(data,nx,ny,nz,fc):
    # fc = fractional coord; trilinear with periodic wrap
    x=fc[0]%1.0*nx; y=fc[1]%1.0*ny; z=fc[2]%1.0*nz
    i=int(x); j=int(y); k=int(z); fx=x-i; fy=y-j; fz=z-k
    i0=i%nx; j0=j%ny; k0=k%nz; i1=(i+1)%nx; j1=(j+1)%ny; k1=(k+1)%nz
    c000=data[i0,j0,k0]; c100=data[i1,j0,k0]; c010=data[i0,j1,k0]; c110=data[i1,j1,k0]
    c001=data[i0,j0,k1]; c101=data[i1,j0,k1]; c011=data[i0,j1,k1]; c111=data[i1,j1,k1]
    return ((1-fx)*(1-fy)*(1-fz)*c000+(fx)*(1-fy)*(1-fz)*c100+(1-fx)*(fy)*(1-fz)*c010+(fx)*(fy)*(1-fz)*c110+
            (1-fx)*(1-fy)*(fz)*c001+(fx)*(1-fy)*(fz)*c101+(1-fx)*(fy)*(fz)*c011+(fx)*(fy)*(fz)*c111)
print("reading cubes...");
g,_,_,_=read_cube(T+"/rho17_gamma.cube")
a,_,_,_=read_cube(T+"/rho17_k221.cube")
b,nx,ny,nz=read_cube(T+"/rho17_k442.cube")
print("%9s %9s %9s %9s %7s %7s"%("rho_gamma","rho_k221","rho_k442","(critic2G)","k221%","k442%"))
r2=[];r4=[]
for x,y,z,rho_c in coords:
    rg=interp(g,nx,ny,nz,(x,y,z)); ra=interp(a,nx,ny,nz,(x,y,z)); rb=interp(b,nx,ny,nz,(x,y,z))
    print("%9.4f %9.4f %9.4f %9.4f %+7.2f %+7.2f"%(rg,ra,rb,rho_c,(ra-rg)/rg*100,(rb-rg)/rg*100))
    r2.append(abs((ra-rg)/rg)); r4.append(abs((rb-rg)/rg))
print("\nI-I BCP rho (trilinear interp at gamma BCP coords): n=%d"%len(r2))
print("  k221 vs gamma: max=%.2f%% mean=%.2f%%"%(max(r2)*100,sum(r2)/len(r2)*100))
print("  k442 vs gamma: max=%.2f%% mean=%.2f%%"%(max(r4)*100,sum(r4)/len(r4)*100))
