import math
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test"
def bcps(path):
    out=[]
    for line in open(path,errors="ignore"):
        if "(3,-1)" not in line or "bond" not in line: continue
        f=line.split()
        try: x,y,z,rho=float(f[4]),float(f[5]),float(f[6]),float(f[9])
        except: continue
        out.append((x,y,z,rho))
    return out
def mind(p,q):
    s=0.0
    for i in range(3):
        d=abs(p[i]-q[i])%1.0
        if d>0.5: d=1.0-d
        s+=d*d
    return math.sqrt(s)
g=bcps(T+"/bcp_gamma_full.out"); a=bcps(T+"/bcp_k221_full.out"); b=bcps(T+"/bcp_k442_full.out")
def match(ref,other,tol=0.03):
    res=[]
    for r in ref:
        best=min(other,key=lambda o:mind(r,o)); dd=mind(r,best)
        if dd<tol: res.append((r[3],best[3],dd))
    return res
ma=match(g,a); mb=match(g,b)
print("Coordinate-matched BCPs (gamma->k221,k442), I-I range rho 0.018-0.072:")
print("%8s %9s %6s   %9s %6s %6s"%("rho_G","rho_k221","d%","rho_k442","d%","dxyz"))
s2=[];s4=[]
for i in range(min(len(ma),len(mb))):
    rg=ma[i][0]
    if not (0.018<=rg<=0.072): continue
    rg2,ra2,dda=ma[i]; _,rb2,ddb=mb[i]
    print("%8.4f %9.4f %+6.2f   %9.4f %+6.2f %6.3f"%(rg,ra2,(ra2-rg)/rg*100,rb2,(rb2-rg)/rg*100,ddb))
    s2.append(abs((ra2-rg)/rg)); s4.append(abs((rb2-rg)/rg))
print("\nI-I BCPs (rho 0.018-0.072): n=%d"%len(s2))
print("  k221 vs gamma: max=%.2f%% mean=%.2f%%"%(max(s2)*100,sum(s2)/len(s2)*100))
print("  k442 vs gamma: max=%.2f%% mean=%.2f%%"%(max(s4)*100,sum(s4)/len(s4)*100))
# also report ALL matched BCPs (full range) for context
s2a=[];s4a=[]
for i in range(min(len(ma),len(mb))):
    rg2,ra2,_=ma[i]; _,rb2,_=mb[i]
    s2a.append(abs((ra2-rg2)/rg2)); s4a.append(abs((rb2-rg2)/rg2))
print("\nALL matched bond CPs: n=%d"%len(s2a))
print("  k221 vs gamma: max=%.2f%% mean=%.2f%%"%(max(s2a)*100,sum(s2a)/len(s2a)*100))
print("  k442 vs gamma: max=%.2f%% mean=%.2f%%"%(max(s4a)*100,sum(s4a)/len(s4a)*100))
