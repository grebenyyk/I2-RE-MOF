import re,glob,os
A="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"
L="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/landscape_3variants"
DD="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/decomposition_Y3"
RYKJ=1312.75; N=4
def E(p):
    if not p.endswith(".out"): p+=".out"
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(p,errors="ignore").read(),re.M); return float(m[-1]) if m else None
def D3(p):
    if not p.endswith(".out"): p+=".out"
    m=re.findall(r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)",open(p,errors="ignore").read()); return float(m[-1]) if m else None
def r2(x): return round(x,2)
gas=E(A+"/rd_gas.out"); d3g=D3(A+"/rd_gas.out")

print("="*72); print("TABLE 1 - 9-system (per I2): ads=round(.out), disp=round(.out), elec=ads-disp"); print("="*72)
print("system             ads      disp      elec")
T1={}
for m in ["Y","Gd","Tb"]:
    for v in ["st","act","vc"]:
        Lp=E(f"{A}/rd_poly_{m}{v}"); Em=E(f"{A}/rd_empty_{m}{v}"); Dp=D3(f"{A}/rd_poly_{m}{v}"); De=D3(f"{A}/rd_empty_{m}{v}")
        ads=r2((Lp-Em-N*gas)/N*RYKJ); disp=r2((Dp-De-N*d3g)/N*RYKJ); elec=r2(ads-disp)
        T1[(m,v)]=(ads,disp,elec)
        print(f"I2@{m}-MOF{v:<6}{ads:>+9.2f}{disp:>+9.2f}{elec:>+9.2f}")

print("\n"+"="*72); print("TABLE 4 - Y Born-Haber: ads=Table1, asm/def=round(.out), anchor=ads-asm-def"); print("="*72)
print("System            asm    anchor      def      ads")
gasD=E(DD+"/rd_gas.out")
for v in ["st","act","vc"]:
    Lp=E(f"{DD}/rd_poly_Y{v}.out"); Em=E(f"{DD}/rd_empty_Y{v}.out"); I=E(f"{DD}/dc_isubl_Y{v}.out")
    Hs=E(f"{DD}/dc_hostld_Yvc.out") if v=="vc" else Em
    asm=r2((I-N*gasD)/N*RYKJ); df=r2((Hs-Em)/N*RYKJ); ads=T1[("Y",v)][0]; anc=r2(ads-asm-df)
    print(f"I2@Y-MOF{v:<6}{asm:>+8.2f}{anc:>+9.2f}{df:>+8.2f}{ads:>+9.2f}")

print("\n"+"="*72); print("TABLE 5 - landscape motif: total/disp=round(.out), elec=total-disp, diff=poly-mol"); print("="*72)
EMPT={v:E(f"{L}/rd_empty_Y{v}.out") for v in ["st","act","vc"]}; DEM={v:D3(f"{L}/rd_empty_Y{v}.out") for v in ["st","act","vc"]}
def cl(tag):
    if re.match(r'rd_poly_Y(act|vc|st)$',tag): return (re.match(r'rd_poly_Y(act|vc|st)',tag).group(1),"CIF")
    if re.match(r'hse_lpoly|hse_lmol',tag): return (re.search(r'Y(act|vc|st)',tag).group(1),"mol" if "lmol" in tag else "poly")
    if re.match(r'hmol_Y(act|vc|st)_S1$',tag): return (re.search(r'Y(act|vc|st)',tag).group(1),"poly")
    if re.match(r'rd_mol_Yact$',tag): return ("act","mol")
    if re.match(r'hmol_Yvc$',tag): return ("vc","mol")
rows=[]
for out in glob.glob(L+"/*.out"):
    tag=os.path.basename(out)[:-4]; c=cl(tag)
    if not c: continue
    ep,dp=E(out),D3(out)
    if ep is None: continue
    v,mt=c; ads=r2((ep-EMPT[v]-N*gas)/N*RYKJ); disp=r2((dp-DEM[v]-N*d3g)/N*RYKJ); elec=r2(ads-disp)
    rows.append((v,mt,ads,disp,elec))
print("Host       lowest poly (tot/disp/elec)        lowest mol                  diff        pref")
for v in ["st","act","vc"]:
    poly=[r for r in rows if r[0]==v and r[1] in ("poly","CIF")]; mol=[r for r in rows if r[0]==v and r[1]=="mol"]
    lp=min(poly,key=lambda r:r[2]); lm=min(mol,key=lambda r:r[2]); d=r2(lp[2]-lm[2])
    fp=lambda r:f"{r[2]:.2f}({r[3]:.2f}/{r[4]:+.2f})"
    print(f"Y-MOF{v:<6}{fp(lp):>30}{fp(lm):>26}{d:>+9.2f}   {'polyiodide' if d<=0 else 'molecular'}")

print("\n"+"="*72); print("TABLE S6 - breathing: st/act per I2=Table1; strain=round(.out); gain=4*(st2-act2); Est-Eact=strain+gain"); print("="*72)
print("System       ads_st    ads_act     strain      gain   Est-Eact")
for m in ["Gd","Tb","Y"]:
    st=T1[(m,"st")][0]; act=T1[(m,"act")][0]
    strain=r2((E(f"{A}/rd_empty_{m}st")-E(f"{A}/rd_empty_{m}act"))*RYKJ)
    gain=r2(4*(st-act)); est=r2(strain+gain)
    print(f"{m+'-MOF':<11}{st:>+9.2f}{act:>+9.2f}{strain:>+10.2f}{gain:>+10.2f}{est:>+10.2f}")
