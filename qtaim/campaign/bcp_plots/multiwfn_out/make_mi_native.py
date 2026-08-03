"""M-I table from NATIVE critic2 output (seed pair dist 4.0) + verify I-I/I-O.
rho/lap from critic2 are a.u. (e/bohr3, e/bohr5) already. Distances in Angstrom.
"""
import sys, os, re, numpy as np
sys.path.insert(0, '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out')
from seed_bcp import load_cube, BOHR
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

ROWRE = re.compile(r"\(\s*3\s*,\s*-1\s*\)\s+bond\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)")
BOHR3, BOHR5 = BOHR**3, BOHR**5
TF = 0.3 * (3.0 * np.pi**2)**(2.0/3.0)
MZ = {'Y-MOFst':39,'Tb-MOFst':65,'Gd-MOFst':64,'Y-MOFact':39,'Tb-MOFact':65,'Gd-MOFact':64,'Y-MOFvc':39,'Tb-MOFvc':65,'Gd-MOFvc':64}
MSYM = {39:'Y',65:'Tb',64:'Gd'}
SYS = ['Y-MOFst','Tb-MOFst','Gd-MOFst','Y-MOFact','Tb-MOFact','Gd-MOFact','Y-MOFvc','Tb-MOFvc','Gd-MOFvc']
ARCH = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study'
SHIFTS = [(s1,s2,s3) for s1 in(-1,0,1) for s2 in(-1,0,1) for s3 in(-1,0,1)]

def gvh(rho,lap): g=TF*rho**(5/3)+lap/6; v=lap/4-2*g; return g,v,g+v

def analyze(sys):
    cube=f'{ARCH}/cubes/{sys}_rho.cube'; out=f'/tmp/bcpd4/{sys}.out'
    origin,H,n,Z,pos,rho=load_cube(cube); Lat=np.array([n[i]*H[:,i] for i in range(3)]); mz=MZ[sys]
    Iidx=[i for i in range(len(Z)) if Z[i]==53]; Ilab={i:k+1 for k,i in enumerate(Iidx)}
    def bcp_atom(c,a): return min(np.linalg.norm(c-pos[a]+s1*Lat[0]+s2*Lat[1]+s3*Lat[2]) for s1,s2,s3 in SHIFTS)*BOHR
    def atom_atom(a,b): return min(np.linalg.norm(pos[a]-pos[b]+s1*Lat[0]+s2*Lat[1]+s3*Lat[2]) for s1,s2,s3 in SHIFTS)*BOHR
    cnt={'I-I':0,'I-O':0,'M-I':0}; mi_reps={}
    for ln in open(out):
        m=ROWRE.search(ln)
        if not m: continue
        rho_b=float(m.group(4)); lap=float(m.group(5))
        if not (0.002<rho_b<0.08) or lap>0.1: continue
        f=np.array([float(m.group(1)),float(m.group(2)),float(m.group(3))]); cart=f@Lat
        ds=sorted([(bcp_atom(cart,a),a) for a in range(len(Z))])[:2]
        if ds[0][0]>2.6 or ds[1][0]/ds[0][0]>2: continue
        zp={Z[ds[0][1]],Z[ds[1][1]]}
        if zp=={53,53} and rho_b>=0.02: cnt['I-I']+=1
        elif zp=={53,8}: cnt['I-O']+=1
        elif zp=={mz,53}:
            cnt['M-I']+=1
            ii=[a for a in (ds[0][1],ds[1][1]) if Z[a]==53][0]
            ma=[a for a in (ds[0][1],ds[1][1]) if Z[a]==mz][0]
            key=min(tuple(np.round(f%1,3)),tuple(np.round((1-f)%1,3)))
            if key not in mi_reps:
                mi_reps[key]=(f'{MSYM[mz]}-I{Ilab[ii]}', atom_atom(ma,ii), rho_b, lap)
    return cnt, list(mi_reps.values())

recs=[]; print('system     I-I  I-O  M-I(loc)  M-I(dedup)')
for sys in SYS:
    cnt,mi=analyze(sys)
    for bond,d,rho,lap in mi:
        g,v,h=gvh(rho,lap)
        recs.append((sys,bond,d,rho,lap,g,v,h))
    print(f'{sys:9s} {cnt["I-I"]:4d} {cnt["I-O"]:4d} {cnt["M-I"]:6d}    {len(mi)}')

recs.sort(key=lambda r:(SYS.index(r[0]),r[2]))
import csv
with open(f'{ARCH}/tables/BCP_MI_native.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['System','Bond','d(A)','rho(e/A3)','lap(e/A5)','rho_BCP(e/bohr3)','G(au)','V(au)','H(au)'])
    for sys,bond,d,rho,lap,g,v,h in recs:
        w.writerow([sys,bond,f'{d:.2f}',f'{rho/BOHR3:.3f}',f'{lap/BOHR5:+.3f}',f'{rho:.4f}',f'{g:.4f}',f'{v:.4f}',f'{h:.4f}'])
doc=Document(); sec=doc.sections[0]; sec.orientation=WD_ORIENT.LANDSCAPE
sec.page_width,sec.page_height=sec.page_height,sec.page_width
for m in ('left_margin','right_margin','top_margin','bottom_margin'): setattr(sec,m,Inches(0.4))
doc.styles['Normal'].font.name='Calibri'; doc.styles['Normal'].font.size=Pt(9)
doc.add_paragraph().add_run('Table S. Y/Ln···I bond critical points (CRITIC2, seed pair dist 4.0 A). Ionic closed-shell contacts. G, V, H from Abramov + local virial (a.u.).').font.size=Pt(10)
hdr=['System','Bond','d (Å)','ρ (e/Å³)','∇²ρ (e/Å⁵)','ρ_BCP (e/bohr³)','G (a.u.)','V (a.u.)','H (a.u.)']
widths=[1.2,0.8,0.5,0.6,0.7,0.9,0.6,0.6,0.6]
def sc(c,t,s=8):
    c.text=''; p=c.paragraphs[0]; p.alignment=1; r=p.add_run(t); r.font.size=Pt(s); r.font.name='Calibri'
t=doc.add_table(rows=1,cols=len(hdr)); t.style='Table Grid'; t.autofit=False
for j,h in enumerate(hdr): sc(t.rows[0].cells[j],h)
for sys,bond,d,rho,lap,g,v,h in recs:
    vals=[sys,bond,f'{d:.2f}',f'{rho/BOHR3:.3f}',f'{lap/BOHR5:+.3f}',f'{rho:.4f}',f'{g:.4f}',f'{v:.4f}',f'{h:.4f}']
    cells=t.add_row().cells
    for j,val in enumerate(vals): sc(cells[j],val)
for j,wd in enumerate(widths):
    for row in t.rows: row.cells[j].width=Inches(wd)
doc.save(f'{ARCH}/tables/BCP_MI_native.docx')
print(f'\nwrote BCP_MI_native.docx/csv: {len(recs)} M-I BCPs (native CRITIC2)')
