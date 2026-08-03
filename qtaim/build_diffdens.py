"""Build empty-framework and I-only-subsystem SCF inputs from 13Y_opt/scf.in,
same cell/60/800/Gamma so cube FFT grids match. Also the pp.x plot_num=21 inputs."""
import re
SRC="13Y_opt/scf.in"; BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim"
t=open(SRC).read()
namelist=t.split("ATOMIC_SPECIES")[0]
cell=t.split("CELL_PARAMETERS")[1].split("ATOMIC_POSITIONS")[0]
apos_block=t.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0]
kpts=t.split("K_POINTS")[1]
species=re.findall(r"(\w+)\s+[\d.]+\s+(\S+\.UPF)",t.split("CELL_PARAMETERS")[0].split("ATOMIC_SPECIES")[1])
# atom lines
atoms=[l for l in apos_block.splitlines() if re.match(r"\s*[A-Z][a-z]?\s+[-\d]",l)]
fw=[l for l in atoms if not l.strip().startswith("I ")]
ionly=[l for l in atoms if l.strip().startswith("I ")]
fw_spec=[(s,p) for s,p in species if s!="I"]
def write_in(path,prefix,sp,atms):
    o=open(path,"w")
    nl=namelist.replace("prefix = '13Y_opt'",f"prefix = '{prefix}'").replace("outdir = './tmp'",f"outdir = './tmp_{prefix}'")
    nl=re.sub(r"nat\s*=\s*\d+",f"nat = {len(atms)}",nl)
    nl=re.sub(r"ntyp\s*=\s*\d+",f"ntyp = {len(sp)}",nl)
    o.write(nl)
    o.write("ATOMIC_SPECIES\n")
    for s,p in sp: o.write(f"{s:<3s} { {'C':12.011,'H':1.008,'O':15.999,'Y':88.906,'I':126.904}[s]:8.3f}  {p}\n")
    o.write("CELL_PARAMETERS"+cell)
    o.write("ATOMIC_POSITIONS crystal\n")
    for l in atms: o.write(l.rstrip()+"\n")
    o.write("K_POINTS"+kpts)
    o.close()
write_in(f"{BASE}/13Y_empty/scf.in","13Y_empty",fw_spec,fw) if __import__("os").path.isdir(f"{BASE}/13Y_empty") else None
import os
os.makedirs(f"{BASE}/13Y_empty",exist_ok=True); os.makedirs(f"{BASE}/13Y_i2only",exist_ok=True)
write_in(f"{BASE}/13Y_empty/scf.in","13Y_empty",fw_spec,fw)
write_in(f"{BASE}/13Y_i2only/scf.in","13Y_i2only",[("I","I.pbe-n-kjpaw_psl.1.0.0.UPF")],ionly)
# pp.in for each (plot_num=21 all-electron, same as loaded)
for pre in ["13Y_empty","13Y_i2only"]:
    open(f"{BASE}/{pre}/pp.in","w").write(f"""&INPUTPP
  prefix = '{pre}'
  outdir = './tmp_{pre}'
  plot_num = 21
  filplot = 'ae_rho'
/
&PLOT
  nfile = 1
  filepp(1) = 'ae_rho'
  weight(1) = 1.0
  iflag = 3
  output_format = 6
  fileout = 'ae_rho.cube'
/
""")
print(f"empty: {len(fw)} fw atoms; species {[s for s,_ in fw_spec]}")
print(f"i2only: {len(ionly)} I atoms")
print("wrote 13Y_empty/{scf.in,pp.in} and 13Y_i2only/{scf.in,pp.in}")
