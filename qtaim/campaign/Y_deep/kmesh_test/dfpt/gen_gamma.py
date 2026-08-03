import re
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"; T=C+"/kmesh_test/dfpt"
raw=open(C+"/scf_loaded.in").read()
# Gamma SCF for DFPT (degauss 0.03 to match FD-Gamma Raman)
s=raw.replace("prefix = 'Y_deep_loaded'","prefix = 'Y_deep_dfptg'")
s=s.replace("outdir = './tmp_loaded'","outdir = './tmp_gamma'")
s=re.sub(r"degauss\s*=\s*0\.\d+","degauss = 0.03",s)
s=re.sub(r"K_POINTS.*?\n.*?\n","K_POINTS automatic\n  1 1 1 0 0 0\n",s,count=1,flags=re.S)
open(T+"/scf_gamma.in","w").write(s)
# ph.x input for Gamma (q=Gamma), prefix dfptg
ph=f"""&inputph
  prefix = 'Y_deep_dfptg'
  outdir = './tmp_gamma'
  fildyn = 'dfptg.dyn'
  epsil = .false.
  amass(1) = 12.011
  amass(2) = 1.008
  amass(3) = 15.999
  amass(4) = 88.906
  amass(5) = 126.904
  tr2_ph = 1.0d-12
  verbosity = 'high'
/
0.0 0.0 0.0
"""
open(T+"/ph_g.in","w").write(ph)
print("wrote scf_gamma.in (Gamma, degauss 0.03) + ph_g.in")
