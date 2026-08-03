import re
C="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"; T=C+"/kmesh_test/dfpt"
import os; os.makedirs(T,exist_ok=True)
raw=open(C+"/scf_loaded.in").read()
s=raw.replace("prefix = 'Y_deep_loaded'","prefix = 'Y_deep_dfpt4'")
s=s.replace("outdir = './tmp_loaded'","outdir = './tmp_4x4x2'")
s=re.sub(r"degauss\s*=\s*0\.\d+","degauss = 0.03",s)
s=re.sub(r"K_POINTS.*?\n.*?\n","K_POINTS automatic\n  4 4 2 0 0 0\n",s,count=1,flags=re.S)
open(T+"/scf_4x4x2.in","w").write(s)
print("wrote scf_4x4x2.in")
# verify
for line in open(T+"/scf_4x4x2.in"):
    if any(k in line for k in ["prefix","K_POINTS","degauss","ecutwfc","conv_thr"]):
        print("  ",line.rstrip())
