#!/usr/bin/env python3
"""Decode converged HSE wfc (.wfc1 direct-access) -> sequential wfc1.dat for
projwfc.x, reusing the PBE template's header+Miller (same cell/ecut/Gamma).
Assumes -np 1 local G-order == kdip order (verified sane on Y_str: HSE I-5p
tracks PBE with the expected gap-induced shift).
usage: decode_hse_pdos.py <variant_dir>   (e.g. .../Y_str)
"""
import struct, numpy as np, shutil, os, sys
D = sys.argv[1].rstrip("/") if len(sys.argv)>1 else "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_str"
B = f"{D}/hse"; name=os.path.basename(D); SAVE=f"{B}/hse_pdos_outdir/{name}_hse.save"
os.makedirs(SAVE, exist_ok=True)
shutil.copy(f"{B}/hse_keep/save/data-file-schema.xml", f"{SAVE}/data-file-schema.xml")
for fn in ["charge-density.dat","paw.txt"]:
    shutil.copy(f"{B}/hse_keep/save/{fn}", f"{SAVE}/{fn}")
for fn in os.listdir(f"{B}/pbe_template_save"):
    if fn.endswith(".UPF"): shutil.copy(f"{B}/pbe_template_save/{fn}", f"{SAVE}/{fn}")

def read_recs(path, n):
    f=open(path,"rb"); recs=[]
    for _ in range(n):
        m=struct.unpack('<i',f.read(4))[0]; d=f.read(m); f.read(4); recs.append(d)
    f.close(); return recs
hdr=read_recs(f"{B}/pbe_template_save/wfc1.dat",4)
ngw,igwx,npol,nbnd=struct.unpack('<4i',hdr[1])
hse=np.fromfile(f"{B}/hse_keep/hse_local.wfc1",dtype='<c16').astype(np.complex128)
assert len(hse)==igwx*nbnd, f"size mismatch {len(hse)} vs {igwx*nbnd}"
hse=hse.reshape(nbnd,igwx)
out=open(f"{SAVE}/wfc1.dat","wb")
def wrec(d): out.write(struct.pack('<i',len(d))); out.write(d); out.write(struct.pack('<i',len(d)))
for r in hdr: wrec(r)
for j in range(nbnd): wrec(np.ascontiguousarray(hse[j],dtype='<c16').tobytes())
out.close()
print(f"[{name}] wrote {SAVE}/wfc1.dat ({os.path.getsize(f'{SAVE}/wfc1.dat')} B)")
