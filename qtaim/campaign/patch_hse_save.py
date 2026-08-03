#!/usr/bin/env python3
"""Salvage HSE PDOS without the (intractable) PAW-EXX post-convergence step.

The HSE SCF converges fine and writes the converged wavefunctions to
{outdir}/{prefix}.wfc1 via disk_io='high' -- these ARE the HSE KS states
(ACE is exact at convergence). The slow step is only the post-SCF EXX *energy*
evaluation + final wfc collect, which we skip.

We build a complete .save dir from a fast PBE scf (same cell/ecut/nbnd), then
patch its data-file-schema.xml: replace PBE eigenvalues/occupations/Fermi with
the HSE values parsed from the HSE scf.out, and swap wfc1.dat for the HSE .wfc1.
projwfc then projects the HSE wfc onto the HSE eigenvalue axis = clean HSE PDOS.

usage: patch_hse_save.py <variant_dir>   (e.g. .../Y_str)
"""
import re, sys, os
HA2EV = 27.211386245988
d = sys.argv[1].rstrip("/")
name = os.path.basename(d)
prefix = f"{name}_hse"
save = f"{d}/tmp_hse/{prefix}.save"
schema = f"{save}/data-file-schema.xml"
hse_wfc = f"{d}/hse/hse_converged.wfc1"

# --- parse HSE eigenvalues + occupations + fermi from scf.out ---
t = open(f"{d}/hse/scf.out", errors="ignore").read()
parts = t.split("bands (ev):")
assert len(parts) >= 2, "no 'bands (ev):' block in scf.out"
eb = parts[-1].split("occupation numbers")[0]
eigs_ev = [float(x) for x in re.findall(r"-?\d+\.\d+", eb)]
op = t.split("occupation numbers")[-1]
occs = []
for ln in op.split("\n"):
    n = re.findall(r"-?\d+\.\d+", ln)
    if n and len(n) <= 8 and not any(c in ln for c in ["Fermi","energy","occup","iter","EXX","ACE","smearing"]):
        occs += [float(x) for x in n]
    elif occs:
        break
assert len(eigs_ev) == 230 and len(occs) == 230, f"got {len(eigs_ev)} eigs, {len(occs)} occs"
fe_ev = float(re.search(r"Fermi energy is\s+(-?\d+\.\d+)", t).group(1))
eigs_ha = [e / HA2EV for e in eigs_ev]
fe_ha = fe_ev / HA2EV

def fmt(x):
    return f"{x:+.15E}"
def block(vals, per=3):
    return "\n".join(" " + " ".join(fmt(v) for v in vals[i:i+per]) for i in range(0, len(vals), per))

s = open(schema).read()
s2 = re.sub(r'(<eigenvalues size="230">\n)(.*?)(\n\s*</eigenvalues>)',
            lambda m: m.group(1) + block(eigs_ha) + m.group(3), s, flags=re.S)
assert s2 != s, "eigenvalues block not found"
s2 = re.sub(r'(<occupations size="230">\n)(.*?)(\n\s*</occupations>)',
            lambda m: m.group(1) + block(occs) + m.group(3), s2, flags=re.S)
s2 = re.sub(r'<fermi_energy>.*?</fermi_energy>',
            f'<fermi_energy>{fe_ha:.15E}</fermi_energy>', s2)
open(schema, "w").write(s2)

# --- swap wfc1.dat for the converged HSE wfc ---
assert os.path.exists(hse_wfc), f"missing {hse_wfc}"
os.replace(hse_wfc, f"{save}/wfc1.dat")

print(f"patched {name}: 230 HSE eigs (Ha) + occs + fermi={fe_ev:.4f} eV")
print(f"  near-EF eigs (eV): {[round(e,3) for e in eigs_ev[189:195]]}")
print(f"  near-EF occs:      {[round(o,3) for o in occs[189:195]]}")
print(f"  swapped wfc1.dat <- HSE converged wfc1 ({os.path.getsize(save+'/wfc1.dat')} bytes)")
