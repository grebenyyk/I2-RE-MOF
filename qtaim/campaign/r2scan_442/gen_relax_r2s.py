#!/usr/bin/env python3
"""Build r2SCAN-NC I-only relax inputs for Yact/Yst from the converged SP inputs.
I-only relax: framework (C/H/O/Y) frozen via if_pos 0 0 0, 8 I free 1 1 1.
Matches production relax_442 philosophy (isolate functional's effect on I positions).
Usage: python3 gen_relax_r2s.py <ecutwfc> <ecutrho>"""
import sys, re
E  = int(sys.argv[1]) if len(sys.argv) > 1 else 70
RHO = int(sys.argv[2]) if len(sys.argv) > 2 else E*8
OUT = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/r2scan_442"
IONS_BLOCK = """&IONS
  ion_dynamics = 'damp'
/
"""
for sysname, sp in [("Yact","Yact_r2s_pdojo.in"), ("Yst","Yst_r2s_pdojo.in")]:
    s = open(f"{OUT}/{sp}").read()
    s = s.replace("calculation = 'scf'", "calculation = 'relax'")
    s = re.sub(r"prefix\s*=\s*'[^']*'", f"prefix = 'rel_{sysname}_r2s_e{E}'", s)
    s = re.sub(r"outdir\s*=\s*'[^']*'", f"outdir = './tmp_rel_{sysname}_r2s_e{E}'", s)
    s = re.sub(r"ecutwfc\s*=\s*\d+", f"ecutwfc = {E}", s)
    s = re.sub(r"ecutrho\s*=\s*\d+", f"ecutrho = {RHO}", s)
    s = re.sub(r"nstep\s*=\s*\d+", "nstep = 200", s)
    if "forc_conv_thr" not in s:
        s = s.replace("etot_conv_thr = 1.0d-5", "etot_conv_thr = 1.0d-5\n  forc_conv_thr = 5.0d-3")
    # inject &IONS (damp) before ATOMIC_SPECIES if absent
    if "&IONS" not in s:
        s = s.replace("ATOMIC_SPECIES", IONS_BLOCK + "ATOMIC_SPECIES", 1)
    # if_pos flags: framework frozen, I free
    m = re.search(r"(ATOMIC_POSITIONS crystal\n)(.*?)(\nK_POINTS)", s, re.S)
    nI = 0
    if m:
        body = m.group(2)
        lines = []
        for ln in body.strip().split("\n"):
            t = ln.split()
            if not t or t[0] not in ("C","H","O","Y","I"):
                lines.append(ln); continue
            if t[0] == "I": nI += 1; flag = "1 1 1"
            else: flag = "0 0 0"
            lines.append(f"{t[0]} {float(t[1]):.12f} {float(t[2]):.12f} {float(t[3]):.12f} {flag}")
        s = s.replace(m.group(0), m.group(1) + "\n".join(lines) + m.group(3))
    fn = f"{OUT}/rel_{sysname}_r2s_e{E}.in"
    open(fn, "w").write(s)
    print(f"wrote {fn}  (ecutwfc={E} ecutrho={RHO}, relax+damp, I-only free [{nI} I], framework frozen)")
