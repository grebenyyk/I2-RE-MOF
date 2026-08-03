#!/usr/bin/env python3
"""Batch HSE06-D3(BJ) evaluation on all landscape config geometries + references,
via the QE dftd3 driver (test.x, func=hse06 version=4). Saves {tag: Edisp_Ry} JSON."""
import re, glob, os, json, subprocess
D = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
def classify(tag):
    m = re.match(r'rd_poly_Y(act|vc|st)$', tag)
    if m: return m.group(1)
    if re.match(r'hse_lpoly|hse_lmol', tag): return 'x'
    if re.match(r'hmol_Y(act|vc|st)_S1$', tag): return 'x'
    if re.match(r'rd_mol_Yact$', tag): return 'x'
    if re.match(r'hmol_Yvc$', tag): return 'x'
    return None
tags = sorted(set(os.path.basename(o)[:-4] for o in glob.glob(f"{D}/*.out") if classify(os.path.basename(o)[:-4])))
tags += ["rd_empty_Yst", "rd_empty_Yact", "rd_empty_Yvc", "rd_gas"]
res = {}; miss = []
for tag in tags:
    inp = f"{D}/{tag}.in"
    if not os.path.exists(inp):
        res[tag] = None; miss.append(tag); continue
    r = subprocess.run(["python3", f"{D}/run_d3.py", inp, "hse06", "4"],
                       capture_output=True, text=True, cwd=D)
    out = (r.stdout + r.stderr).strip().splitlines()
    out = out[-1] if out else ""
    try:
        res[tag] = float(out)
    except ValueError:
        res[tag] = f"FAIL:{out[-300:]}"
json.dump(res, open(f"{D}/hse06d3_energies.json", "w"), indent=1)
ok = sum(1 for v in res.values() if isinstance(v, float))
print(f"HSE06-D3 computed: {ok}/{len(tags)} ok; missing .in: {miss}")
for t in ["rd_poly_Yst","rd_poly_Yact","rd_poly_Yvc","rd_empty_Yst","rd_empty_Yact","rd_empty_Yvc","rd_gas"]:
    print(f"  {t}: {res.get(t)}")
