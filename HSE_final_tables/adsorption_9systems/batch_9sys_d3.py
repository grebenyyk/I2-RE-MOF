#!/usr/bin/env python3
"""Batch D3 evaluation (pbe + hse06, version=4 BJ) on the 9-system table geometries
(adsorption_9systems: rd_poly/rd_empty/rd_gas). Saves {tag: {pbe,hse06}} JSON."""
import os, json, subprocess
A9 = "/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/adsorption_9systems"
RUN = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/run_d3.py"
ins = sorted(f for f in os.listdir(A9) if f.endswith(".in") and (f.startswith("rd_poly_") or f.startswith("rd_empty_") or f == "rd_gas.in"))
res = {}
for fn in ins:
    tag = fn[:-3]
    row = {}
    for func in ("pbe", "hse06"):
        r = subprocess.run(["python3", RUN, os.path.join(A9, fn), func, "4"],
                           capture_output=True, text=True)
        out = (r.stdout + r.stderr).strip().splitlines()
        out = out[-1] if out else ""
        try: row[func] = float(out)
        except ValueError: row[func] = "FAIL:" + out[-200:]
    res[tag] = row
json.dump(res, open(f"{A9}/d3_pbe_vs_hse06.json", "w"), indent=1)
print(f"{'tag':20s} {'D3_pbe':>12s} {'D3_hse06':>12s}")
for t, r in sorted(res.items()):
    p = r["pbe"]; h = r["hse06"]
    print(f"{t:20s} {p if isinstance(p,float) else 0:12.6f} {h if isinstance(h,float) else 0:12.6f}")
