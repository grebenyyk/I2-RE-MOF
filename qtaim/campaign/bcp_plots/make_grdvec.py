"""Generate critic2 GRDVEC inputs for QTAIM plane contour plots.

Given a structure dir + 3 atom indices defining a plane, emit a .cri that
runs critic2 grdvec (2D density contours + bond paths + CPs in the plane)
and writes a JSON with the plane definition for the matplotlib renderer.

Plane convention (critic2): plane r0 r1 r2 -> parallelogram corners
r0 (origin), r1 (x-end), r2 (y-end); the plot spans r0 + u*(r1-r0) + v*(r2-r0).
"""
import os, json, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cubeio import read_cube_header, frac_from_cart, BOHR

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
BOHR = 0.5291772109

# ---- panel definitions: name -> (struct, [at# r0, at# r1, at# r2], scale, title) ----
PANELS = {
    "Y_act_cov": {
        "struct": "Y_act", "atoms": [76, 77, 6], "scale": (1.8, 1.9),
        "title": "Y-MOF activated: covalent I-I (I4-I5, 2.69 A)",
        "field": "val.cube",
    },
    "Y_deep_chain": {
        "struct": "Y_deep", "atoms": [76, 75, 73], "scale": (1.35, 1.45),
        "title": "Y-MOF deep: uniform polyiodide chain + I2 (I4-I2-I3 ... I1)",
        "field": "val.cube",
    },
    "Y_str_poly": {
        "struct": "Y_str", "atoms": [76, 75, 73], "scale": (1.30, 1.40),
        "title": "Y-MOF solvent-templated: branched polyiodide (I4-I2-I3 + I1)",
        "field": "val.cube",
    },
}


def emit(name, spec, outdir):
    struct = spec["struct"]; a0, a1, a2 = spec["atoms"]
    field = spec["field"]
    cell, atoms, _nv = read_cube_header(f"{BASE}/{struct}/{field}")
    f0 = frac_from_cart(cell, atoms[a0 - 1][1])
    f1 = frac_from_cart(cell, atoms[a1 - 1][1])
    f2 = frac_from_cart(cell, atoms[a2 - 1][1])
    sx, sy = spec["scale"]

    sdir = f"{BASE}/{struct}"
    outroot = os.path.join(outdir, name)  # absolute; outputs -> bcp_plots/
    cri = []
    cri.append(f"crystal {sdir}/{field}")
    cri.append(f"load {sdir}/{field} id rho")
    cri.append("auto")
    cri.append("grdvec")
    cri.append(f" root {outroot}")
    cri.append(f" plane {f0[0]:.8f} {f0[1]:.8f} {f0[2]:.8f}  "
               f"{f1[0]:.8f} {f1[1]:.8f} {f1[2]:.8f}  "
               f"{f2[0]:.8f} {f2[1]:.8f} {f2[2]:.8f}")
    cri.append(f" scale {sx} {sy}")
    cri.append(f" outcp 1.15 1.15")
    cri.append(f" hmax 0.25")          # bohr; CP/atom inclusion thickness
    cri.append(" bcpall")              # trace all bond paths (BCP->nuclei)
    # total density, 401x401 grid, 18 log contours, bond-resolving range (e/bohr3)
    cri.append(" contour rho 401 401 log 18 1.0e-4 5.0e-1")
    cri.append("end")
    crip = os.path.join(outdir, f"{name}.cri")
    open(crip, "w").write("\n".join(cri) + "\n")

    # plane JSON for matplotlib: fractional plane pts + cell + atom list
    js = {
        "struct": struct, "field": field, "title": spec["title"],
        "plane_atoms": [a0, a1, a2],
        "plane_atoms_sym": [atoms[a0 - 1][0], atoms[a1 - 1][0], atoms[a2 - 1][0]],
        "frac": {"r0": f0.tolist(), "r1": f1.tolist(), "r2": f2.tolist()},
        "scale": [sx, sy],
        "cell_A": cell.tolist(),
        "atoms": [[s, r.tolist()] for s, r in atoms],
    }
    jsonp = os.path.join(outdir, f"{name}_plane.json")
    open(jsonp, "w").write(json.dumps(js, indent=1))
    print(f"[{name}] {struct} plane atoms {spec['plane_atoms' ] if False else spec['atoms']}"
          f" -> {crip}")
    print(f"       r0={atoms[a0-1][0]} frac=({f0[0]:.4f},{f0[1]:.4f},{f0[2]:.4f})")
    print(f"       r1={atoms[a1-1][0]} frac=({f1[0]:.4f},{f1[1]:.4f},{f1[2]:.4f})")
    print(f"       r2={atoms[a2-1][0]} frac=({f2[0]:.4f},{f2[1]:.4f},{f2[2]:.4f})  scale={sx},{sy}")
    return crip


if __name__ == "__main__":
    outdir = os.path.dirname(__file__)
    names = sys.argv[1:] or ["Y_act_cov", "Y_deep_chain"]
    for n in names:
        emit(n, PANELS[n], outdir)
