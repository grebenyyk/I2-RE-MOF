"""QTAIM 2D contour renderer for the I2@MOF BCP publication figures.

Self-contained: slices the cube (trilinear) for the plane electron density,
parses critic2 bcp.out for critical points (authoritative), traces bond paths
by 2D gradient ascent, and renders a publication-style figure with atoms,
bond paths, BCPs and ring CPs overlaid on log-spaced density contours.

Plane = orthonormal frame through 3 atoms, so plot axes are true Angstrom.
"""
import os, re, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
sys.path.insert(0, os.path.dirname(__file__))
from cubeio import read_cube_grid, interp_frac, frac_from_cart, cart_from_frac

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
ELCOL = {"I": "#7b3f9e", "Y": "#22e22c", "Tb": "#2cbfd8", "Gd": "#14c0c0",
         "O": "#d62728", "C": "#404040", "H": "#b0b0b0", "N": "#1f77b4"}
ELRAD = {"I": 0.50, "Y": 0.36, "Tb": 0.38, "Gd": 0.38, "O": 0.24,
         "C": 0.24, "H": 0.16, "N": 0.24}  # Å, draw scale


# -------------------------- bcp.out parsing ---------------------------------
def parse_bcp(path):
    """Parse critic2 bcp.out -> list of CP dicts: {n,type,frac,rho,lap}."""
    txt = open(path, errors="ignore").read()
    parts = re.split(r"\+ Critical point no\. (\d+)", txt)
    cps = []
    for i in range(1, len(parts), 2):
        b = parts[i + 1]
        typ = re.search(r"Type :\s*([^\n]+)", b)
        crys = re.search(r"Crystallographic coordinates:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", b)
        f = re.search(r"Field value \(f\):\s*(-?\d+\.\d+E[+-]\d+)", b)
        lap = re.search(r"Laplacian \(del2 f\):\s*(-?\d+\.\d+E[+-]\d+)", b)
        if not (typ and crys and f and lap):
            continue
        t = typ.group(1).strip()
        cps.append({
            "n": int(parts[i]), "type": t,
            "frac": np.array([float(crys.group(1)), float(crys.group(2)), float(crys.group(3))]),
            "rho": float(f.group(1)), "lap": float(lap.group(1)),
        })
    return cps


def cp_signature(t):
    """Map critic2 type string -> (rank,sig) tag."""
    t = t.strip().lower()
    m = re.match(r"\(\s*3\s*,\s*(-?\d+)\s*\)", t)
    if m:
        return int(m.group(1))
    return {"nucleus": -3, "bond": -1, "ring": 1, "cage": 3}.get(t.split()[0], None)


# -------------------------- plane geometry ----------------------------------
class Plane:
    """Orthonormal plot frame through fractional points r0,r1,r2 (3 atoms)."""
    def __init__(self, cell, f0, f1, f2, upad=1.6, vpad=1.6, N=420):
        self.cell = cell
        c0, c1, c2 = cart_from_frac(cell, f0), cart_from_frac(cell, f1), cart_from_frac(cell, f2)
        e1 = c1 - c0
        self.uhat = e1 / np.linalg.norm(e1)
        e2 = c2 - c0
        e2p = e2 - np.dot(e2, self.uhat) * self.uhat
        self.vhat = e2p / np.linalg.norm(e2p)
        self.nhat = np.cross(self.uhat, self.vhat)
        self.c0 = c0
        # determine plot window from the 3 corners (parallelogram) + pad
        corners = [c0, c1, c2, c1 + (c2 - c0)]
        us = [np.dot((c - c0), self.uhat) for c in corners]
        vs = [np.dot((c - c0), self.vhat) for c in corners]
        self.umin, self.umax = min(us) - upad, max(us) + upad
        self.vmin, self.vmax = min(vs) - vpad, max(vs) + vpad
        self.Nu = N
        self.Nv = max(int(round(N * (self.vmax - self.vmin) / (self.umax - self.umin))), 120)

    def project(self, cart):
        """cart (Å) or fractional-handled outside -> (u,v,h) Å in plane frame."""
        d = cart - self.c0
        return (np.dot(d, self.uhat), np.dot(d, self.vhat), np.dot(d, self.nhat))

    def field(self, grid):
        """Sample cube field on the plane -> 2D array (Nv,Nu), value in e/bohr3."""
        us = np.linspace(self.umin, self.umax, self.Nu)
        vs = np.linspace(self.vmin, self.vmax, self.Nv)
        U, V = np.meshgrid(us, vs)
        # plane points in cart Å
        pts = self.c0[None, None, :] + U[..., None] * self.uhat + V[..., None] * self.vhat
        flat = pts.reshape(-1, 3)
        # -> fractional (periodic), MIC wrap into [0,1)
        inv = np.linalg.inv(self.cell.T)
        frac = (flat @ inv.T) % 1.0
        vals = np.array([interp_frac(grid, f) for f in frac])
        return vals.reshape(self.Nv, self.Nu), us, vs

    def grad_uv(self, F, us, vs):
        """Gradient of F in (u,v) Å space. Returns (dFdu, dFdv)."""
        du = us[1] - us[0]; dv = vs[1] - vs[0]
        dFdu, dFdv = np.gradient(F, du, dv)
        return dFdu, dFdv


# -------------------------- bond path tracing -------------------------------
def trace_paths(plane, F, us, vs, bcp_uv, nuc_uv, nuc_cart, hmax=0.45):
    """From each BCP (u,v), gradient-ascend on the plane field toward the 2
    nearest projected nuclei. Returns list of polylines [(u,v),...]."""
    dFdu, dFdv = plane.grad_uv(F, us, vs)
    du = us[1] - us[0]; dv = vs[1] - vs[0]

    def grad_at(u, v):
        iu = (u - us[0]) / du; iv = (v - vs[0]) / dv
        if not (0 <= iu < plane.Nu - 1 and 0 <= iv < plane.Nv - 1):
            return None
        i0, j0 = int(iu), int(iv)
        fu, fv = iu - i0, iv - j0
        gu = (dFdu[j0, i0] * (1 - fu) * (1 - fv) + dFdu[j0, i0 + 1] * fu * (1 - fv) +
              dFdu[j0 + 1, i0] * (1 - fu) * fv + dFdu[j0 + 1, i0 + 1] * fu * fv)
        gv = (dFdv[j0, i0] * (1 - fu) * (1 - fv) + dFdv[j0, i0 + 1] * fu * (1 - fv) +
              dFdv[j0 + 1, i0] * (1 - fu) * fv + dFdv[j0 + 1, i0 + 1] * fu * fv)
        return np.array([gu, gv])

    paths = []
    for (bu, bv) in bcp_uv:
        # 2 nearest nuclei in plane
        dn = sorted(((bu - nu) ** 2 + (bv - nv) ** 2, k)
                    for k, (nu, nv) in enumerate(nuc_uv))[:2]
        for _, k in dn:
            nu, nv = nuc_uv[k]
            dirv = np.array([nu - bu, nv - bv]); dirv /= (np.linalg.norm(dirv) + 1e-12)
            r = np.array([bu, bv]) + dirv * 0.04
            line = [r.copy()]
            for _ in range(400):
                g = grad_at(r[0], r[1])
                if g is None:
                    break
                gn = np.linalg.norm(g)
                if gn < 1e-9:
                    break
                r = r + (g / gn) * 0.035
                line.append(r.copy())
                if (r[0] - nu) ** 2 + (r[1] - nv) ** 2 < 0.09:  # within 0.3 A of nucleus
                    break
            if len(line) > 3:
                paths.append(np.array(line))
    return paths


# -------------------------- panel render ------------------------------------
def render_panel(name, struct, plane_atoms, field="val.cube", hmax=0.55,
                 contour_field_levels=None, ax=None, title=None, show_lap=False):
    sdir = f"{BASE}/{struct}"
    cell, atoms, nvox, grid = read_cube_grid(f"{sdir}/{field}")
    cps = parse_bcp(f"{sdir}/bcp.out")
    f0 = frac_from_cart(cell, atoms[plane_atoms[0] - 1][1])
    f1 = frac_from_cart(cell, atoms[plane_atoms[1] - 1][1])
    f2 = frac_from_cart(cell, atoms[plane_atoms[2] - 1][1])
    plane = Plane(cell, f0, f1, f2)
    F, us, vs = plane.field(grid)
    F = np.maximum(F, 1e-9)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 6.2))
    # --- density: filled log contours + line contours ---
    if contour_field_levels is None:
        lev = np.array([2, 4, 8])[:, None] * 10.0 ** np.arange(-4, 0)  # 2,4,8 x 1e-4..1e-1
        lev = np.unique(np.sort(lev.ravel()))
        lev = lev[lev < F.max()]
    else:
        lev = contour_field_levels
    cmap = plt.get_cmap("Greys")
    cf = ax.contourf(us, vs, F, levels=np.concatenate([[lev[0] * 0.3], lev, [lev[-1] * 4]]),
                     cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=lev[0] * 0.3, vmax=lev[-1] * 4),
                     extend="max")
    cs = ax.contour(us, vs, F, levels=lev, colors="k", linewidths=0.45, alpha=0.55)
    # --- Laplacian zero-line (covalent vs closed-shell boundary) ---
    if show_lap:
        # recompute laplacian via grad of grad (coarse); optional, off by default
        pass

    # --- project atoms (within hmax of plane) ---
    proj = []
    for s, r in atoms:
        u, v, h = plane.project(r)
        if (plane.umin - 0.3 <= u <= plane.umax + 0.3 and
                plane.vmin - 0.3 <= v <= plane.vmax + 0.3 and abs(h) <= hmax):
            proj.append((s, u, v, h))
    for s, u, v, h in proj:
        ax.add_patch(Circle((u, v), ELRAD.get(s, 0.25), facecolor=ELCOL.get(s, "#888"),
                            edgecolor="k", lw=0.6, zorder=5))
        if s in ("I", "Y", "Tb", "Gd", "Cd"):
            ax.text(u, v, s, ha="center", va="center", fontsize=8.5, fontweight="bold",
                    color="white", zorder=6,
                    path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=1.6, foreground=ELCOL.get(s, "k"))])

    nuc_uv = [(plane.project(r)[0], plane.project(r)[1]) for s, r in atoms]
    nuc_frac = [frac_from_cart(cell, r) for _, r in atoms]

    def mic_d(f1, f2):
        d = (f1 - f2); d -= np.round(d)
        return np.linalg.norm(cell.T @ d)

    # --- BCPs (3,-1) and ring CPs (3,+1), with artifact filter ---
    bcp_uv, bcp_info = [], []
    for c in cps:
        sig = cp_signature(c["type"])
        if sig is None:
            continue
        cart = cart_from_frac(cell, c["frac"])
        u, v, h = plane.project(cart)
        if not (plane.umin - 0.2 <= u <= plane.umax + 0.2 and plane.vmin - 0.2 <= v <= plane.vmax + 0.2):
            continue
        if abs(h) > hmax:
            continue
        if sig == -1:  # bond CP
            # drop near-core spurious BCPs (heavy-atom cores, esp. I, spawn these):
            # a real BCP sits ~between 2 atoms at comparable distance; artifacts are
            # one-sided (very close to ONE nucleus). Filter: too close to nearest,
            # or 2nd-nearest much farther than 1st (ratio > 2.5).
            dn = sorted((mic_d(c["frac"], nf), k) for k, nf in enumerate(nuc_frac))
            r1, r2 = dn[0][0], dn[1][0]
            if r1 < 0.6 or r2 > 2.0 * r1:
                continue
            lab = "".join(sorted([atoms[dn[0][1]][0], atoms[dn[1][1]][0]]))
            bcp_uv.append((u, v)); bcp_info.append((u, v, c["rho"], c["lap"], lab))
            ax.plot(u, v, marker="D", ms=4.4, mfc="#ffd700", mec="k", mew=0.6, zorder=8)
        elif sig == 1:  # ring CP
            ax.plot(u, v, marker="o", ms=3.4, mfc="#1f77b4", mec="k", mew=0.4, zorder=7)
        elif sig == 3:  # cage CP
            ax.plot(u, v, marker="s", ms=2.8, mfc="#2ca02c", mec="k", mew=0.3, zorder=7)

    # --- bond paths ---
    paths = trace_paths(plane, F, us, vs, bcp_uv, nuc_uv, [r for _, r in atoms], hmax)
    for p in paths:
        ax.plot(p[:, 0], p[:, 1], color="#b8860b", lw=1.1, zorder=6, alpha=0.9)

    ax.set_aspect("equal"); ax.set_xlim(plane.umin, plane.umax); ax.set_ylim(plane.vmin, plane.vmax)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
    if title:
        ax.set_title(title, fontsize=10)
    return {"plane": plane, "F": F, "bcp_info": bcp_info, "atoms_proj": proj}


# -------------------------- two-panel figure --------------------------------
def two_panel():
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    r1 = render_panel("Y_act_cov", "Y_act", [76, 77, 6], field="rho17.cube",
                      ax=axes[0], title="Activated Y-MOF — covalent I–I (2.69 Å, ∇²ρ < 0)")
    r2 = render_panel("Y_deep_chain", "Y_deep", [76, 75, 73], field="rho17.cube",
                      ax=axes[1], title="Y-MOF (deep) — uniform polyiodide chain + I₂ (2.94 Å, ∇²ρ > 0)")
    # annotate the I-I BCPs on each panel with their rho/lap (total dens, e/bohr3)
    def label_ii(ax, info, pick="covalent"):
        ii = [t for t in info if t[4] == "II"]
        if not ii:
            return
        if pick == "covalent":
            tgt = min(ii, key=lambda t: t[3])  # most negative lap = most covalent
        else:
            tgt = max(ii, key=lambda t: t[3])  # closed-shell polyiodide
        u, v, rho, lap, _ = tgt
        ax.annotate(f"I–I BCP\nρ = {rho:.4f}\n∇²ρ = {lap:+.3f}", (u, v), fontsize=8,
                    xytext=(10, 10), textcoords="offset points", color="k",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#b8860b", alpha=0.85))
    label_ii(axes[0], r1["bcp_info"], "covalent")
    label_ii(axes[1], r2["bcp_info"], "closedshell")
    fig.suptitle("QTAIM topology of I–I interactions in Y-MOF", fontsize=12, y=0.99)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "qtaim_2panel_Y.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
    print("wrote", out)
    print("  panel0 I-I BCPs:", [(round(t[2], 4), round(t[3], 4)) for t in r1["bcp_info"] if t[4] == "II"])
    print("  panel1 I-I BCPs:", [(round(t[2], 4), round(t[3], 4)) for t in r2["bcp_info"] if t[4] == "II"])


if __name__ == "__main__":
    two_panel()
