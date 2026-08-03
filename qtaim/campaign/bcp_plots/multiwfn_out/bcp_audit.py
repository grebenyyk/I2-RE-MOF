"""Audit located vs missing BCPs per dir: list geometric I-I contacts and mark
which have a located (3,-1) BCP. Proves whether short missing bonds are
location-failures (no CP found) or filter-casualties (CP found, later dropped).
"""
import sys, re, numpy as np, glob, os

BOHR = 0.5291772109
ROWRE = re.compile(
    r"\(\s*3\s*,\s*(-?\d)\s*\)\s+(nucleus|bond|ring|cage)\s+"
    r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+"
    r"[0-9.Ee+-]+\s+([0-9.Ee+-]+)")


def cube_cell(cube):
    L = open(cube).readlines()
    nat = int(L[2].split()[0])
    A = []
    for i in (3, 4, 5):
        p = L[i].split(); n = int(p[0]); A.append([n * float(p[1]), n * float(p[2]), n * float(p[3])])
    A = np.array(A)
    Z = []; pos = []
    for k in range(6, 6 + nat):
        p = L[k].split(); Z.append(int(float(p[0]))); pos.append([float(p[2]), float(p[3]), float(p[4])])
    return A, np.array(Z), np.array(pos)


def parse_bonds(bcpout):
    out = []; in_tbl = False
    for ln in open(bcpout, errors="ignore"):
        s = ln.strip()
        if s.startswith("# ncp") and "cryst" in s: in_tbl = True; continue
        if s.startswith("# ncp"): in_tbl = False
        if not in_tbl: continue
        m = ROWRE.search(s)
        if m and int(m.group(1)) == -1 and m.group(2) == "bond":
            f = np.array([float(m.group(3)), float(m.group(4)), float(m.group(5))])
            out.append((f, float(m.group(6)), float(m.group(7))))
    return out


def mic(f1, f2, A):
    d = f1 - f2; d -= np.round(d)
    return np.linalg.norm(d @ A) * BOHR


def two_nearest(frac, Ain, nucl_frac):
    d = frac - nucl_frac; d -= np.round(d)
    ds = np.linalg.norm(d @ Ain, axis=1) * BOHR
    o = np.argsort(ds)[:2]
    return o[0], o[1], ds[o[0]], ds[o[1]]


def audit(bcpout, cube, dcut=3.6):
    A, Z, pos = cube_cell(cube)
    Ain = np.linalg.inv(A)
    nucl = pos @ Ain
    Iidx = np.where(Z == 53)[0]
    Ilab = {i: k + 1 for k, i in enumerate(Iidx)}      # I1..I8 in cube order
    # geometric I-I contacts (MIC)
    contacts = []
    for a in range(len(Iidx)):
        for b in range(a + 1, len(Iidx)):
            d = mic(nucl[Iidx[a]], nucl[Iidx[b]], A)
            if d < dcut:
                contacts.append((Ilab[Iidx[a]], Ilab[Iidx[b]], d))
    contacts.sort(key=lambda t: t[2])
    # located I-I BCPs
    bcps = []
    for f, rho, lap in parse_bonds(bcpout):
        i, j, di, dj = two_nearest(f, Ain, nucl)
        if i in set(Iidx) and j in set(Iidx):
            bcps.append((min(Ilab[i], Ilab[j]), max(Ilab[i], Ilab[j]), rho, lap, di + dj))
    found_pairs = {(a, b) for a, b, *_ in bcps}
    return contacts, bcps, found_pairs


if __name__ == "__main__":
    dirs = sorted(glob.glob("Y_*")) if len(sys.argv) < 2 else sys.argv[1:]
    for d in dirs:
        bcp = os.path.join(d, "bcp.out"); cube = os.path.join(d, "val.cube")
        if not (os.path.exists(bcp) and os.path.exists(cube)):
            print(f"[{d}] missing bcp.out/val.cube"); continue
        contacts, bcps, found = audit(bcp, cube)
        print(f"\n===== {d} =====")
        print(f"I-I located BCPs: {len(bcps)} | geometric contacts <{3.6}Å: {len(contacts)}")
        print("  geometric I-I contacts (I_a-I_b  d/Å  BCP?) :")
        for a, b, dd in contacts:
            mark = f"FOUND rho={next((r for x,y,r,*_ in bcps if (x,y)==(a,b)),None):.4f}" if (a,b) in found else "*** MISSING ***"
            flag = "  <- short" if dd < 3.0 else ""
            print(f"    I{a}-I{b}  {dd:5.2f}  {mark}{flag}")
