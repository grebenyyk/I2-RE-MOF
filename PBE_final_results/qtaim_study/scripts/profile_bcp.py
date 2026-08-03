"""Recover I-I BCP rho/nabla2rho via the 1-D density profile along each bond.

For each I-I contact, scan rho(t) along the bond (t=0..1, MIC path) and take the
minimum -> the BCP density (robust; no Newton convergence needed). nabla2rho at
that point from metric-aware grid finite differences. This is the same
'bond-path rho(r) profile' used as the independent cross-check in the methods.
"""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from seed_bcp import load_cube, fields, Interp, BOHR


def profile_scan(cube, dcut=3.5):
    origin, H, n, Z, pos, rho = load_cube(cube)
    Ain = np.linalg.inv(np.array([n[i] * H[:, i] for i in range(3)]).T)
    A = np.array([n[i] * H[:, i] for i in range(3)]).T
    Hinv = np.linalg.inv(H)
    _, _, lap = fields(rho, H)
    ip = Interp(origin, Hinv, n, dict([('rho', rho), ('lap', lap)]))
    nuc = pos @ Ain
    Iidx = np.where(Z == 53)[0]
    Ilab = {i: k + 1 for k, i in enumerate(Iidx)}
    rows = []
    for a in range(len(Iidx)):
        for b in range(a + 1, len(Iidx)):
            ia, ib = Iidx[a], Iidx[b]
            disp = nuc[ib] - nuc[ia]; disp -= np.round(disp)   # ia -> ib (MIC)
            dist = np.linalg.norm(disp @ A) * BOHR
            if dist > dcut:
                continue
            ts = np.linspace(0.3, 0.7, 80)
            xs = np.array([(nuc[ia] + t * disp) @ A for t in ts])
            rhos = np.array([ip('rho', x) for x in xs])
            k = int(np.argmin(rhos))
            xbcp = xs[k]; rho_bcp = rhos[k]; lap_bcp = ip('lap', xbcp)
            rows.append((Ilab[ia], Ilab[ib], dist, rho_bcp, lap_bcp))
    rows.sort(key=lambda r: r[2])
    print(f"{cube.split('/')[-2]}: I-I BCP via profile minimum")
    print(f"{'Ia-Ib':>7} {'d(A)':>5} {'rho(au)':>8} {'rho(e/A3)':>9} {'lap(au)':>9} {'lap(e/A5)':>9}")
    for ia, ib, dist, rho_bcp, lap_bcp in rows:
        print(f"I{ia}-I{ib:<3} {dist:5.2f} {rho_bcp:8.4f} {rho_bcp/0.14818:9.3f} {lap_bcp:9.4f} {lap_bcp/0.041494:+9.3f}")


if __name__ == "__main__":
    for c in sys.argv[1:]:
        profile_scan(c)
