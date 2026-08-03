"""Seed-refine I-I BCPs on a density cube (recover BCPs the auto-scan missed).

Loads rho17/val.cube, computes rho, gradient, Hessian and Laplacian on the
native (non-orthogonal) grid via metric-aware finite differences, then Newton-
refines grad(rho)=0 from the geometric midpoint of every short I-I contact.
Reports rho, lap, Hessian signature -> (3,-1) confirmation. VALIDATED first
against bonds critic2 already found.
"""
import sys, numpy as np

BOHR = 0.5291772109


def load_cube(path):
    L = open(path).readlines()
    nat = int(L[2].split()[0]); origin = np.array([float(x) for x in L[2].split()[1:4]])
    hv = []
    for i in (3, 4, 5):
        p = L[i].split(); hv.append([float(p[1]), float(p[2]), float(p[3])]); n = int(p[0])
    n1 = int(L[3].split()[0]); n2 = int(L[4].split()[0]); n3 = int(L[5].split()[0])
    H = np.array(hv).T                    # columns = voxel vectors h1,h2,h3
    Z = []; pos = []
    for k in range(6, 6 + nat):
        p = L[k].split(); Z.append(int(float(p[0]))); pos.append([float(p[2]), float(p[3]), float(p[4])])
    data = np.array("".join(L[6 + nat:]).split(), dtype=float).reshape(n1, n2, n3)  # k-fastest
    return origin, H, (n1, n2, n3), np.array(Z), np.array(pos), data


def fields(rho, H):
    """rho, gradient(cart), Hessian(cart,6), lap on the native grid (periodic)."""
    r = rho
    gu = np.stack([(np.roll(r, -1, a) - np.roll(r, 1, a)) / 2.0 for a in (0, 1, 2)])  # (3,n1,n2,n3)
    d2 = [np.roll(r, -1, a) - 2 * r + np.roll(r, 1, a) for a in (0, 1, 2)]
    def mix(a, b):
        return (np.roll(np.roll(r, -1, a), -1, b) - np.roll(np.roll(r, -1, a), 1, b)
                - np.roll(np.roll(r, 1, a), -1, b) + np.roll(np.roll(r, 1, a), 1, b)) / 4.0
    Hu = np.array([[d2[0], mix(0, 1), mix(0, 2)],
                   [mix(0, 1), d2[1], mix(1, 2)],
                   [mix(0, 2), mix(1, 2), d2[2]]])     # (3,3,n1,n2,n3)
    Hinv = np.linalg.inv(H); HinvT = Hinv.T
    N = r.size
    gu_v = gu.reshape(3, N)
    gx = (HinvT @ gu_v).reshape(3, *r.shape)            # cartesian gradient
    Hu_v = Hu.reshape(3, 3, N)
    Hx_v = np.einsum('ab,bcN,cd->adN', HinvT, Hu_v, Hinv)   # cartesian Hessian
    Hx = Hx_v.reshape(3, 3, *r.shape)
    lap = np.trace(Hx, axis1=0, axis2=1)                # (n1,n2,n3)
    return gx, Hx, lap


class Interp:
    def __init__(self, origin, Hinv, n, arrays):
        self.origin = origin; self.Hinv = Hinv; self.n = n; self.arrays = arrays
    def __call__(self, key, x):
        v = self.Hinv @ (x - self.origin)
        a = self.arrays[key]
        return self._tri(a, v)
    def _tri(self, a, v):
        n1, n2, n3 = self.n
        fi, fj, fk = v[0] % n1, v[1] % n2, v[2] % n3
        i0 = int(fi); j0 = int(fj); k0 = int(fk)
        dfi, dfj, dfk = fi - i0, fj - j0, fk - k0
        i1 = (i0 + 1) % n1; j1 = (j0 + 1) % n2; k1 = (k0 + 1) % n3
        c000 = a[i0, j0, k0]; c100 = a[i1, j0, k0]; c010 = a[i0, j1, k0]; c110 = a[i1, j1, k0]
        c001 = a[i0, j0, k1]; c101 = a[i1, j0, k1]; c011 = a[i0, j1, k1]; c111 = a[i1, j1, k1]
        c00 = c000 * (1 - dfi) + c100 * dfi; c10 = c010 * (1 - dfi) + c110 * dfi
        c01 = c001 * (1 - dfi) + c101 * dfi; c11 = c011 * (1 - dfi) + c111 * dfi
        c0 = c00 * (1 - dfj) + c10 * dfj; c1 = c01 * (1 - dfj) + c11 * dfj
        return c0 * (1 - dfk) + c1 * dfk


def newton(seed, ip, n=80, tol=1e-4):
    x = seed.copy()
    for _ in range(n):
        g = np.array([ip('g' + a, x) for a in 'xyz'])
        if np.linalg.norm(g) < tol:
            break
        Hx = np.array([[ip('h' + a + b, x) for b in 'xyz'] for a in 'xyz'])
        try:
            dx = np.linalg.solve(Hx, g)
        except np.linalg.LinAlgError:
            break
        step = np.linalg.norm(dx)
        if step > 0.5:
            dx *= 0.5 / step
        x = x - dx
    return x, np.linalg.norm(g)


def classify(ip, x):
    Hx = np.array([[ip('h' + a + b, x) for b in 'xyz'] for a in 'xyz'])
    ev = np.linalg.eigvalsh(Hx)
    neg = int((ev < 0).sum())
    sig_val = (3 - neg) - neg            # QTAIM signature = n_pos - n_neg
    sig = {-3: '(3,-3) NCP', -1: '(3,-1) BCP', 1: '(3,+1) ring', 3: '(3,+3) cage'}.get(sig_val, f'? ev={ev}')
    return sig, ev


def main(cube):
    origin, H, n, Z, pos, rho = load_cube(cube)
    Hinv = np.linalg.inv(H)
    Ain = np.linalg.inv(np.array([n[i] * H[:, i] for i in range(3)]).T)  # fractional lattice
    gx, Hx, lap = fields(rho, H)
    Hxs = {a + b: Hx[i, j] for i, a in enumerate('xyz') for j, b in enumerate('xyz')}
    ip = Interp(origin, Hinv, n, dict([('rho', rho), ('lap', lap)] +
                                      [('g' + a, gx[i]) for i, a in enumerate('xyz')] +
                                      [('h' + k, v) for k, v in Hxs.items()]))
    Iidx = np.where(Z == 53)[0]
    Ilab = {i: k + 1 for k, i in enumerate(Iidx)}
    nuc_f = pos @ Ain
    print(f"{cube}: {len(Iidx)} I atoms; metric-aware grid fields built")
    rows = []
    for a in range(len(Iidx)):
        for b in range(a + 1, len(Iidx)):
            ia, ib = Iidx[a], Iidx[b]
            d = nuc_f[ia] - nuc_f[ib]; d -= np.round(d)
            dist = np.linalg.norm(d @ np.array([n[i] * H[:, i] for i in range(3)]).T) * BOHR
            if dist > 3.5:
                continue
            mid = (nuc_f[ia] + d / 2.0) @ np.array([n[i] * H[:, i] for i in range(3)]).T  # bohr cart
            x, gnorm = newton(mid, ip)
            rho_bcp = ip('rho', x); lap_bcp = ip('lap', x)
            sig, ev = classify(ip, x)
            # 2 nearest nuclei
            df = x @ Ain - pos @ Ain; df -= np.round(df)
            nd = np.linalg.norm(df @ np.array([n[i] * H[:, i] for i in range(3)]).T, axis=1) * BOHR
            o = np.argsort(nd)[:2]
            rows.append((Ilab[ia], Ilab[ib], dist, sig, rho_bcp, lap_bcp * 1 / BOHR**5, gnorm,
                         Ilab.get(o[0], '?'), Ilab.get(o[1], '?'), nd[o[0]], nd[o[1]]))
    rows.sort(key=lambda r: r[2])
    print(f"{'Ia-Ib':>7} {'d':>5} {'sig':>10} {'rho(au)':>8} {'lap(e/A5)':>9} {'rho(e/A3)':>9} {'conv':>4}  near1 near2")
    for r in rows:
        conv = "ok" if r[6] < 1e-3 else "**"
        rhoA3 = r[4] / BOHR**3
        print(f"I{r[0]}-I{r[1]:<3} {r[2]:5.2f} {r[3]:>10} {r[4]:8.4f} {r[5]:+9.4f} {rhoA3:9.3f} {conv:>4}  I{r[7]}({r[9]:.2f}) I{r[8]}({r[10]:.2f})")


if __name__ == "__main__":
    main(sys.argv[1])
