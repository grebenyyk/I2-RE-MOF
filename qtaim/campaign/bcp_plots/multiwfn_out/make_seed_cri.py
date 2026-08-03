"""Write a critic2 input seeding the MIC-midpoint of every short I-I contact.

critic2's `CP x y z` relaxes a Newton CP-search from fractional coords (x,y,z);
seeding at bond midpoints recovers shallow BCPs the `auto` scan misses.
"""
import sys, numpy as np

BOHR = 0.5291772109


def cube_parse(cube):
    L = open(cube).readlines()
    nat = int(L[2].split()[0])
    hv = [np.array([float(x) for x in L[i].split()[1:4]]) for i in (3, 4, 5)]
    n = [int(L[i].split()[0]) for i in (3, 4, 5)]
    A = np.array([n[i] * hv[i] for i in range(3)]).T            # cols = a1,a2,a3 (bohr)
    Ain = np.linalg.inv(A)
    Z = []; frac = []
    for k in range(6, 6 + nat):
        p = L[k].split(); Z.append(int(float(p[0])))
        frac.append(np.array([float(p[2]), float(p[3]), float(p[4])]) @ Ain)
    return A, np.array(Z), np.array(frac)


cube = sys.argv[1]; dcut = float(sys.argv[2]) if len(sys.argv) > 2 else 3.4
A, Z, frac = cube_parse(cube)
Iidx = np.where(Z == 53)[0]
Ilab = {i: k + 1 for k, i in enumerate(Iidx)}
print("crystal " + cube)
print("load " + cube + " id rho")
print("cprelax 1.0d-9 1.0d-9")
for a in range(len(Iidx)):
    for b in range(a + 1, len(Iidx)):
        ia, ib = Iidx[a], Iidx[b]
        d = frac[ia] - frac[ib]; d -= np.round(d)              # MIC
        dist = np.linalg.norm(d @ A) * BOHR
        if dist > dcut:
            continue
        mid = frac[ia] + d / 2.0
        print(f"  # I{Ilab[ia]}-I{Ilab[ib]} d={dist:.2f}")
        print(f"CP {mid[0]:.6f} {mid[1]:.6f} {mid[2]:.6f}")
