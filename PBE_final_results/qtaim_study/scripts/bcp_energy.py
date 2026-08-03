"""Estimate per-bond interaction energies from BCP density descriptors.

Those four literature correlations (Espinosa-style + iodine-specific) need the
LOCAL ENERGY DENSITIES at the BCP:
    V(r)  potential energy density  [Hartree/bohr^3]
    G(r)  kinetic   energy density  [Hartree/bohr^3]
critic2/bcp.out gives only rho and lap(=nabla^2 rho), so V and G are recovered
via the Abramov kinetic-energy-density approximation + the local virial theorem
(the SAME route critic2/Multiwfn use internally when only a rho cube exists):

    G = (3/10)(3 pi^2)^(2/3) rho^(5/3) + (1/6) nabla^2 rho          [a.u.]
    V = (1/4) nabla^2 rho - 2 G                                     [a.u.]
    H = G + V           (total energy density; sign -> covalent/ionic)

All quantities in ATOMIC UNITS (rho in e/bohr^3, nabla^2 rho in e/bohr^5).
Input lap may be given in e/Å^5 -- set LAP_IN_ANG=True and it is converted by
*(0.5291772109)^5 = 0.041494.

Correlations (E in Hartree when the property is in a.u.; x2625.5 -> kJ/mol):
    E1 = -V/2          Espinosa-Molins-Lecomte  (general / H-bonds)
    E2 = 0.429*G       Abramov/G-based           (general / H-bonds)
    E3 = 0.68*(-V)     V-based, iodine-specific
    E4 = 0.67*G        G-based, iodine-specific
"""
import numpy as np

BOHR = 0.5291772109
ANG5_TO_BOHR5 = BOHR**5                 # 0.041494  (e/Å^5 -> e/bohr^5)
HA_TO_KJMOL = 2625.4996                 # Hartree -> kJ/mol
TF = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)   # 2.87124  Thomas-Fermi prefactor


def energy_densities(rho_au, lap_au):
    """Abramov G, virial V, total H (a.u.) from rho & lap (a.u.)."""
    g = TF * rho_au**(5.0 / 3.0) + lap_au / 6.0
    v = lap_au / 4.0 - 2.0 * g
    return g, v, g + v


def e_int(rho_au, lap_au):
    """Four literature interaction-energy estimates, in kJ/mol."""
    g, v, h = energy_densities(rho_au, lap_au)
    return {
        "E_ML=-V/2":  -0.5 * v,          # Espinosa-Molins-Lecomte
        "E=0.429G":    0.429 * g,        # general G-based
        "E=0.68|V|":   0.68 * (-v),      # iodine V-based
        "E=0.67G":     0.67 * g,         # iodine G-based
    }


def report(rows, lap_in_ang=True):
    """rows = list of (label, rho_e/bohr3, lap_value)."""
    conv = ANG5_TO_BOHR5 if lap_in_ang else 1.0
    print(f"{'bond':16s}{'rho':>8s}{'lap(au)':>9s}{'G':>8s}{'V':>9s}"
          f"{'H':>8s}{'-V/2':>7s}{'.429G':>7s}{'.68|V|':>7s}{'.67G':>7s}")
    print("-" * 86)
    for lab, rho, lap in rows:
        lap_au = lap * conv
        g, v, h = energy_densities(rho, lap_au)
        e = e_int(rho, lap_au)
        print(f"{lab:16s}{rho:8.4f}{lap_au:9.4f}{g:8.4f}{v:9.4f}{h:8.4f}"
              f"{e['E_ML=-V/2']*HA_TO_KJMOL:7.1f}{e['E=0.429G']*HA_TO_KJMOL:7.1f}"
              f"{e['E=0.68|V|']*HA_TO_KJMOL:7.1f}{e['E=0.67G']*HA_TO_KJMOL:7.1f}")
    print("-" * 86)
    print("energies in kJ/mol; rho/G/V/H in a.u. (rho e/bohr3, G V H Hartree/bohr3)")


if __name__ == "__main__":
    # Representative bonds from the user's tables: (label, rho e/bohr3, lap e/Å^5)
    demo = [
        ("Y_act I5-I6",  0.0695, -0.211),   # covalent-ish I2, 2.74 Å, shared-shell
        ("Y_act I3-I4",  0.0533, +0.459),   # polyiodide, 2.89 Å
        ("Y_solv I7-I8", 0.0493, +0.602),   # polyiodide, 2.93 Å
        ("Y_vc I5-I6",   0.0213, +0.762),   # long weak I-I, 3.41 Å
        ("Y_solv I5-O",  0.0245, +1.317),   # I..O halogen bond, 2.87 Å
        ("Y_act I7-O",   0.0104, +0.776),   # weak I..O, 3.32 Å
    ]
    report(demo, lap_in_ang=True)
