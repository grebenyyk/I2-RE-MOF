"""Build the QTAIM BCP energy Word table from tables_input.txt.

Parses the two user tables (I-I, I-O) verbatim, adds Abramov-derived local
energy densities G, V, H (a.u.) and the four literature interaction-energy
estimates (kJ/mol), and writes a PLAIN (unformatted, Table-Grid) .docx plus a
.csv for convenience. ρ is taken in e/bohr^3 (a.u.); ∇²ρ in e/Å^5 is converted
to a.u. An internal check flags any row where ρ(e/bohr^3) != ρ(e/Å^3)*0.1482.
"""
import re
import csv
import numpy as np
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

BOHR = 0.5291772109
ANG3_TO_BOHR3 = BOHR**3        # 0.14818  (e/Å^3 -> e/bohr^3)
ANG5_TO_BOHR5 = BOHR**5        # 0.04150  (e/Å^5 -> e/bohr^5)
HA_TO_KJMOL = 2625.4996
TF = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)   # 2.87124

BOND_RE = re.compile(r"^I\d+-[IO]\d+\s")


def gv_h(rho_au, lap_au):
    g = TF * rho_au**(5.0 / 3.0) + lap_au / 6.0
    v = lap_au / 4.0 - 2.0 * g
    return g, v, g + v


def e_int(rho_au, lap_au):
    g, v, _ = gv_h(rho_au, lap_au)
    return [-0.5 * v, 0.429 * g, 0.68 * (-v), 0.67 * g]


def parse(path):
    tables = {"I-I": [], "I-O": []}
    cur_tbl, cur_sys = None, None
    for ln in open(path):
        s = ln.strip()
        if s.startswith("###TABLE"):
            cur_tbl = "I-I" if "I-I" in s else "I-O"
            continue
        if "-MOF" in s and s and not BOND_RE.match(s):
            cur_sys = s
            continue
        if cur_tbl and cur_sys and BOND_RE.match(s):
            t = s.split()
            # bond d rhoA3 lapA5 rhoBCP_b3 rhoPROF_b3
            tables[cur_tbl].append({
                "sys": cur_sys, "bond": t[0], "d": float(t[1]),
                "rhoA3": float(t[2]), "lapA5": float(t[3]),
                "rho_b3": float(t[4]), "rho_prof": float(t[5]),
            })
    return tables


def main():
    tables = parse("tables_input.txt")

    # consistency check + compute
    warnings = []
    for kind, rows in tables.items():
        for r in rows:
            chk = r["rhoA3"] * ANG3_TO_BOHR3
            if abs(chk - r["rho_b3"]) > 0.0006:           # ~1% tolerance
                warnings.append(f"  {kind} {r['sys']} {r['bond']}: "
                                f"rho {r['rho_b3']} vs {chk:.4f} from e/A3")
            lap_au = r["lapA5"] * ANG5_TO_BOHR5
            g, v, h = gv_h(r["rho_b3"], lap_au)
            r.update(lap_au=lap_au, G=g, V=v, H=h)
            r["E"] = [e * HA_TO_KJMOL for e in e_int(r["rho_b3"], lap_au)]

    # ---- stdout TSV + stats for verification ----
    print("=== consistency warnings (should be empty) ===")
    print("\n".join(warnings) if warnings else "  none")
    for kind, rows in tables.items():
        es = np.array([r["E"] for r in rows])           # (n,4)
        print(f"\n=== {kind}: {len(rows)} bonds ===")
        print(f"  E(-V/2)   {es[:,0].min():.1f} .. {es[:,0].max():.1f} kJ/mol")
        print(f"  E(0.429G) {es[:,1].min():.1f} .. {es[:,1].max():.1f}")
        print(f"  E(0.68|V|){es[:,2].min():.1f} .. {es[:,2].max():.1f}")
        print(f"  E(0.67G)  {es[:,3].min():.1f} .. {es[:,3].max():.1f}")

    # ---- CSV ----
    cols = ["system", "bond", "d(A)", "rho(e/A3)", "lap(e/A5)",
            "rhoBCP(e/bohr3)", "rho_prof(e/bohr3)", "G(au)", "V(au)", "H(au)",
            "E=-V/2(kJ/mol)", "E=0.429G", "E=0.68|V|", "E=0.67G"]
    with open("BCP_energies.csv", "w", newline="") as f:
        w = csv.writer(f)
        for kind, rows in tables.items():
            w.writerow([f"### {kind} bonds ###"])
            w.writerow(cols)
            for r in rows:
                w.writerow([r["sys"], r["bond"], f'{r["d"]:.2f}', f'{r["rhoA3"]:.3f}',
                            f'{r["lapA5"]:+.3f}', f'{r["rho_b3"]:.4f}', f'{r["rho_prof"]:.4f}',
                            f'{r["G"]:.4f}', f'{r["V"]:.4f}', f'{r["H"]:.4f}']
                           + [f'{e:.1f}' for e in r["E"]])
            w.writerow([])
    print("\nwrote BCP_energies.csv")

    # ---- Word .docx (plain, landscape) ----
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Inches(0.3))

    # default font
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(9)

    h = ["System", "Bond", "d (Å)", "ρ (e/Å³)", "∇²ρ (e/Å⁵)",
         "ρ_BCP (e/bohr³)", "ρ_prof (e/bohr³)", "G (a.u.)", "V (a.u.)",
         "H (a.u.)", "−V/2", "0.429G", "0.68|V|", "0.67G"]
    # column widths (inches); numeric cols narrower
    widths = [1.45, 0.62, 0.42, 0.55, 0.62, 0.72, 0.78, 0.5, 0.5, 0.5,
              0.42, 0.5, 0.5, 0.5]

    def setcell(cell, text, size=7, hdr=False):
        cell.text = ""
        p = cell.paragraphs[0]
        p.alignment = 1  # center
        run = p.add_run(text)
        run.font.size = Pt(size)
        run.font.name = "Calibri"

    title = doc.add_paragraph()
    tr = title.add_run("QTAIM bond-critical-point descriptors and interaction-energy estimates")
    tr.bold = False
    tr.font.size = Pt(11)
    note = doc.add_paragraph()
    n = note.add_run(
        "Local energy densities from the Abramov approximation + local virial: "
        "G = (3/10)(3π²)^(2/3)ρ^(5/3) + (1/6)∇²ρ,  V = (1/4)∇²ρ − 2G,  H = G + V "
        "(a.u. = Hartree/bohr³; ρ in e/bohr³, ∇²ρ converted from e/Å⁵). "
        "Interaction energies: −V/2 and 0.429G (general, H-bond calibrated); "
        "0.68|V| and 0.67G (iodine-specific). All energies in kJ/mol.")
    n.font.size = Pt(8)

    for kind, rows in tables.items():
        cap = doc.add_paragraph()
        cr = cap.add_run(f"Table: {kind} bonds")
        cr.font.size = Pt(9)
        t = doc.add_table(rows=1, cols=len(h))
        t.style = "Table Grid"
        t.autofit = False
        for j, hd in enumerate(h):
            setcell(t.rows[0].cells[j], hd, size=7)
        for r in rows:
            cells = t.add_row().cells
            vals = [r["sys"], r["bond"], f'{r["d"]:.2f}', f'{r["rhoA3"]:.3f}',
                    f'{r["lapA5"]:+.3f}', f'{r["rho_b3"]:.4f}', f'{r["rho_prof"]:.4f}',
                    f'{r["G"]:.4f}', f'{r["V"]:.4f}', f'{r["H"]:.4f}'] \
                   + [f'{e:.1f}' for e in r["E"]]
            for j, val in enumerate(vals):
                setcell(cells[j], val, size=7)
        for j, wd in enumerate(widths):
            for row in t.rows:
                row.cells[j].width = Inches(wd)
        doc.add_paragraph()

    out = "BCP_energies.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
