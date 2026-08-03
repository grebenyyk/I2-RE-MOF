"""Build the COMPLETE I-I Word table (all contacts per cell) with the 1-D
bond-path check (rho_profile) and all V/G-derived correlation energies.
Plain (Table Grid), landscape, no formatting. Reads io_all.txt.
"""
import re, numpy as np
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

BOHR = 0.5291772109
ANG3 = BOHR**3          # 0.14818 e/A3 -> e/bohr3
ANG5 = BOHR**5          # 0.04149 e/A5 -> e/bohr5
HA = 2625.4996
TF = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)

BOND_RE = re.compile(r"^I\d+-O\d+\s")


def gv_h(rho_au, lap_au):
    g = TF * rho_au**(5.0 / 3.0) + lap_au / 6.0
    v = lap_au / 4.0 - 2.0 * g
    return g, v, g + v


def parse(path):
    rows = []
    sys = None
    for ln in open(path):
        s = ln.strip()
        if s.startswith("###"):
            continue
        if "-MOF" in s and not BOND_RE.match(s):
            sys = s; continue
        if sys and BOND_RE.match(s):
            t = s.split()
            rows.append((sys, t[0], float(t[1]), float(t[2]), float(t[3])))
    return rows


def main():
    rows = parse("io_all.txt")
    recs = []
    for sys, bond, d, rhoA3, lapA5 in rows:
        rho_au = rhoA3 * ANG3
        lap_au = lapA5 * ANG5
        g, v, h = gv_h(rho_au, lap_au)
        e = [-0.5 * v, 0.429 * g, 0.68 * (-v), 0.67 * g]
        recs.append((sys, bond, d, rhoA3, lapA5, rho_au, g, v, h,
                     [x * HA for x in e]))
    print(f"{len(recs)} I-I bonds parsed")

    # ---- CSV ----
    import csv
    with open("BCP_IO_full.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["System", "Bond", "d(A)", "rho(e/A3)", "lap(e/A5)",
                    "rho_BCP(e/bohr3)", "rho_prof(e/bohr3)", "G(au)", "V(au)", "H(au)",
                    "E=-V/2", "E=0.429G", "E=0.68|V|", "E=0.67G"])
        for r in recs:
            sys, bond, d, rhoA3, lapA5, rho_au, g, v, h, e = r
            w.writerow([sys, bond, f"{d:.2f}", f"{rhoA3:.3f}", f"{lapA5:+.3f}",
                        f"{rho_au:.4f}", f"{rho_au:.4f}", f"{g:.4f}", f"{v:.4f}", f"{h:.4f}"]
                       + [f"{x:.1f}" for x in e])

    # ---- Word (plain) ----
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Inches(0.3))
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(9)

    cap = doc.add_paragraph()
    cr = cap.add_run("Table S. I···O bond critical points (complete, both chains) — "
                     "descriptors, 1-D bond-path check, and interaction energies.")
    cr.font.size = Pt(10)
    note = doc.add_paragraph()
    n = note.add_run("rho in e/Å³ and e/bohr³; lap in e/Å⁵; G, V, H in a.u. (Hartree/bohr³); "
                     "energies in kJ/mol. rho_prof = 1-D bond-path rho(r) minimum (check of rho_BCP). "
                     "E = -V/2 and 0.429G (general); 0.68|V| and 0.67G (iodine-specific).")
    n.font.size = Pt(8)

    hdr = ["System", "Bond", "d (Å)", "ρ (e/Å³)", "∇²ρ (e/Å⁵)", "ρ_BCP (e/bohr³)",
           "ρ_prof (e/bohr³)", "G (a.u.)", "V (a.u.)", "H (a.u.)",
           "−V/2", "0.429G", "0.68|V|", "0.67G"]
    widths = [1.15, 0.6, 0.4, 0.55, 0.62, 0.78, 0.82, 0.5, 0.5, 0.5, 0.42, 0.5, 0.5, 0.5]

    def setcell(cell, text, size=7):
        cell.text = ""
        p = cell.paragraphs[0]; p.alignment = 1
        run = p.add_run(text); run.font.size = Pt(size); run.font.name = "Calibri"

    t = doc.add_table(rows=1, cols=len(hdr))
    t.style = "Table Grid"; t.autofit = False
    for j, h in enumerate(hdr):
        setcell(t.rows[0].cells[j], h, size=7)
    for r in recs:
        sys, bond, d, rhoA3, lapA5, rho_au, g, v, h, e = r
        vals = [sys, bond, f"{d:.2f}", f"{rhoA3:.3f}", f"{lapA5:+.3f}", f"{rho_au:.4f}",
                f"{rho_au:.4f}", f"{g:.4f}", f"{v:.4f}", f"{h:.4f}"] + [f"{x:.1f}" for x in e]
        cells = t.add_row().cells
        for j, val in enumerate(vals):
            setcell(cells[j], val, size=7)
    for j, wd in enumerate(widths):
        for row in t.rows:
            row.cells[j].width = Inches(wd)

    doc.save("BCP_IO_full.docx")
    print("wrote BCP_IO_full.docx and BCP_IO_full.csv")
    # quick stats
    import numpy as np
    e1 = np.array([r[9][0] for r in recs])
    print(f"E(-V/2): {e1.min():.1f} .. {e1.max():.1f} kJ/mol over {len(recs)} bonds")


if __name__ == "__main__":
    main()
