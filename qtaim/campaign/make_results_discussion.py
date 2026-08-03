#!/usr/bin/env python3
"""Rebuild results_discussion_draft.docx with the corrected QTAIM understanding
(3c-4e strained polyiodide; asymmetric activated; ELF reconciled; I-O verified
with same rigor as I-I; EML bond energies dropped). Overwrites the draft (old
backed up as results_discussion_draft_prev.docx). Sub/superscripts via ~x~/^x^."""
import re
from docx import Document
from docx.shared import Pt
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"

doc=Document()
st=doc.styles["Normal"]; st.font.name="Liberation Sans"; st.font.size=Pt(10.5)

def add_runs(p,text):
    for seg in re.split(r"(~[^~]+~|\^[^^]+\^)",text):
        if not seg: continue
        if seg.startswith("~") and seg.endswith("~"):
            r=p.add_run(seg[1:-1]); r.font.subscript=True
        elif seg.startswith("^") and seg.endswith("^"):
            r=p.add_run(seg[1:-1]); r.font.superscript=True
        else:
            p.add_run(seg)

def H1(t):
    p=doc.add_paragraph(); p.style=doc.styles["Heading 1"]; r=p.add_run(t); r.font.name="Liberation Sans"
def H2(t):
    p=doc.add_paragraph(); p.style=doc.styles["Heading 2"]; r=p.add_run(t); r.font.name="Liberation Sans"
def P(t):
    p=doc.add_paragraph(); add_runs(p,t)

H1("Results and Discussion")

H2("Polyiodide chain topology (QTAIM bond critical points)")
P("Bond critical points (BCPs) were located on the valence charge density for the I~2~-loaded "
  "structures of Y-, Tb-, and Gd-MOFs in both the strained (guest-imposed) and activated "
  "(relaxed-empty) pore geometries. Every I–I contact along the adsorbed iodine chain carries a "
  "BCP with substantial electron density (ρ~BCP~ = 0.025–0.075 e·bohr^−3^, i.e. 0.17–0.50 e·Å^−3^), "
  "well above the van der Waals contact threshold (≲0.005 e·bohr^−3^) and confirming a connected, "
  "chemically bonded polyiodide rather than a chain of weakly interacting I~2~ molecules. The bond "
  "character, classified by the sign of the Laplacian, depends on pore geometry. In the strained "
  "pores every I–I BCP is closed-shell (∇^2^ρ~BCP~ > 0, ρ~BCP~ ≈ 0.037–0.049 e·bohr^−3^) at "
  "near-uniform distances of 2.9–3.1 Å. In the activated pores each chain contains one or two short "
  "shared-shell or borderline contacts (∇^2^ρ~BCP~ ≲ 0, ρ~BCP~ ≈ 0.07 e·bohr^−3^ ≈ 0.43–0.50 e·Å^−3^, "
  "d ≈ 2.69–2.79 Å) — the σ-bonded, near-gas-phase I~2~ cores — plus longer closed-shell links "
  "(∇^2^ρ~BCP~ > 0, ρ~BCP~ ≈ 0.04–0.05 e·bohr^−3^) at 2.9–3.1 Å that complete the chain. Within each "
  "system ρ~BCP~ decreases monotonically with distance, as required for a single atom pair.")
P("The uniformly closed-shell, near-symmetric strained chain is the QTAIM signature of a "
  "three-centre–four-electron (3c–4e) hypervalent interaction — the same bonding found in isolated "
  "I~3~^−^/I~5~^−^, whose I–I BCPs likewise carry ∇^2^ρ > 0 with ρ~BCP~ ≈ 0.05 e·bohr^−3^. The "
  "activated chain is the asymmetric [I–I···I]^−^ limit, retaining localized I~2~ cores alongside "
  "weaker donor–acceptor links. Pore geometry thus does not switch the polyiodide on or off; it "
  "tunes how far the adsorbate has progressed from neutral, molecular I~2~ toward a fully "
  "delocalized polyiodide.")
P("BCPs were located with critic2 on the PAW-reconstructed valence density (plot_num = 17). A "
  "small number of intra-atomic (3,−1) critical points in iodine’s 5p-shell region — identifiable "
  "by their gradient paths both terminating at the same nucleus (or at non-nuclear maxima) and by "
  "ρ ≈ 0.10 e·bohr^−3^ outside the bonding range — were excluded as topological artefacts of the "
  "PAW reconstruction. The reported ρ~BCP~ and ∇^2^ρ values are reproduced on a clean, "
  "PAW-augmentation-free valence field (plot_num = 0) and by an independent bond-path ρ(r) profile, "
  "confirming they are intrinsic to the bonds (see Supporting Information).")

H2("Electron localization function (ELF)")
P("The ELF is fully consistent with the QTAIM picture and, importantly, is not in conflict with "
  "the closed-shell Laplacian of the strained chain. Continuous localization basins (ELF > 0.5) "
  "run along every I–I contact in both pore geometries, with the localization graded — strongest "
  "(ELF ≈ 0.9) on the short σ-bonded I~2~ cores of the activated chains and weaker but "
  "unambiguously present (ELF ≈ 0.65–0.75) along the closed-shell links and across the entire "
  "strained chain. The coexistence of ELF > 0.5 with ∇^2^ρ~BCP~ > 0 is the expected signature of a "
  "delocalized 3c–4e / donor–acceptor bond: ELF > 0.5 reports a localized electron pair (a "
  "Pauli-repulsion deficit) between the atoms, whereas a positive Laplacian reports that charge is "
  "locally concentrated into the atomic basins rather than at the bond critical point. The two "
  "descriptors measure different quantities, and their simultaneous satisfaction is precisely what "
  "distinguishes hypervalent and charge-transfer bonding from both fully covalent (ELF > 0.5, "
  "∇^2^ρ < 0) and pure van der Waals (ELF < 0.5, ∇^2^ρ > 0) interactions.")
P("No ELF basin (ELF > 0.5) appears between iodine and the framework oxygens, consistent with the "
  "closed-shell (∇^2^ρ > 0) nature of the I···O contacts: charge transfer through these contacts "
  "proceeds via a donor–acceptor mechanism (n~O~ → σ*~I~), not through a shared electron pair.")

H2("Charge transfer (Bader differential density)")
P("The charge-transfer density δρ = ρ(MOF+I~2~) − ρ(MOF) − ρ(chain ref), integrated over Bader "
  "basins (Σδρ = 0 to numerical precision, confirming exact partitioning), quantifies the "
  "electronic reorganization upon adsorption. Iodine gains 0.3–0.9 e^−^ per unit cell across the "
  "six systems, with the charge donated almost entirely by the organic linker: carboxylate O "
  "(−0.16 to −0.45 e^−^) and linker C (−0.10 to −0.37 e^−^). The rare-earth cation is a "
  "charge-transfer spectator (|δq~M~| < 0.04 e^−^ per atom) for Y, Tb, and Gd alike. Strained-pore "
  "geometries transfer systematically more charge than activated ones (e.g. Y: +0.74 vs +0.50 "
  "e^−^/cell; Gd: +0.93 vs +0.32 e^−^/cell), tracking the stronger strained-pore sorption energies.")
P("This trend is mirrored in the bond topology: the more strongly reduced strained polyiodide is "
  "the fully delocalized, near-symmetric 3c–4e form (all I–I bonds equalized at ~3 Å, "
  "∇^2^ρ > 0 throughout), whereas the less-reduced activated polyiodide retains localized I~2~-like "
  "σ-bonded cores. The degree of charge transfer thus tracks how far the adsorbate has progressed "
  "from neutral, molecular I~2~ toward a delocalized polyiodide — pore geometry, not metal "
  "identity, is the controlling variable across the isostructural rare-earth series.")

H2("Iodine–framework anchoring (I···O bond critical points)")
P("Bond critical points between iodine and the framework, located with the same topological "
  "rigour applied to the I–I bonds (clean-field re-location, rejection of intra-atomic 5p-shell "
  "and near-nuclear artefacts, and confirmation by an independent bond-path profile), occur "
  "predominantly at the linker carboxylate (–COO^−^) oxygens, with a single μ~3~-OH (node "
  "hydroxyl) contact in Gd-activated. These I···O contacts are all closed-shell (∇^2^ρ > 0) with "
  "ρ~BCP~ ≈ 0.02–0.31 e·Å^−3^ (0.003–0.046 e·bohr^−3^) over d ≈ 2.5–4.0 Å, decreasing "
  "monotonically with distance, characteristic of halogen-bond / n~O~ → σ*~I~ donor–acceptor "
  "interactions. The activated pore presents more I···O contacts (e.g. Y: 8 vs 3), but each is "
  "individually weaker, and the net charge transfer is lower — more anchor points do not equate to "
  "more transfer, because the transfer is gated by the iodine acceptor level rather than by "
  "contact count.")

H2("Orbital mechanism of charge transfer (PDOS)")
P("Projected densities of states resolve the orbital mechanism of the O → I charge transfer. The "
  "iodine 5p manifold presents a low-lying σ* acceptor peak immediately above the Fermi level "
  "(E − E~F~ ≈ 0.1–0.3 eV), well separated from the higher 5p antibonding states at 2–3 eV. Its "
  "proximity to E~F~ identifies iodine as the readily reducible charge-transfer acceptor. In the "
  "strained pore this σ* acceptor lies lower (0.13 eV in Y) and is denser than in the activated "
  "pore (0.24 eV), so charge donation from the carboxylate O-2p valence states is more "
  "favourable — consistent with the larger Bader charge transfer in the strained pore and with its "
  "more delocalized 3c–4e polyiodide character.")
P("The rare-earth d states (Y-4d, Tb-5d, Gd-5d) peak at +4 to +5 eV above E~F~ — deep in the "
  "conduction band — with negligible density (<6% of peak) near E~F~. The 4f states (Tb, Gd) are "
  "in the pseudopotential core and even more deeply localized. The absence of metal d/f states "
  "near E~F~ rules out a through-metal charge-transfer pathway (O → RE → I): the RE^3+^ cations "
  "are redox-inert hard Lewis acids with no accessible electronic states at the energy of the "
  "O-2p donor or the I-5p σ* acceptor. Charge therefore transfers directly from the carboxylate "
  "lone pair to the iodine σ* antibonding orbital, bypassing the metal entirely.")

H2("Summary: the charge-transfer polyiodide picture")
P("The convergent picture from QTAIM, ELF, Bader charge analysis, and PDOS is as follows. The "
  "I~2~ adsorbate forms a polyiodide chain within the MOF pore whose bonding character is tuned by "
  "pore geometry. In the strained (more strongly reducing) pore the chain is a near-symmetric, "
  "delocalized I~3~^−^/I~5~^−^-type species with three-centre–four-electron I–I bonds (uniform "
  "~3 Å contacts, ∇^2^ρ > 0 throughout, continuous ELF > 0.5 basins); in the activated pore it "
  "adopts the asymmetric [I–I···I]^−^ form, retaining short σ-bonded I~2~ cores alongside weaker "
  "closed-shell links. The chain is anchored to the framework at carboxylate oxygens via "
  "closed-shell I···O (n~O~ → σ*~I~) contacts, verified to be genuine BCPs (ρ~BCP~ confirmed by "
  "an independent bond-path profile). Charge transfers directly from the carboxylate lone pair "
  "(O-2p, at E~F~) into the iodine σ* acceptor (I-5p, just above E~F~), as confirmed by PDOS "
  "energy overlap and Bader basin integration. The rare-earth cation does not participate "
  "electronically — it is redox-inert with d/f states far from E~F~ — but serves a structural "
  "role by coordinating the carboxylate oxygens into the geometry that enables the I···O contact. "
  "The degree of charge transfer is gated by the iodine 5p σ* acceptor level, which is lower and "
  "denser in the strained pore, explaining both the greater charge transfer and the more "
  "delocalized (3c–4e) polyiodide character there.")

doc.save(f"{BASE}/results_discussion_draft.docx")
print("wrote results_discussion_draft.docx (old backed up as results_discussion_draft_prev.docx)")
