# Raman discussion (manuscript draft)

Folds the band→contact assignment (Table `tab:raman_assignment`,
`raman_band_assignment.tex`; Fig. `raman_Y_st_act_vc`) into the argument.
Tone per author guidance: no "fingerprint", no "falsifiable predictions";
explain why each band is where it is, especially why vc is highest.

v2: only chemically bonded I–I pairs (≤3.25 Å, the coordination gap) enter the
polarizability model; non-bonded vdw contacts excluded (the st 103 cm⁻¹
"cross-link" band was an artifact of the earlier 4.0 Å cutoff and is removed).

---

We computed the iodine-sublattice vibrational spectrum by finite differences
(central differences, Δ = 0.01 Å, on the 8 I atoms with the Y–O framework held
fixed; PBE+D3/PAW). Applying the same protocol to gas-phase I₂ returns the
experimental stretching frequency (214.5 vs 213 cm⁻¹), so the MOF frequencies
are reported on an absolute scale. Every calculated band lies below the
gas-phase I₂ line (Fig. `raman_Y_st_act_vc`), within the polyiodide range
(I₃⁻ ≈ 110, I₅⁻ ≈ 145–170 cm⁻¹), consistent with the charge-transfer /
polyiodide picture established by the PDOS and bond-critical-point analysis.

To connect each band to its structural origin we projected every normal mode
onto its constituent I–I bond-stretch coordinates and weighted the contacts by
a bond-polarizability derivative (dα/dr ∝ r⁻⁴, r₀ = 2.67 Å for gas I₂). Only
chemically bonded pairs enter this sum: the I–I distance spectrum shows a clean
coordination gap with no contact between 3.12 and 3.40 Å, so we cut at 3.25 Å
and exclude the longer through-space contacts (e.g. the 3.6 Å
branching-iodine···in-chain distance), which carry no bond polarizability —
including them would leak spurious intensity into chain-bending modes. The
frequencies are Hessian eigenvalues; the polarizability model governs only the
intensities and the band→contact split, so the assignment is robust to the
model while the per-contact percentages are semi-quantitative
(Table `tab:raman_assignment`). The three motifs differ structurally and their
spectra follow those differences directly. In Y-MOF_act the intense 150–152 cm⁻¹
band is dominated (82 %) by the short 2.74 Å I₂-like pairs, red-shifted from
free I₂ by coupling to the neighbouring chain; its 160 cm⁻¹ band belongs to the
2.89 Å chain bonds and its 144–145 cm⁻¹ band to the longer 3.0 Å secondary
bonds. In Y-MOF_st the principal 152–156 cm⁻¹ band is a delocalised stretch
over the 2.93–2.99 Å chain (no single dominant contact), with a second band at
126 cm⁻¹ from the 2.99 Å contacts.

Y-MOF_vc carries the highest-frequency band (173–176 cm⁻¹), and the assignment
shows why: this band is entirely (100 %) the stretch of the uniform ~2.94 Å
chain — every contributing contact is the same length. A uniform chain has no
soft or elongated bonds to localise modes at lower frequency, so its stretching
modes are maximally coupled and near-degenerate and accumulate at the top of the
range. This is not a shortest-bond effect: Y-MOF_act actually contains the
shortest contacts (2.74 Å) yet its bands sit lower, because its Raman intensity
is distributed over a spread of bond lengths (2.74–3.04 Å) whose collective
modes vibrate at lower frequencies. The ordering vc > st ≈ act therefore
reflects the uniformity of the iodine chain rather than the shortest I–I
distance, and traces the same motif dependence seen in the charge-transfer and
PDOS trends.
