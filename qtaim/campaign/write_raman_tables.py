#!/usr/bin/env python3
"""Build manuscript-ready Raman band->contact assignment tables for the 3 Y-MOF
variants, from the saved eigenvectors (evecs.npy) + reference I geometry
(I_frac_ref.npy, cell.npy). NO new SCFs: this just un-sums the bond-polarizability
activity that raman_driver.py already computes as a per-contact sum.

Outputs (in final_activated_results/):
  raman_band_assignment.csv   full per-mode per-contact fractions (SI / data)
  raman_contact_geometry.csv  I-I contacts <4 A per variant, length + motif
  raman_band_assignment.md    condensed peak table + geometry table (manuscript)
  raman_band_assignment.tex   LaTeX condensed peak table

Activity model (identical to raman_driver.py:529):
  dadr_b = (r0/r_b)^4 , r0 = 2.67 A (gas I2)
  contrib_{b,k} = (dadr_b * s_b)^2 ,  s_b = (e_j - e_i).u_hat / sqrt(2)
  band activity_k = sum_b contrib_{b,k}
  fraction f_{b,k} = contrib_{b,k} / activity_k  (this table)
"""
import os, numpy as np

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
OUT  = f"{BASE}/final_activated_results"
os.makedirs(OUT, exist_ok=True)
R0, P = 2.67, 4.0
CLUSTER_CM = 4.0          # merge modes within this into one spectral peak
MIN_PEAK_ACT = 10.0       # CSV: keep a mode if its stick activity >= this
MIN_BROADENED = 5.0       # table: report a band if broadened height >= this %
MIN_REPORT_FRAC = 4.0     # report a contact in a peak only if >= this %

VAR = [("Y_str","st","Y-MOF$_{\\rm st}$","Y-MOFst"),
       ("Y_actc","act","Y-MOF$_{\\rm act}$","Y-MOFact"),
       ("Y_deep","vc","Y-MOF$_{\\rm vc}$","Y-MOFvc")]

# Bond vs van-der-Waals split. The pooled I-I distance spectrum has a clean
# 0.28 A coordination gap with NOTHING between 3.121 and 3.402 A, so any cutoff
# in (3.12, 3.40) is objective. Contacts >= this are through-space / vdw
# (e.g. st's 3.61 A branching-I1 ... in-chain-I4) and carry no bond
# polarizability -> excluded from the Raman intensity model. Listed in the
# geometry table for completeness only.
BOND_CUTOFF = 3.25

def motif(r):
    if r < 2.80: return "I2-like pair"
    if r < 3.00: return "chain bond"
    if r < BOND_CUTOFF: return "secondary bond"
    return "vdw contact (non-bonded)"

def load(sys, cutoff=BOND_CUTOFF):
    r = f"{BASE}/{sys}/raman"
    freqs=np.load(f"{r}/freqs.npy"); evecs=np.load(f"{r}/evecs.npy")
    act=np.load(f"{r}/raman_activities.npy")
    cell=np.load(f"{r}/cell.npy"); If=np.load(f"{r}/I_frac_ref.npy")
    bonds=[]
    for i in range(8):
        for j in range(i+1,8):
            df=If[i]-If[j]; df-=np.round(df)
            rr=np.linalg.norm(df@cell)
            if rr<=cutoff: bonds.append((i+1,j+1,rr,(df@cell)/rr))
    return freqs,evecs,act,bonds

def per_mode_frac(freqs,evecs,act,bonds):
    """Recompute per-mode Raman activity AND per-contact fractions from the
    (filtered) bond set, so the displayed activities reflect the BOND_CUTOFF
    (the saved raman_activities.npy used a 4.0 A cutoff and would keep the
    spurious vdw contributions). Activity normalized to per-variant max=100.
    Returns list of (freq, act, {contact(i,j,r): frac%})."""
    dadr=np.array([(R0/r)**P for (_,_,r,_) in bonds])
    raw=[]
    for ix in range(len(freqs)):
        e=evecs[:,ix]; tot=0.0; s={}
        for kb,(i,j,r,u) in enumerate(bonds):
            di=e[3*(i-1):3*(i-1)+3]; dj=e[3*(j-1):3*(j-1)+3]
            amp=(dadr[kb]*np.dot(dj-di,u)/np.sqrt(2))**2
            s[(i,j,r)]=amp; tot+=amp
        raw.append((freqs[ix], tot, s))
    mx=max((t for _,t,_ in raw),default=0.0) or 1.0
    out=[]
    for f,t,s in raw:
        a=100.0*t/mx
        frac={k:(v/t*100.0 if t>0 else 0.0) for k,v in s.items()}
        out.append((f,a,frac))
    return out

def cluster_peaks(rows):
    """Group modes (freq desc) within CLUSTER_CM into peaks. Return list of
    (freq_lo, freq_hi, centroid, stick_max, n_modes, merged_frac)."""
    sig=[r for r in rows if r[0]>30]                      # drop imag/far
    sig=sorted(sig,key=lambda r:-r[0])
    peaks=[]; cur=[]
    for r in sig:
        if cur and abs(cur[-1][0]-r[0])<=CLUSTER_CM: cur.append(r)
        else:
            if cur: peaks.append(cur)
            cur=[r]
    if cur: peaks.append(cur)
    merged=[]
    for grp in peaks:
        ws=sum(rr[1] for rr in grp)                       # total stick activity
        cf={}
        for f,a,frac in grp:
            for k,v in frac.items(): cf[k]=cf.get(k,0.0)+a*v
        cf={k:v/ws for k,v in cf.items()}                 # activity-weighted
        flo=min(rr[0] for rr in grp); fhi=max(rr[0] for rr in grp)
        cent=sum(rr[0]*rr[1] for rr in grp)/ws            # activity-weighted centroid
        amax=max(rr[1] for rr in grp)
        merged.append((flo,fhi,cent,amax,len(grp),cf))
    return merged

# Lorentzian broadening, identical to plot_raman_Y_motifs.py (W=5 cm-1). The
# broadened curve is what the spectrum shows; degenerate sticks add up, so a
# 4-mode cluster (e.g. vc 173-176) towers over a single mode of the same stick
# activity (e.g. vc 159). We report the broadened peak height (apparent
# intensity in the spectrum), not the bare per-mode stick, so the table matches
# the figure. n_modes makes the degeneracy explicit.
W_BR = 5.0
def broadened_curve(rows, xmin=30.0, xmax=300.0):
    f=np.array([r[0] for r in rows]); a=np.array([r[1] for r in rows])
    x=np.linspace(xmin,xmax,int(9*(xmax-xmin)))
    y=np.zeros_like(x)
    for fi,ai in zip(f,a): y+=ai*W_BR*W_BR/((x-fi)**2+W_BR*W_BR)
    return x, y/(y.max() or 1.0)

# ---------------- write full CSV ----------------
with open(f"{OUT}/raman_band_assignment.csv","w") as fcsv:
    fcsv.write("variant,freq_cm1,rel_activity,contact,bond_ang,motif,pct_of_band\n")
    for sys,key,_,short in VAR:
        freqs,evecs,act,bonds=load(sys)
        rows=per_mode_frac(freqs,evecs,act,bonds)
        for fr,a,frac in sorted(rows,key=lambda r:-r[0]):
            if fr<30 or a<1.0: continue
            for (i,j,r),v in sorted(frac.items(),key=lambda x:-x[1]):
                if v<MIN_REPORT_FRAC: continue
                fcsv.write(f"{short},{fr:.1f},{a:.1f},I{i:02d}-I{j:02d},{r:.2f},"
                           f"{motif(r).replace(' ','_')},{v:.1f}\n")
print("wrote", f"{OUT}/raman_band_assignment.csv")

# ---------------- contact geometry CSV ----------------
with open(f"{OUT}/raman_contact_geometry.csv","w") as fg:
    fg.write("variant,contact,bond_ang,motif\n")
    for sys,_,_,short in VAR:
        _,_,_,bonds=load(sys, cutoff=4.0)   # all contacts incl vdw, for completeness
        for i,j,r,_ in sorted(bonds,key=lambda b:b[2]):
            fg.write(f"{short},I{i:02d}-I{j:02d},{r:.2f},{motif(r).replace(' ','_')}\n")
print("wrote", f"{OUT}/raman_contact_geometry.csv")

# ---------------- condensed peak table (MD + LaTeX) ----------------
def peak_desc(cf, fmt="md"):
    """Lead with motif-binned fraction (the headline), then top contacts.
    cf values are percentages (sum ~100). fmt: 'md' or 'tex'."""
    mf={}
    for (i,j,r),v in cf.items(): mf[motif(r)]=mf.get(motif(r),0.0)+v
    lead_m,lead_p=max(mf.items(),key=lambda x:x[1])
    items=sorted(cf.items(),key=lambda x:-x[1])[:3]
    parts=[]
    for (i,j,r),v in items:
        if v<MIN_REPORT_FRAC: break
        if fmt=="tex":
            parts.append(rf"I{i:02d}--I{j:02d} ({r:.2f}\,\AA; {v:.0f}\%)")
        else:
            parts.append(f"I{i:02d}-I{j:02d} ({r:.2f} Å; {v:.0f}%)")
    contacts=", ".join(parts) if parts else "delocalized"
    if fmt=="tex":
        return rf"{lead_m} ({lead_p:.0f}\%): {contacts}", lead_p
    return f"{lead_m} ({lead_p:.0f}%): {contacts}", lead_p

md=[]; md.append("# Raman band → I–I contact assignment (I sublattice, PBE+D3)\n")
md.append("Frequencies from the finite-difference Hessian (framework frozen, "
          "Δ=0.01 Å). Relative activities and contact percentages from the "
          "bond-polarizability model (dα/dr ∝ r⁻⁴, r₀=2.67 Å gas I₂), applied "
          "only to chemically bonded I–I pairs (≤3.25 Å — the natural "
          "coordination gap; see geometry table). Non-bonded van-der-Waals "
          "contacts (≥3.25 Å, e.g. the 3.61 Å branching-I···in-chain distance) "
          "carry no bond polarizability and are excluded: they would otherwise "
          "leak spurious intensity into chain-bending modes. Mode→contact "
          "fractions are read off the saved eigenvectors; the spectrum and this "
          "table share the same eigenvectors/model.\n")

# geometry table
md.append("\n## I–I contact geometry (< 4.0 Å, min-image). ≥3.25 Å = non-bonded vdw (excluded from Raman model)\n")
md.append("| variant | contact | r (Å) | motif |")
md.append("|---|---|---|---|")
for sys,key,lbl,short in VAR:
    _,_,_,bonds=load(sys, cutoff=4.0)   # all contacts incl vdw, listed for completeness
    for i,j,r,_ in sorted(bonds,key=lambda b:b[2]):
        md.append(f"| {short} | I{i:02d}-I{j:02d} | {r:.2f} | {motif(r)} |")

# band table
md.append("\n## Band assignment (spectral peaks; broadened intensity + motif)\n")
md.append("`rel. int.` = broadened peak height as a % of each variant's tallest "
          "peak (i.e. exactly the height in Fig. `raman_Y_st_act_vc`). `n` = "
          "number of near-degenerate modes in the peak — a 4-mode cluster towers "
          "over a single mode of equal per-mode activity.\n")
md.append("| variant | ν (cm⁻¹) | n | rel. int. | dominant motif + top contacts |")
md.append("|---|---|---|---|---|")
tex=[]; tex.append(r"\begin{table}[ht]")
tex.append(r"\centering")
tex.append(r"\caption{Assignment of calculated Raman bands of the I sublattice "
           r"(PBE+D3, finite-difference Hessian, framework frozen) to I--I "
           r"contacts for the three Y-MOF variants. Frequencies from the "
           r"Hessian eigenvalues; per-contact percentages from a "
           r"bond-polarizability model ($d\alpha/dr\propto r^{-4}$, "
           r"$r_0=2.67$~\AA\ gas I$_2$), summed only over chemically bonded "
           r"I--I pairs ($\le 3.25$~\AA\ -- the coordination gap; non-bonded "
           r"van-der-Waals contacts excluded). Near-degenerate modes within "
           r"4~cm$^{-1}$ are merged. \emph{rel.\ int.} is the Lorentzian-broadened "
           r"peak height ($w=5$~cm$^{-1}$, as plotted) as a \% of each variant's "
           r"tallest peak, NOT the bare per-mode stick activity; $n$ is the "
           r"number of degenerate modes in the peak, so a 4-mode cluster "
           r"(e.g.\ vc 173--176) rises well above a single mode of equal stick "
           r"activity (e.g.\ vc 159). The descriptor leads with the dominant "
           r"contact motif (fraction of the band's intensity) and lists the top "
           r"individual contacts (length; \% of band).}")
tex.append(r"\label{tab:raman_assignment}")
tex.append(r"\begin{tabular}{llccl}")
tex.append(r"\hline")
tex.append(r"variant & $\tilde\nu$ (cm$^{-1}$) & $n$ & rel.\ int. & dominant motif + top contacts ($r$/\AA; \% of band) \\")
tex.append(r"\hline")
for sys,key,lbl,short in VAR:
    freqs,evecs,act,bonds=load(sys)
    rows=per_mode_frac(freqs,evecs,act,bonds)
    x,yn=broadened_curve(rows)
    peaks=cluster_peaks(rows)
    out=[]
    for flo,fhi,cent,amax,nmodes,cf in peaks:
        h=float(np.interp(cent,x,yn))*100.0
        if h<MIN_BROADENED: continue
        rng = f"{fhi:.0f}" if abs(fhi-flo)<1 else f"{flo:.0f}--{fhi:.0f}"
        out.append((cent,rng,nmodes,h,cf))
    for cent,rng,nmodes,h,cf in sorted(out,key=lambda z:-z[0]):   # high freq first
        desc_md,_=peak_desc(cf,"md"); desc_tex,_=peak_desc(cf,"tex")
        md.append(f"| {short} | {rng} | {nmodes} | {h:.0f} | {desc_md} |")
        tex.append(rf"{short} & {rng} & {nmodes} & {h:.0f} & {desc_tex} \\")
tex.append(r"\hline")
tex.append(r"\end{tabular}")
tex.append(r"\end{table}")

with open(f"{OUT}/raman_band_assignment.md","w") as f:
    f.write("\n".join(md)+"\n")
print("wrote", f"{OUT}/raman_band_assignment.md")
with open(f"{OUT}/raman_band_assignment.tex","w") as f:
    f.write("\n".join(tex)+"\n")
print("wrote", f"{OUT}/raman_band_assignment.tex")
print("\ndone.")
