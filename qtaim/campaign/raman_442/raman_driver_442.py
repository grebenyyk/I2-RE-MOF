#!/usr/bin/env python3
"""[442 variant] Adapted 2026-11-10 from the Gamma-era driver (qtaim/campaign/raman_driver.py,
gas-I2-validated): k-mesh 4x4x2, degauss 0.02 FD set (whole-set 0.03 fallback if deadlocks),
mixed-degauss guard, per-SCF timeout 7200 s. RY_TO_CM 3634.77 and BOND_CUTOFF 3.25 unchanged.
Requires <sys>/scf_loaded.in at the tight re-relaxed geometry (prefix raman442_<sys>_loaded).
Run from the raman_442 dir. RAMAN_GPUS=<gpu> confines a system to one GPU (two systems parallel).
"""
"""
Finite-difference phonon + Raman driver for the I sublattice of I2@Y-MOF.

Strategy
--------
- 8 I atoms (QE indices 73-80, 1-based) x 3 (x,y,z) x 2 (+/-) = 48 displaced SCFs per system.
- Framework atoms (1-72) HELD FIXED: never displaced.
- Each displaced SCF has its OWN outdir (./tmp_disp_<tag>) seeded from tmp_loaded/<prefix>.save
  so it restarts from the converged PBE wavefunctions.
- Forces on the 8 I atoms are parsed from each SCF output.
- 24x24 force-constant block on I DOF: Phi_ij = -(F_i(+d_j) - F_i(-d_j)) / (2*delta).
- Mass-weight, diagonalize -> 24 frequencies (cm-1).
- Raman: bond-polarizability model (relative intensities) by DEFAULT. The script also has a
  hook to parse lepsilon polarizability tensors if ./eps_calcs were run (not done here; full
  epsilon adds ~2x cost per displacement and 48 extra SCFs per system).

Run as:
  python3 raman_driver.py gen      <sys>     # generate 48 displaced inputs + seed outdirs
  python3 raman_driver.py run      <sys>     # queue SCFs on 2 GPUs (restartable via .done)
  python3 raman_driver.py status   <sys>     # how many done / running
  python3 raman_driver.py hessian  <sys>     # parse forces, build 24x24 Phi, diagonalize
  python3 raman_driver.py raman    <sys>     # bond-polarizability Raman activities
  python3 raman_driver.py all      <sys>     # gen -> run -> hessian -> raman (single system)
"""
import os, sys, re, time, subprocess, shutil, signal
import numpy as np

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442"
QE = "/workspace/qe-7.5"
MPI = "/workspace/openmpi"
NVHPC = "/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
CUDA = "13.0"

# env for QE GPU
def qe_env():
    e = dict(os.environ)
    e["PATH"] = MPI + "/bin:" + e.get("PATH","")
    e["LD_LIBRARY_PATH"] = (MPI + "/lib:" + NVHPC + "/compilers/lib:" +
                            NVHPC + f"/math_libs/{CUDA}/lib64:" +
                            NVHPC + f"/cuda/{CUDA}/lib64")
    return e

# I atoms are the LAST 8 atoms (QE indices 73..80 in 1-based list of nat=80)
N_I = 8
I_FIRST = 73   # 1-based QE index of first I atom
DELTA_A = 0.01 # Angstrom (finite-displacement step)
M_I = 126.904  # amu (iodine mass)
# cm-1 conversion: nu = sqrt(eig [Ha/bohr^2 / amu]) * (0.5/pi/c)
# Standard QE factor (sqrt(Ha/bohr^2/amu) -> cm-1): 5140.4843
# Derivation: omega = sqrt(k/m). k in Ha/bohr^2, m in amu.
#   1 Ha/bohr^2 = 4.3597447222071e-18 J / (5.29177210903e-11 m)^2 = 1.556893e9 N/m
#   1 amu = 1.66053906660e-27 kg
# Frequency conversion. QE forces are Ry/bohr, so the central-difference Hessian
# Phi = dF/dx is in Ry/bohr^2; mass-weighted by amu -> eigenvalues in (Ry/bohr^2)/amu.
# nu(cm-1) = sqrt(eig) * factor.
#   1 Ry/bohr^2 = 778.409 N/m ; 1 amu = 1.660539e-27 kg
#   factor = sqrt(778.409/1.660539e-27) / (2*pi*2.99792458e10) = 3634.8 cm-1
# Validated on gas-phase I2 (same pipeline): I-I stretch = 214.5 cm-1 vs exp 213.
# NOTE: 5140.4843 is the HARTREE/bohr^2/amu factor; applying it to a Ry-valued Phi
# overestimates every frequency by sqrt(2). Do NOT use 5140 here.
RY_TO_CM = 3634.77
# Phi is in (Ry/bohr)/bohr = Ry/bohr^2 -> use RY_TO_CM directly.
# Bond vs vdw-contact split for the Raman bond-polarizability model. The pooled
# I-I distance spectrum has a clean 0.28 A coordination gap with nothing between
# 3.121 and 3.402 A, so <=3.25 A = bonded (primary + secondary polyiodide bond,
# carries bond polarizability) and >3.25 A = non-bonded vdw contact (e.g. st's
# 3.61 A branching-I ... in-chain-I), which carries NO bond polarizability and
# must be excluded -- otherwise it leaks spurious intensity into chain-bending
# modes (the st 103 cm-1 band was such an artifact: act 12 -> 2.5 once excluded).
BOND_CUTOFF = 3.25
# per-system FD smearing (env override; act_mol_rd runs at 0.03 after its SCF floor)
FD_DEGAUSS = float(__import__("os").environ.get("FD_DEGAUSS", "0.02"))

def sys_dir(sys): return f"{BASE}/{sys}"
def raman_dir(sys): return f"{BASE}/{sys}/raman"

# ---------- geometry parsing ----------
def parse_scf_loaded(sys):
    """Return (control_dict_text, system_text, electrons_text, species_text, cell(3x3 A),
       positions list of (el, frac[3]), kpoints_text). All as raw text blocks for re-emit."""
    f = f"{sys_dir(sys)}/scf_loaded.in"
    with open(f) as fh: txt = fh.read()
    # split into namelist chunks by lines
    lines = txt.split("\n")
    blocks = {}
    cur = None
    for ln in lines:
        s = ln.strip()
        m = re.match(r'^&(\w+)', s)
        if m:
            cur = m.group(1); blocks[cur] = []
            continue
        if cur and s == '/':
            blocks[cur].append(ln); cur = None; continue
        if cur is not None:
            blocks[cur].append(ln)
    # extract cell
    cell = []; in_cell=False
    for ln in lines:
        s=ln.strip()
        if s.startswith("CELL_PARAMETERS"):
            in_cell=True; continue
        if in_cell:
            if s=="" or s.startswith("ATOMIC"): in_cell=False
            else:
                parts=s.split()
                if len(parts)>=3:
                    try: cell.append([float(x) for x in parts[:3]])
                    except: pass
    cell = np.array(cell)
    # atomic positions
    pos=[]; in_pos=False
    for ln in lines:
        s=ln.strip()
        if s.startswith("ATOMIC_POSITIONS"):
            in_pos=True; continue
        if in_pos:
            if s.startswith("K_POINTS") or s.startswith("ATOMIC"): in_pos=False
            else:
                parts=s.split()
                if len(parts)>=4 and parts[0] in ("C","H","O","Y","Tb","Gd","I"):
                    pos.append((parts[0], np.array([float(x) for x in parts[1:4]])))
    assert len(pos)==80, f"expected 80 atoms, got {len(pos)}"
    assert pos[72][0]=="I" and pos[79][0]=="I", "atoms 73,80 not I!"
    # capture ATOMIC_SPECIES card (verbatim lines including header) -- it's a card, not a namelist
    species_lines=[]; in_sp=False
    for ln in lines:
        s=ln.strip()
        if s.startswith("ATOMIC_SPECIES"):
            in_sp=True; species_lines.append("ATOMIC_SPECIES"); continue
        if in_sp:
            if s.startswith("CELL_PARAMETERS") or s.startswith("ATOMIC_POSITIONS") or s.startswith("K_POINTS"):
                in_sp=False
            else:
                species_lines.append(ln)
    blocks["ATOMIC_SPECIES"] = species_lines
    return blocks, cell, pos

def prefix_for(sys): return f"{sys}_loaded"  # match scf_loaded.in prefix

# ---------- input generation ----------
def gen_displaced(sys):
    rdir = raman_dir(sys); os.makedirs(rdir, exist_ok=True)
    blocks, cell, pos = parse_scf_loaded(sys)
    I_pos_frac = np.array([pos[72+i][1] for i in range(N_I)])  # (8,3) fractional
    # cartesian for delta displacement (delta is in Angstrom)
    # delta_cart -> delta_frac via cell^-1
    invcell = np.linalg.inv(cell)
    # dump reference geometry for later analysis
    np.save(f"{rdir}/cell.npy", cell)
    np.save(f"{rdir}/I_frac_ref.npy", I_pos_frac)
    # save full positions for any post-analysis
    with open(f"{rdir}/pos_ref.txt","w") as fh:
        for el,fc in pos: fh.write(f"{el} {fc[0]:.12f} {fc[1]:.12f} {fc[2]:.12f}\n")
    # build the 48 displaced inputs. Each tag entry: (i_atom, axis, sign_str, tag, sign_num)
    tags=[]
    for i in range(N_I):       # I atom 0..7
        for a in range(3):      # x,y,z
            for s,sign in [(1,'p'),(-1,'n')]:
                tag = f"I{73+i:03d}_{'xyz'[a]}{sign}"
                tags.append((i,a,sign,tag,s))
    # assemble namelists from blocks (verbatim, but override prefix/outdir in CONTROL)
    def emit_block(name, fh, overrides=None):
        overrides = overrides or {}
        skip_keys = set(overrides.keys())
        fh.write(f"&{name}\n")
        for ln in blocks[name]:
            ss=ln.strip()
            if ss=='/' or ss=='': continue
            key = ss.split('=')[0].strip().lower()
            if key in skip_keys: continue
            fh.write(ln+"\n")
        for k,v in overrides.items():
            if isinstance(v,bool):
                fh.write(f"  {k} = .{('true' if v else 'false')}.\n")
            elif isinstance(v,str):
                fh.write(f"  {k} = '{v}'\n")
            else:
                fh.write(f"  {k} = {v}\n")
        fh.write("/\n")
    # the displaced position writer: modify pos[72+i] along cartesian axis 'a' by sign*delta
    for (i,a,sign,tag,snum) in tags:
        dcart = np.zeros(3); dcart[a] = snum*DELTA_A
        dfrac = dcart @ invcell
        newpos = [list(p) for p in pos]  # copy
        newpos[72+i] = (pos[72+i][0], pos[72+i][1] + dfrac)
        with open(f"{rdir}/disp_{tag}.in","w") as fh:
            emit_block("CONTROL", fh, {
                "prefix": prefix_for(sys),
                "outdir": f"./tmp_disp_{tag}",
                "calculation":"scf",
                "pseudo_dir":"/workspace/qe-7.5/pseudo",
                "tprnfor": True,
                "tstress": False,   # only forces needed
                # IMPORTANT: do NOT use restart_mode='restart'. For this metallic I2@MOF system,
                # restarting from the converged undisplaced wavefunction causes the SCF to settle
                # back into the UNDISPLACED basin (verified: gives identical forces for +/- disp,
                # zero Hessian). from_scratch (default) with mixing_beta=0.1 converges reliably
                # in ~42 iters and correctly resolves the 0.01 A displacement.
            })
            emit_block("SYSTEM", fh, {
                "degauss": FD_DEGAUSS,   # CRITICAL: degauss=0.01 (loaded value) causes Fermi-level bisection
                                    # DEADLOCK on displacements that collapse the Gamma-point gap below
                                    # the smearing width (verified: I074_zp, I075_yn, I075_zp hang at 0%
                                    # CPU in QE's Fermi search). degauss=0.03 cures it (converges in ~40
                                    # iters). For metallic polyiodide, 0.03 Ry smearing barely affects forces.
            })
            emit_block("ELECTRONS", fh, {
                "mixing_beta": 0.1,   # 0.2 (loaded) stalls at ~0.02 Ry; 0.1 converges cleanly to 1e-8.
                                       # (NB: mixing_beta is NOT the cure for the Fermi deadlock -- degauss is.)
            })
            # atomic species / cell / positions
            fh.write("\n".join(blocks.get("ATOMIC_SPECIES",[]))+"\n")
            fh.write("CELL_PARAMETERS angstrom\n")
            fh.write(f"  {cell[0,0]:.10f} {cell[0,1]:.10f} {cell[0,2]:.10f}\n")
            fh.write(f"  {cell[1,0]:.10f} {cell[1,1]:.10f} {cell[1,2]:.10f}\n")
            fh.write(f"  {cell[2,0]:.10f} {cell[2,1]:.10f} {cell[2,2]:.10f}\n")
            fh.write("ATOMIC_POSITIONS crystal\n")
            for el,fc in newpos:
                fh.write(f"{el:<2s} {fc[0]:.12f} {fc[1]:.12f} {fc[2]:.12f}\n")
            fh.write("K_POINTS automatic\n  4 4 2 0 0 0\n")   # 442 production mesh (442 campaign)
    # write the tag list
    with open(f"{rdir}/tags.txt","w") as fh:
        for t in tags: fh.write(f"{t[0]} {t[1]} {t[2]} {t[3]}\n")
    print(f"[{sys}] wrote {len(tags)} displaced inputs to {rdir}/")
    return tags

# ---------- outdir seeding ----------
def seed_outdirs(sys):
    """Copy tmp_loaded/<prefix>.save into each tmp_disp_<tag>/ AND patch the IONIC
    atomic_positions block in data-file-schema.xml to the DISPLACED geometry.

    CRITICAL GOTCHA (verified): with restart_mode='restart', QE reads atomic positions
    from the .save dir and DISCARDS the input ATOMIC_POSITIONS card. So we must overwrite
    the ionic geometry in data-file-schema.xml with the displaced positions BEFORE running,
    else all 48 SCFs run on the undisplaced geometry and the Hessian is identically zero.

    The schema has TWO <atomic_positions> blocks:
      (1) ionic positions (in <atomic_structure nat= alat=>)     -- PATCH THIS
      (2) atomic wavefunction positions (num_of_atomic_wfc=)     -- LEAVE ALONE
    Patching (2) corrupts the wavefunction basis (verified: 4 eV error). The saved
    wfc1.dat is the converged PBE wavefunction; QE projects it onto the displaced
    geometry's basis automatically. Converges in ~40 iters vs 150+ for from_scratch.
    Schema positions are cartesian in BOHR (atomic units)."""
    rdir = raman_dir(sys)
    pf = prefix_for(sys)
    src = f"{sys_dir(sys)}/tmp_loaded/{pf}.save"
    assert os.path.isdir(src), f"missing wavefunction save: {src}"
    bohr = 0.529177210563  # A -> bohr
    blocks, cell, pos = parse_scf_loaded(sys)
    invcell = np.linalg.inv(cell)
    tags = [l.split() for l in open(f"{rdir}/tags.txt") if l.strip()]  # i,a,sign,tag[,snum]
    n_new=0; n_patched=0
    for r in tags:
        i=int(r[0]); a=int(r[1]); sign_str=r[2]; tag=r[3]; snum=int(r[5]) if len(r)>5 else (1 if sign_str=='p' else -1)
        dst_dir = f"{sys_dir(sys)}/tmp_disp_{tag}"
        dst_save = f"{dst_dir}/{pf}.save"
        os.makedirs(dst_dir, exist_ok=True)
        if not os.path.isdir(dst_save):
            shutil.copytree(src, dst_save)
            n_new+=1
        # patch XML: compute displaced cartesian (bohr) for all 80 atoms (only I atom i moves)
        dcart = np.zeros(3); dcart[a] = snum*DELTA_A
        dfrac = dcart @ invcell
        # build per-atom displaced cartesian in bohr
        schema_pos = []
        for k,(el,fc) in enumerate(pos):
            disp_fc = fc.copy()
            if k == 72+i:  # the displaced I atom (0-indexed 72+i = QE 73+i)
                disp_fc = fc + dfrac
            cart_A = disp_fc @ cell
            cart_bohr = cart_A / bohr
            schema_pos.append((el, k+1, cart_bohr))
        # rewrite the IONIC atomic_positions block in data-file-schema.xml.
        # The QE schema has TWO <atomic_positions> sections:
        #   (1) <atomic_structure nat="80" alat=...>        <- ionic positions (QE reads on restart)
        #   (2) <atomic_structure nat="80" num_of_atomic_wfc=...>  <- atomic wfc starting positions
        # We patch ONLY (1). Patching (2) corrupts the wavefunction basis (verified: gives 4 eV
        # error). Leaving (2) at the undisplaced positions is fine -- the saved wfc1.dat is the
        # converged PBE wavefunction, projected onto the displaced geometry's basis by QE.
        xml_path = f"{dst_save}/data-file-schema.xml"
        with open(xml_path) as fh: xml = fh.read()
        # build replacement lines (one per atom, in input order)
        new_lines=[]
        for el,idx,cb in schema_pos:
            new_lines.append(f'        <atom name="{el}" index="{idx}">'
                             f'{cb[0]:.15E}  {cb[1]:.15E}  {cb[2]:.15E}</atom>')
        new_block = "\n".join(new_lines)
        # replace ONLY THE FIRST contiguous run of <atom ...>...</atom> lines (the ionic block)
        new_xml, nsubs = re.subn(
            r'(\s*<atom name="\w+" index="\d+">[^\n]*</atom>\n?)+',
            new_block + "\n", xml, count=1)
        if nsubs == 0:
            print(f"[warn] {tag}: no atom block replaced in XML!")
        else:
            n_patched+=1
            with open(xml_path,"w") as fh: fh.write(new_xml)
    print(f"[{sys}] seeded {n_new} fresh outdirs, patched XML in {n_patched}/{len(tags)} (of {len(tags)} total) from {src}")

# ---------- SCF running ----------
def is_done(sys, tag):
    return os.path.exists(f"{raman_dir(sys)}/disp_{tag}.done")
def is_conv(sys, tag):
    p = f"{raman_dir(sys)}/disp_{tag}.out"
    if not os.path.exists(p): return False
    return "convergence has been achieved" in open(p,errors="ignore").read()

def queue_scfs(sys, poll=20, max_runtime_s=None):
    max_runtime_s = max_runtime_s or int(__import__("os").environ.get("MAX_SCF_S", "7200"))
    """Run all not-yet-done displaced SCFs, one per GPU. Restartable.
    GPU set is controlled by env var RAMAN_GPUS (comma-separated, default '0,1').
    Set RAMAN_GPUS=1 to confine to GPU1 only (e.g. when GPU0 is shared with another job)."""
    rdir = raman_dir(sys)
    pf = prefix_for(sys)
    GPUS = [int(x) for x in os.environ.get("RAMAN_GPUS","0,1").split(",") if x.strip()]
    tags = [l.split()[3] for l in open(f"{rdir}/tags.txt") if l.strip()]
    todo = [t for t in tags if not (is_done(sys,t) or is_conv(sys,t))]
    # also clear any partial outputs that lack convergence (so we re-run cleanly)
    redo=[]
    for t in tags:
        if not is_done(sys,t) and is_conv(sys,t):
            open(f"{rdir}/disp_{t}.done","w").write("ok")  # mark already-converged
    todo = [t for t in tags if not is_done(sys,t)]
    print(f"[{sys}] {len(todo)} SCFs to run (of {len(tags)})", flush=True)
    if not todo: return
    env = qe_env()
    running = {}  # gpu -> (tag, popen, t0)
    log = open(f"{rdir}/scf_queue.log","a")
    def launch(gpu):
        t = todo.pop(0)
        d = sys_dir(sys)
        # handle __robust suffix: input is disp_<base>_robust.in, output goes to disp_<base>.out
        # (so the normal parser finds it). Mark with a .robust-running sentinel.
        if t.endswith("__robust"):
            base = t[:-len("__robust")]
            in_file = f"{rdir}/disp_{base}_robust.in"
            out_file = f"{rdir}/disp_{base}.out"
            # clear any stale non-converged output
            if os.path.exists(out_file): os.remove(out_file)
        else:
            base = t
            in_file = f"{rdir}/disp_{t}.in"
            out_file = f"{rdir}/disp_{t}.out"
        # from_scratch mode: ensure a clean outdir exists (no save seeding needed).
        outdir = f"{d}/tmp_disp_{base}"
        robust_outdir = f"{d}/tmp_robust_{base}"
        if t.endswith("__robust"):
            os.makedirs(robust_outdir, exist_ok=True)
        else:
            os.makedirs(outdir, exist_ok=True)
        in_rel = in_file.replace(f"{d}/","")  # make relative to d (where pw.x runs)
        out_rel = out_file.replace(f"{d}/","")
        cmd = (f"cd {d} && CUDA_VISIBLE_DEVICES={gpu} "
               f"mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x "
               f"-nk 1 -nd 1 -nb 1 -nt 1 -inp {in_rel} "
               f"> {out_rel} 2>&1")
        running[gpu] = (t, subprocess.Popen(cmd, shell=True, env=env, preexec_fn=os.setsid), time.time())
        msg = f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {sys}/{t}"
        print(msg, flush=True); log.write(msg+"\n"); log.flush()
    # verify GPUs are free before launching
    def gpu_free(gpu):
        try:
            # query a specific GPU's memory via index in the CSV (avoid --id flag quirks)
            out = subprocess.check_output(
                ["nvidia-smi","--query-gpu=index,memory.used",
                 "--format=csv,noheader,nounits"], timeout=20).decode().strip()
            for line in out.splitlines():
                parts=line.split(",")
                if len(parts)>=2 and int(parts[0].strip())==gpu:
                    return int(parts[1].strip()) < 500
            return True
        except Exception as e:
            print(f"[warn] nvidia-smi check failed: {e}", flush=True)
            return True  # be permissive
    while todo or running:
        for gpu in GPUS:
            if gpu not in running and todo:
                if not gpu_free(gpu):
                    continue  # try again next loop
                launch(gpu)
        time.sleep(poll)
        for gpu in list(running):
            t,p,t0 = running[gpu]
            rc = p.poll()
            if rc is not None:
                # map __robust sentinel back to base tag for done/conv checking
                base = t[:-len("__robust")] if t.endswith("__robust") else t
                ok = is_conv(sys, base)
                dt = time.time()-t0
                if ok:
                    open(f"{rdir}/disp_{base}.done","w").write("ok")
                    # disk hygiene: a completed displaced SCF's .save is pure scratch
                    # (~4 GB at 442); the Hessian parses forces from the .out text.
                    import shutil as _sh
                    for od in (f"{sys_dir(sys)}/tmp_disp_{base}", f"{sys_dir(sys)}/tmp_robust_{base}"):
                        _sh.rmtree(od, ignore_errors=True)
                else:
                    # SCF failed to converge (Fermi bisection deadlock or metallic stall).
                    # Auto-retry ONCE with larger degauss (0.05 Ry) which cures the Fermi deadlock.
                    robust_marker = f"{rdir}/disp_{base}.robust"
                    if not t.endswith("__robust") and not os.path.exists(robust_marker):
                        rb_in = f"{rdir}/disp_{base}_robust.in"
                        if not os.path.exists(rb_in):
                            txt = open(f"{rdir}/disp_{base}.in").read()
                            txt = re.sub(r"degauss\s*=\s*0\.\d+", "degauss = 0.05", txt)
                            txt = txt.replace(f"tmp_disp_{base}", f"tmp_robust_{base}")
                            open(rb_in,"w").write(txt)
                        open(robust_marker,"w").write("robust")
                        open(f"{rdir}/../MIXED_DEGAUSS_{sys}","w").write(
                            f"{time.strftime('%Y-%m-%d %H:%M')} {base} needed degauss=0.05 - "
                            "Hessian would MIX 0.02 and 0.05 displacements: REGENERATE the whole "
                            "set at degauss=0.03 before trusting frequencies!")
                        todo.append(base + "__robust")
                        msg = f"[{time.strftime('%H:%M:%S')}] GPU{gpu} FAIL {sys}/{base} rc={rc} ({dt:.0f}s) -> re-queued as robust (degauss=0.05)"
                    else:
                        msg = f"[{time.strftime('%H:%M:%S')}] GPU{gpu} FAIL {sys}/{base} rc={rc} ({dt:.0f}s) robust also failed -- manual review needed"
                    print(msg, flush=True); log.write(msg+"\n"); log.flush()
                del running[gpu]
            elif time.time()-t0 > max_runtime_s:
                # kill runaway
                try: os.killpg(p.pid, signal.SIGTERM)
                except: pass
                msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} TIMEOUT {sys}/{t} >{max_runtime_s}s, killed"
                print(msg,flush=True); log.write(msg+"\n"); log.flush()
                del running[gpu]
    # final summary
    ndone = sum(1 for t in tags if is_done(sys,t))
    print(f"[{sys}] queue complete: {ndone}/{len(tags)} done", flush=True)

# ---------- force parsing ----------
FORCE_RE = re.compile(r'^\s*atom\s+(\d+)\s+type.*?force\s+=\s+(-?\S+)\s+(-?\S+)\s+(-?\S+)', re.MULTILINE)
def parse_forces(sys, tag):
    """Return (8,3) array of forces (Ry/bohr) on the 8 I atoms from disp_<tag>.out."""
    p = f"{raman_dir(sys)}/disp_{tag}.out"
    txt = open(p, errors="ignore").read()
    # find all "Total force = ... " blocks; we want the LAST one (final SCF forces)
    # Actually each "ATOMIC_FORCES" block lists per-atom. QE writes:
    #     atom    1 type  1   force =    ...
    # We take the LAST such block.
    # split by "ATOMIC_FORCES" if present else by "Forces acting on atoms"
    blocks = re.split(r'(?:ATOMIC_FORCES|Forces acting on atoms)', txt)
    last = blocks[-1]
    # collect atom force lines within ~120 lines
    forces = {}
    for ln in last.splitlines()[:200]:
        m = re.match(r'\s*atom\s+(\d+)\s+type\s+\S+\s+force\s+=\s+(-?\S+)\s+(-?\S+)\s+(-?\S+)', ln)
        if m:
            idx = int(m.group(1)); f = [float(m.group(2)),float(m.group(3)),float(m.group(4))]
            forces[idx]=f
    # build (8,3) for I atoms 73..80
    F = np.zeros((N_I,3))
    for i in range(N_I):
        idx = 73+i
        if idx not in forces:
            raise RuntimeError(f"force on atom {idx} missing in {p}")
        F[i] = forces[idx]
    return F

# ---------- Hessian assembly ----------
def build_hessian(sys):
    rdir = raman_dir(sys)
    # tags.txt format: "i_atom axis sign_str tag"  (i, axis are ints)
    raw = [l.split() for l in open(f"{rdir}/tags.txt") if l.strip()]
    tags = [(int(r[0]), int(r[1]), r[2], r[3]) for r in raw]
    # Phi is 24x24 in Ry/bohr^2  (force Ry/bohr per bohr displacement)
    Phi = np.zeros((3*N_I, 3*N_I))
    bohr = 0.529177210563  # A -> bohr
    delta_bohr = DELTA_A/bohr
    # build dict (i,a,sign) -> parsed force row vector (24,)
    F = {}
    for i,a,sign,tag in tags:
        if not is_done(sys, tag):
            print(f"[warn] {tag} not done; skipping (Hessian will be incomplete)")
            continue
        F[(i,a,sign)] = parse_forces(sys, tag).flatten()  # (24,)
    # Phi[col] = -(F[+,col] - F[-,col]) / (2*delta) ; F[+,col] = force when atom col_i moved +delta along col_a
    for i,a,_,_ in tags:
        col = 3*i+a
        if (i,a,'p') in F and (i,a,'n') in F:
            dF = F[(i,a,'p')] - F[(i,a,'n')]
            Phi[:,col] = -dF / (2*delta_bohr)   # Ry/bohr^2
        else:
            print(f"[warn] column {col} (atom {73+i}, axis {a}) missing +/- pair")
    np.save(f"{rdir}/Phi.npy", Phi)
    # mass-weight
    masses = np.array([M_I]*3*N_I)
    D = Phi / np.sqrt(masses[:,None]*masses[None,:])
    evals, evecs = np.linalg.eigh(D)
    # frequencies: nu = sign-aware sqrt(|eig|) * factor, imaginary if eig<0
    freqs = np.sign(evals)*np.sqrt(np.abs(evals))*RY_TO_CM
    np.save(f"{rdir}/evals.npy", evals)
    np.save(f"{rdir}/evecs.npy", evecs)
    np.save(f"{rdir}/freqs.npy", freqs)
    # report
    order = np.argsort(freqs)[::-1]  # highest first
    print(f"\n=== {sys}: 24 I-sublattice phonon frequencies ===")
    print(f"{'mode':>4} {'freq(cm-1)':>12} {'eig(Ry/b^2/amu)':>18}")
    for k,ix in enumerate(order):
        flag = " (imag)" if freqs[ix] < 0 else ""
        print(f"{k+1:>4} {freqs[ix]:>12.2f} {evals[ix]:>18.6e}{flag}")
    print(f"\nPositive (real) modes: {(freqs>0).sum()}")
    print(f"Imaginary modes: {(freqs<0).sum()}")
    print(f"Near-zero (|f|<10 cm-1): {(np.abs(freqs)<10).sum()}")
    return freqs, evecs

# ---------- mode -> bond stretching analysis ----------
def mode_bond_character(sys, freqs, evecs, topK=12):
    """For each mode, project eigenvector onto I-I bond-stretch coordinates and report
    which bonds it stretches. Bond set = all I-I <= BOND_CUTOFF (min-image;
    vdw contacts excluded from the polarizability model)."""
    rdir = raman_dir(sys)
    cell = np.load(f"{rdir}/cell.npy")
    Ifrac = np.load(f"{rdir}/I_frac_ref.npy")  # (8,3)
    Icart = Ifrac @ cell
    # bond list: pairs with min-image dist <= BOND_CUTOFF (vdw contacts excluded)
    bonds=[]
    for i in range(N_I):
        for j in range(i+1,N_I):
            df = Ifrac[i]-Ifrac[j]; df-=np.round(df)
            r = np.linalg.norm(df@cell)
            if r <= BOND_CUTOFF:
                u = (df@cell)/r  # unit vector i->j
                bonds.append((i,j,r,u))
    print(f"\n=== {sys}: bond stretch character of top {topK} modes ===")
    print(f"Bond set ({len(bonds)} bonds <4A):")
    for i,j,r,u in sorted(bonds,key=lambda b:b[2]):
        print(f"  I{i+1}-I{j+1}: {r:.3f} A")
    # for each mode eigenvector e (24,), the displacement of atom i is e[3i:3i+3]
    # bond stretch amplitude = (u_i . u) where u_i is the relative displacement along bond
    # stretch coord = (disp_j - disp_i) . u / sqrt(2)  (mass-weighted already)
    order = np.argsort(freqs)[::-1]
    rows=[]
    for k,ix in enumerate(order[:topK]):
        e = evecs[:,ix]
        stretches=[]
        for (i,j,r,u) in bonds:
            di = e[3*i:3*i+3]; dj = e[3*j:3*j+3]
            s = np.dot(dj-di, u)/np.sqrt(2)
            stretches.append((abs(s),i,j,r,s))
        stretches.sort(reverse=True)
        top = stretches[:3]
        desc = ", ".join([f"I{i+1}-I{j+1}({r:.2f}A){'+' if s>=0 else '-'}{abs(s):.2f}"
                          for _,i,j,r,s in top])
        rows.append((k+1, freqs[ix], desc))
        print(f"  mode {k+1:>2}  {freqs[ix]:>7.2f} cm-1  : {desc}")
    return bonds, rows

# ---------- bond-polarizability Raman model ----------
def bond_pol_raman(sys, freqs, evecs, bonds):
    """Relative Raman activities from I-I bond-polarizability derivative model.
    Activity_n ~ sum_b [ (d alpha_b / d r_b)^2 * (dQ_n/d r_b)^2 ]
    Bond polarizability derivative scales with bond ORDER (shorter => stiffer => larger dadr):
    use empirical dadr_b = A * (r0/r_b)^p with r0=2.67 A (I2 gas), p=4.
    Stochastic ratio: alpha_parallel derivative dominates (~10x alpha_perp). Use scalar."""
    rdir = raman_dir(sys)
    r0 = 2.67; p = 4.0; A=1.0
    dadr = np.array([A*(r0/r)**p for (_,_,r,_) in bonds])
    activities = np.zeros(len(freqs))
    for ix in range(len(freqs)):
        e = evecs[:,ix]
        s = np.zeros(len(bonds))
        for kb,(i,j,r,u) in enumerate(bonds):
            di = e[3*i:3*i+3]; dj = e[3*j:3*j+3]
            s[kb] = np.dot(dj-di, u)/np.sqrt(2)
        # Raman activity ~ sum_b (dadr_b * s_b)^2  (parallel bond-polarizability derivative)
        activities[ix] = np.sum((dadr*s)**2)
    # normalize to max=100 for relative reporting
    if activities.max()>0:
        activities = 100*activities/activities.max()
    np.save(f"{rdir}/raman_activities.npy", activities)
    return activities

# ---------- dispatch ----------
def main():
    if len(sys.argv)<3:
        print(__doc__); sys.exit(1)
    cmd = sys.argv[1]; sysname = sys.argv[2]
    assert sysname in ("st_poly_lS5","vc_poly_lS3","act_poly_lS2","act_mol_rd","adduct"), "unknown sys (442 set)"
    if cmd=="gen":
        gen_displaced(sysname)
        # NOTE: seed_outdirs is no longer needed (from_scratch mode ignores saved wavefunctions).
    elif cmd=="run":
        queue_scfs(sysname)
    elif cmd=="status":
        rdir = raman_dir(sysname)
        tags=[l.split()[3] for l in open(f"{rdir}/tags.txt") if l.strip()]
        ndone=sum(1 for t in tags if is_done(sysname,t))
        print(f"{sysname}: {ndone}/{len(tags)} SCFs done")
        running = [t for t in tags if (not is_done(sysname,t)) and os.path.exists(f"{rdir}/disp_{t}.out") and not is_conv(sysname,t)]
        if running: print("in-progress (no .done yet):", running[:5], "..." if len(running)>5 else "")
    elif cmd=="hessian":
        f,e = build_hessian(sysname)
        mode_bond_character(sysname, f, e)
    elif cmd=="raman":
        rdir = raman_dir(sysname)
        f = np.load(f"{rdir}/freqs.npy"); e = np.load(f"{rdir}/evecs.npy")
        cell = np.load(f"{rdir}/cell.npy"); Ifrac = np.load(f"{rdir}/I_frac_ref.npy")
        bonds=[]
        for i in range(N_I):
            for j in range(i+1,N_I):
                df=Ifrac[i]-Ifrac[j]; df-=np.round(df)
                r=np.linalg.norm(df@cell)
                if r<=BOND_CUTOFF: bonds.append((i,j,r,(df@cell)/r))
        act = bond_pol_raman(sysname, f, e, bonds)
        order=np.argsort(f)[::-1]
        print(f"\n=== {sysname}: predicted Raman (bond-pol model) ===")
        print(f"{'mode':>4} {'freq':>8} {'activity':>10}")
        for k,ix in enumerate(order):
            print(f"{k+1:>4} {f[ix]:>8.2f} {act[ix]:>10.2f}")
    elif cmd=="all":
        gen_displaced(sysname)
        queue_scfs(sysname)
        f,e = build_hessian(sysname)
        mode_bond_character(sysname, f, e)
    else:
        print(__doc__); sys.exit(1)

if __name__=="__main__":
    main()
