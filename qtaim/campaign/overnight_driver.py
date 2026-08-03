#!/usr/bin/env python3
"""Overnight driver: poll empty_vcr -> compute binding corrections -> launch HSE."""
import subprocess, re, os, time, json

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
QE = "/workspace/qe-7.5/bin/pw.x"
MPI = "/workspace/openmpi"
NVHPC = "/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
CUDA = "13.0"
I2_REF = -3031.19956
RY2KJ = 1312.75

def env_setup():
    env = os.environ.copy()
    env["PATH"] = f"{MPI}/bin:" + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/{CUDA}/lib64:{NVHPC}/cuda/{CUDA}/lib64"
    return env

def Etot(path):
    try:
        t = open(path, errors="ignore").read()
        m = re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)", t)
        return float(m[-1]) if m else None
    except:
        return None

def is_done(path):
    try:
        return "JOB DONE" in open(path, errors="ignore").read()
    except:
        return False

def gpu_mem():
    r = subprocess.run("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits", shell=True, capture_output=True, text=True)
    return [int(x.strip()) for x in r.stdout.strip().split("\n")]

def launch_scf(sys, inp, gpu):
    d = f"{BASE}/{sys}"
    env = env_setup()
    subprocess.Popen(f"cd {d} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE} -nk 1 -nd 1 -nb 1 -nt 1 -inp {inp} > {inp.replace('.in','.out')} 2>&1", shell=True, env=env)
    print(f"[{time.strftime('%H:%M:%S')}] launched {sys}/{inp} on GPU{gpu}")

def wait_gpu(timeout=3600):
    """Wait for a GPU to free."""
    start = time.time()
    while time.time() - start < timeout:
        mems = gpu_mem()
        for g, m in enumerate(mems):
            if m < 500:
                return g
        time.sleep(30)
    return None

def main():
    log = open(f"{BASE}/overnight_log.txt", "w")
    def L(msg):
        print(msg); log.write(msg + "\n"); log.flush()

    SYSTEMS = ["Gd_act", "Y_actc", "Tb_act"]
    L(f"[{time.strftime('%H:%M:%S')}] overnight driver started")

    # Phase 1: wait for all empty_vcr RELAX + fix03 SCF to finish
    L("Phase 1: waiting for empty_vcr + fix03 calcs...")
    while True:
        all_done = True
        status = []
        for sys in SYSTEMS:
            for calc in ["scf_empty_vcr", "scf_empty_fix03"]:
                path = f"{BASE}/{sys}/{calc}.out"
                done = is_done(path) or (Etot(path) is not None and "fix03" in calc)
                if calc == "scf_empty_vcr" and not done:
                    done = is_done(path)
                status.append(f"{sys}/{calc}: {'done' if done else 'running'}")
                if not done:
                    all_done = False
        L(f"[{time.strftime('%H:%M:%S')}] " + " | ".join(status))
        if all_done:
            break
        time.sleep(120)

    # Phase 2: compute binding corrections
    L("\nPhase 2: computing binding corrections...")
    results = []
    for sys in SYSTEMS:
        loaded_fixed = Etot(f"{BASE}/{sys}/vc-relax.out")  # will get LAST = vc-relaxed
        # actually need FIRST (fixed) and LAST (vc) from vc-relax
        vc_lines = open(f"{BASE}/{sys}/vc-relax.out", errors="ignore").read().splitlines()
        vc_energies = [float(l.split("=")[1].split()[0]) for l in vc_lines if re.match(r"\s*!\s*total energy", l)]
        loaded_fixed = vc_energies[0]  # first = fixed-cell start
        loaded_vc = vc_energies[-1]    # last = vc-relaxed

        empty_fix03 = Etot(f"{BASE}/{sys}/scf_empty_fix03.out")
        empty_vcr = Etot(f"{BASE}/{sys}/scf_empty_vcr.out")

        binding_fixed = (loaded_fixed - empty_fix03 - I2_REF) * RY2KJ / 4
        binding_vc = (loaded_vc - empty_vcr - I2_REF) * RY2KJ / 4
        correction = binding_vc - binding_fixed

        L(f"\n{sys}:")
        L(f"  loaded_fixed={loaded_fixed:.5f} loaded_vc={loaded_vc:.5f}")
        L(f"  empty_fixed={empty_fix03:.5f} empty_vc={empty_vcr:.5f}")
        L(f"  binding_fixed(0.03) = {binding_fixed:.1f} kJ/I2")
        L(f"  binding_vc(0.03)     = {binding_vc:.1f} kJ/I2")
        L(f"  CORRECTION           = {correction:+.1f} kJ/I2")
        results.append({"sys": sys, "binding_fixed": binding_fixed, "binding_vc": binding_vc, "correction": correction})

    # write results
    with open(f"{BASE}/vc_relax_binding_corrections.json", "w") as f:
        json.dump(results, f, indent=2)
    L(f"\nResults written to vc_relax_binding_corrections.json")

    # Phase 3: launch HSE (Tier 1: degauss 0.03)
    L("\nPhase 3: launching HSE (degauss 0.03)...")
    gpu = wait_gpu()
    if gpu is None:
        L("No free GPU for HSE. Exiting.")
        return

    # build HSE input from Y_str/scf_loaded.in
    s = open(f"{BASE}/Y_str/scf_loaded.in").read()
    cards = s[s.index("ATOMIC_SPECIES"):]
    hse_input = """&CONTROL
  calculation = 'scf'
  prefix = 'Y_hse_dg03'
  outdir = './tmp_hse'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  verbosity = 'high'
/
&SYSTEM
  ibrav = 0
  nat = 80
  ntyp = 5
  ecutwfc = 50
  ecutrho = 400
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.03
  nspin = 1
  input_dft = 'hse'
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .true.
  exx_fraction = 0.25
  screening_parameter = 0.106
/
&ELECTRONS
  conv_thr = 1.0d-6
  diagonalization = 'david'
  diago_david_ndim = 2
  mixing_mode = 'local-TF'
  mixing_beta = 0.1
  mixing_ndim = 20
  electron_maxstep = 600
/
""" + cards + "K_POINTS gamma\n"
    open(f"{BASE}/Y_str/hse_dg03.in", "w").write(hse_input)
    L(f"HSE input written to Y_str/hse_dg03.in (degauss=0.03, ecut=50/400)")

    d = f"{BASE}/Y_str"
    env = env_setup()
    subprocess.Popen(f"cd {d} && rm -rf tmp_hse && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE} -nk 1 -nd 1 -nb 1 -nt 1 -inp hse_dg03.in > hse_dg03.out 2>&1", shell=True, env=env)
    L(f"HSE launched on GPU{gpu}")

    # Phase 4: monitor HSE, launch projwfc when done
    L("Phase 4: monitoring HSE...")
    hse_timeout = 7200  # 2 hours
    start = time.time()
    while time.time() - start < hse_timeout:
        hse_out = f"{BASE}/Y_str/hse_dg03.out"
        if is_done(hse_out):
            L(f"HSE DONE!")
            e = Etot(hse_out)
            L(f"HSE energy = {e}")
            # check if .save exists
            save_dir = f"{BASE}/Y_str/tmp_hse/Y_hse_dg03.save"
            dat_files = glob.glob(f"{save_dir}/*.dat") if os.path.exists(save_dir) else []
            L(f".save dat files: {len(dat_files)}")
            if dat_files:
                # launch projwfc
                L("Launching projwfc...")
                proj_input = f"""&PROJWFC
  prefix = 'Y_hse_dg03'
  outdir = './tmp_hse'
  filpdos = 'hse_pdos'
  Emin = -20
  Emax = 10
  DeltaE = 0.01
  degauss = 0.03
  ngauss = 0
/
"""
                open(f"{d}/hse_proj.in", "w").write(proj_input)
                subprocess.run(f"cd {d} && {QE.replace('pw.x','projwfc.x')} -inp hse_proj.in > hse_proj.out 2>&1", shell=True, env=env)
                L("projwfc done!")
            break
        # check for hang (no output change for 30 min)
        time.sleep(120)

    L(f"\n[{time.strftime('%H:%M:%S')}] overnight driver complete.")
    log.close()

if __name__ == "__main__":
    import glob
    main()
