"""PTY driver for Multiwfn: generate a QTAIM gradient+contour plane map from a cube.

Multiwfn is a menu-driven Fortran program; its prompts don't flush to a plain pipe,
so we attach it to a pseudo-terminal and respond by matching the prompt text.
"""
import os, sys, pty, time, select, re

BIN = "/workspace/dimitry/calc/qe/i2_mof/tools/multiwfn/Multiwfn_2026.7.11_bin_Linux_noGUI/Multiwfn_noGUI"
OUTDIR = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out"


def drive(cube, atoms, out_prefix, save_post_cmds=("4", "1", "8", "3", "-8", "2", "0"),
          map_type="6", log=None):
    os.environ.setdefault("MultiwfnPATH",
        "/workspace/dimitry/calc/qe/i2_mof/tools/multiwfn/Multiwfn_2026.7.11_bin_Linux_noGUI")
    os.environ["OMP_STACKSIZE"] = "256m"
    os.chdir(OUTDIR)
    pid, fd = pty.fork()
    if pid == 0:  # child
        os.execv(BIN, [BIN])
    buf = ""
    logf = open(f"{OUTDIR}/{out_prefix}.ptylog", "w") if log else None

    def read_idle(t=2.0):
        nonlocal buf
        deadline = time.time() + t
        while time.time() < deadline:
            r, _, _ = select.select([fd], [], [], 0.3)
            if r:
                try:
                    chunk = os.read(fd, 4096).decode("utf-8", "ignore")
                except OSError:
                    break
                buf += chunk
                if logf:
                    logf.write(chunk); logf.flush()
                deadline = time.time() + 0.8  # extend if data flowing
            else:
                if buf and time.time() > deadline - 0.3:
                    break

    def send(s):
        os.write(fd, (s + "\n").encode())

    atoms_str = ",".join(str(a) for a in atoms)
    steps = [
        ("input file path", cube),          # load cube
        ("Press ENTER", ""),                 # non-orthogonal warning
        ("Main function menu", "4"),         # -> function 4 (plane map)
        ("real space function", "-1"),       # -1 = trilinear interp of loaded grid
        ("Select the type of map", map_type),  # 6 = gradient + contour
        ("Number of grids", ""),             # default grid resolution
        ("If using", "4"),                   # 4 = define plane by 3 atoms  (fallback match)
        ("three atoms", atoms_str),          # atom indices
    ]
    # Map idleness-based matching: read until idle, match a keyword, respond.
    stepi = 0
    last_sent = None
    # custom flow: we match the tail of buf against known prompt signatures in priority order
    rules = [
        ("input file path", lambda: send(cube)),
        ("Press ENTER button", lambda: send("")),
        ("e.g. 3 (Excluding", lambda: send("0")),       # how many files (if appears)
        ("Main function menu", lambda: send("4")),
        ("real space function", lambda: send("-1")),
        ("type of map", lambda: send(map_type)),
        ("Number of grids in", lambda: send("")),
        ("If using X, Y, Z", lambda: send("")),          # grid ext prompt variants
        ("method to define", lambda: send("4")),
        ("Define the plane", lambda: send("4")),
        ("indices of the three", lambda: send(atoms_str)),
        ("nuclear positions of the three", lambda: send(atoms_str)),
        ("three atoms", lambda: send(atoms_str)),
        # post-processing menu (appears after the map is generated):
        ("Return to upper-level menu", None),            # marker: we are in post-menu
    ]
    started = time.time()
    post_seen = False
    post_cmd_iter = iter(save_post_cmds)
    n_post = 0
    while True:
        if time.time() - started > 280:
            print("[timeout]"); break
        read_idle(2.5)
        # child exited?
        try:
            wpid, status = os.waitpid(pid, os.WNOHANG)
            if wpid != 0:
                break
        except ChildProcessError:
            break
        tail = buf[-1200:]
        matched = False
        for sig, action in rules:
            if sig.lower() in tail.lower():
                if action is None:
                    post_seen = True
                    # feed post-menu commands one per idle
                    nxt = next(post_cmd_iter, None)
                    if nxt is None:
                        send("q"); time.sleep(1); break
                    send(nxt); matched = True; time.sleep(0.5); break
                else:
                    action(); matched = True; time.sleep(0.5); break
        if not matched and post_seen:
            # in post-menu sub-prompts (e.g. color choice after enabling labels)
            nxt = next(post_cmd_iter, None)
            if nxt is None:
                send("q"); time.sleep(1)
            else:
                send(nxt); time.sleep(0.4)
        # clear matched tail to avoid re-matching the same prompt
        if matched:
            buf = buf[len(buf)//2:]
    try:
        os.close(fd)
    except OSError:
        pass
    try:
        os.waitpid(pid, 0)
    except ChildProcessError:
        pass
    if logf:
        logf.close()


if __name__ == "__main__":
    cube = sys.argv[1]
    atoms = [int(x) for x in sys.argv[2].split(",")]
    prefix = sys.argv[3]
    drive(cube, atoms, prefix, log=True)
    print("done", prefix)
