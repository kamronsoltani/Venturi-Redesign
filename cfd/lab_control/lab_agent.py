#!/usr/bin/env python3
"""Lab agent for the TuriTap 3-D CFD campaign (runs inside WSL Ubuntu-22.04).

Every 30 s it reads the state of every case (no writes to case folders), and
every 2 min it publishes that state as status.json to the GitHub branch
"lab-status" (a single commit, force-pushed; main is never touched). The page
docs/control/ on GitHub Pages reads it.

Commands come from GitHub Issues whose title starts with "[lab]" and that were
opened by the repository owner (checked against the issue author). Allowed:
  action: run       assembly: <name>  flows: <q q ...>  level: medium|coarse|fine
  action: campaign  flows: <q q ...>  assemblies: all | <name name ...>
  action: pull      (git pull --ff-only, only when nothing is running)
  action: stop      target: current | all
  action: cancel    job: <id>
  action: clear-queue
Runs are queued and started one at a time (one sweep at a time on this
machine), each in its own visible PowerShell window, exactly like the manual
commands in LAB_SETUP.md. The agent never changes physics, mesh or solver
settings, never commits, and never pushes to main.

Run:  wsl -d Ubuntu-22.04 -- python3 /mnt/c/cfd/lab_agent/lab_agent.py
"""
import csv, json, os, re, subprocess, sys, time, glob, urllib.request, datetime

OWNER, REPO_NAME = "kamronsoltani", "Venturi-Redesign"
REPO = "/mnt/c/cfd/Venturi-Redesign"
WIN_REPO = r"C:\cfd\Venturi-Redesign"
CFD = REPO + "/cfd"
HOME = os.path.expanduser("~")
TC = HOME + "/turitap_cfd"
STATE_DIR = "/mnt/c/cfd/lab_agent/state"
JOBS_DIR = "/mnt/c/cfd/jobs"
STATE_FILE = STATE_DIR + "/state.json"
STATUS_FILE = STATE_DIR + "/status.json"
BRANCH = "lab-status"
DISTRO = "Ubuntu-22.04"
LOOP_S, ISSUE_S, PUSH_S = 30, 90, 120
VERSION = "1.0"

ASSEMBLIES = ["asm34_current", "asm34_A", "asm34_Bvane",
              "asm1e_V1", "asm1e_V2", "asm1e_V3", "asm1e_V4",
              "asm1lr_V1", "asm1lr_V2", "asm1lr_V3", "asm1lr_V4"]
UP = {"asm34": "3/4in E-Sig + 3/4in threaded elbow",
      "asm1e": "1in E-Sig + McMaster 4880K464 1in elbow",
      "asm1lr": "1in long-radius reducing elbow"}
VEN = {"current": "current Venturi", "A": "Venturi A (14 mm throat, 19 mm exit)",
       "Bvane": "Venturi B (14 mm throat, 26 mm exit) + cross-vane",
       "V1": "V1 (14 mm throat, 26 mm exit) + vane", "V2": "V2 (16 mm throat, 26 mm exit) + vane",
       "V3": "V3 (18 mm throat, 26 mm exit) + vane", "V4": "V4 (14 mm throat, 30 mm exit) + 6-blade vane"}
LEVELS = ("medium", "coarse", "fine")
BUSY_PAT = "run_sweep.sh|run_campaign.sh|snappyHexMesh|foamRun"


def label(a):
    up, v = a.split("_", 1)
    return f"{UP[up]} + {VEN.get(v, v)}"


def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg):
    print(f"{datetime.datetime.now():%H:%M:%S} {msg}", flush=True)


def sh(cmd, inp=None, timeout=60):
    try:
        r = subprocess.run(cmd, input=inp, capture_output=True, timeout=timeout)
        return r.returncode, r.stdout.decode(errors="replace"), r.stderr.decode(errors="replace")
    except Exception as e:
        return 99, "", str(e)


def read(p, default=""):
    try:
        with open(p, errors="replace") as f:
            return f.read()
    except Exception:
        return default


def tail(p, n):
    return [l for l in read(p).splitlines() if l.strip()][-n:]


# ---------------------------------------------------------------- state
def load_state():
    try:
        s = json.load(open(STATE_FILE))
    except Exception:
        s = {}
    s.setdefault("since", now_iso())          # ignore issues older than the first start
    s.setdefault("processed", [])
    s.setdefault("queue", [])
    s.setdefault("next_id", 1)
    s.setdefault("cmdlog", [])
    s.setdefault("launched", None)
    return s


def save_state(s):
    os.makedirs(STATE_DIR, exist_ok=True)
    tmp = STATE_FILE + ".tmp"
    json.dump(s, open(tmp, "w"), indent=1)
    os.replace(tmp, STATE_FILE)


def cmdlog(s, issue, text, ok=True):
    s["cmdlog"].append({"t": now_iso(), "issue": issue, "ok": ok, "msg": text})
    s["cmdlog"] = s["cmdlog"][-40:]
    log(("OK  " if ok else "REJ ") + f"#{issue}: {text}")


# ---------------------------------------------------------------- status
def processes():
    rc, out, _ = sh(["pgrep", "-af", BUSY_PAT])
    lines = [l for l in out.splitlines() if "pgrep" not in l and "lab_agent" not in l]
    return lines


def csv_rows():
    rows = {}
    try:
        for r in csv.DictReader(open(CFD + "/results/assembly_3d_stage5b.csv")):
            try:
                rows[(r["assembly"], float(r["Q_Lpm"]), r.get("mesh_level", "medium"))] = r
            except Exception:
                pass
    except Exception:
        pass
    return rows


def fnum(x, nd=2):
    try:
        return round(float(x), nd)
    except Exception:
        return None


def case_entry(a, d, rows):
    m = re.match(r"^Q([\d.]+)(?:_(coarse|fine))?$", os.path.basename(d))
    if not m:
        return None
    q, lev = float(m.group(1)), m.group(2) or "medium"
    e = {"q": q, "level": lev, "dir": os.path.basename(d),
         "status": read(d + "/STATUS").strip() or None,
         "mesh": read(d + "/mesh_summary.txt").strip() or None}
    rj = None
    try:
        rj = json.load(open(d + "/results.json"))
    except Exception:
        pass
    if rj:
        e["result"] = {
            "level_inlet_cm": fnum(rj.get("level_inlet_cm")),
            "level_venturi_in_cm": fnum(rj.get("level_venturi_in_cm")),
            "static_inlet_cm": fnum(rj.get("static_inlet_cm")),
            "static_inlet_negative": bool(rj.get("static_inlet_negative_FLAG")),
            "port_p_mean_kPa": fnum((rj.get("port_wall_p_kPa") or {}).get("mean"), 3),
            "loss_venturi_cm": fnum((rj.get("loss_cm") or {}).get("venturi_total")),
            "mass_imbalance_pct": fnum(rj.get("mass_imbalance_pct"), 3),
            "flow_in_Lpm": fnum(rj.get("flow_in_Lpm")),
        }
    r = rows.get((a, q, lev))
    if r:
        e["drift_last500_cm"] = fnum(r.get("inlet_level_drift_last500_cm"), 3)
        e["final_res_p"] = r.get("final_res_p")
        e["iterations"] = r.get("iterations")
    # convergence verdict per LAB_SETUP B5
    if e.get("status", "") and e["status"].startswith("done") and rj:
        bad = []
        dr = e.get("drift_last500_cm")
        if dr is not None and abs(dr) >= 0.5:
            bad.append(f"drift {dr} cm")
        mi = e["result"]["mass_imbalance_pct"]
        if mi is not None and abs(mi) >= 0.5:
            bad.append(f"mass imbalance {mi}%")
        e["converged"] = not bad
        e["concerns"] = bad
    return e


_prev_it = {}


def current_case():
    w = read(TC + "/CURRENT_CASE").strip()
    if not w:
        return None
    c = {"work": w, "assembly": os.path.basename(os.path.dirname(w)), "case": os.path.basename(w),
         "stage": read(w + "/STATUS").strip() or None}
    try:
        c["stage_since_min"] = round((time.time() - os.path.getmtime(w + "/STATUS")) / 60, 1)
    except Exception:
        pass
    m = re.search(r"^endTime\s+(\d+)", read(w + "/system/controlDict"), re.M)
    c["end"] = int(m.group(1)) if m else None
    res = sorted(glob.glob(w + "/postProcessing/residuals*/*/residuals.dat"))
    if res:
        t = tail(res[-1], 1)
        if t and not t[0].startswith("#"):
            f = t[0].split()
            try:
                c["iteration"] = int(float(f[0]))
                c["res_p"] = float(f[1])
            except Exception:
                pass
        try:
            c["output_age_min"] = round((time.time() - os.path.getmtime(res[-1])) / 60, 1)
        except Exception:
            pass
        it = c.get("iteration")
        p = _prev_it.get(w)
        if it is not None:
            if p and it > p[0]:
                rate = (it - p[0]) / (time.time() - p[1]) * 60
                c["it_per_min"] = round(rate, 1)
                if c["end"]:
                    c["eta_min"] = round((c["end"] - it) / rate, 0)
            if not p or it != p[0]:
                _prev_it[w] = (it, time.time())
    fi = sorted(glob.glob(w + "/postProcessing/inletFlow/*/*.dat"))
    fo = sorted(glob.glob(w + "/postProcessing/outletFlow/*/*.dat"))
    try:
        c["flow_in_Lpm"] = round(-float(tail(fi[-1], 1)[0].split()[1]) * 120000, 2)
        c["flow_out_Lpm"] = round(float(tail(fo[-1], 1)[0].split()[1]) * 120000, 2)
    except Exception:
        pass
    snap = w + "/log.snappyHexMesh"
    if c["stage"] and "snappy" in c["stage"] and os.path.exists(snap):
        txt = read(snap)
        c["shell_refinement_iters"] = len(re.findall(r"Shell refinement iteration", txt))
        last = [l for l in txt.splitlines() if re.match(r"^(Surface refinement|Shell refinement|Feature refinement|Snapping|Adding layers|Layer addition|Finished|Morphing|Smoothing)", l)]
        c["snappy_phase"] = last[-1] if last else None
        warn = c["assembly"].startswith("asm1lr") and (c.get("stage_since_min", 0) > 35 or c["shell_refinement_iters"] > 20)
        c["mesh_stall_warning"] = bool(warn)
    return c


def build_status(s, procs):
    rows = csv_rows()
    grid = {}
    for a in ASSEMBLIES:
        cases = []
        for d in sorted(glob.glob(f"{CFD}/{a}/Q*")):
            e = case_entry(a, d, rows)
            if e:
                cases.append(e)
        grid[a] = {"label": label(a), "cases": cases,
                   "sweep_log": tail(f"{CFD}/{a}/sweep.log", 4)}
    running = bool(procs)
    cur = current_case() if running else None
    if cur:
        a = cur["assembly"]
        if a in grid:
            for e in grid[a]["cases"]:
                if e["dir"] == cur["case"]:
                    e["running"] = True
    sweep = read(TC + "/SWEEP").strip().split("|") if read(TC + "/SWEEP").strip() else None
    camp = read(TC + "/CAMPAIGN").strip() if running else ""
    meshcheck = []
    for l in read(CFD + "/assembly_previews/meshcheck_medium.txt").splitlines():
        p = [x.strip() for x in l.split("|")]
        if len(p) >= 3:
            meshcheck.append({"assembly": p[0], "summary": p[2]})
    rc, head, _ = sh(["git.exe", "-C", WIN_REPO, "rev-parse", "--short", "HEAD"], timeout=30)
    return {
        "generated": now_iso(), "agent_version": VERSION, "host": "lab PC",
        "running": running, "processes": [l[:140] for l in procs][:8],
        "current": cur,
        "sweep": {"assembly": os.path.basename(sweep[0]), "flows": sweep[1], "level": sweep[2]} if sweep and len(sweep) >= 3 else None,
        "campaign": camp or None,
        "campaign_log": tail(CFD + "/results/campaign_stage5b.log", 8),
        "grid": grid, "meshcheck": meshcheck,
        "queue": s["queue"], "cmdlog": s["cmdlog"][-15:][::-1],
        "repo_head": head.strip() or None,
        "last_push_ok": s.get("last_push_ok"),
    }


def push_status(st, s):
    data = json.dumps(st, indent=1).encode()
    os.makedirs(STATE_DIR, exist_ok=True)
    open(STATUS_FILE, "wb").write(data)
    g = ["git.exe", "-C", WIN_REPO]
    rc, blob, err = sh(g + ["hash-object", "-w", "--stdin"], inp=data)
    if rc:
        return False, "hash-object: " + err.strip()[:200]
    blob = blob.strip()
    rc, tree, err = sh(g + ["mktree"], inp=f"100644 blob {blob}\tstatus.json\n".encode())
    if rc:
        return False, "mktree: " + err.strip()[:200]
    rc, commit, err = sh(g + ["commit-tree", tree.strip(), "-m", "lab status " + st["generated"]])
    if rc:
        return False, "commit-tree: " + err.strip()[:200]
    rc, _, err = sh(g + ["push", "-q", "-f", "origin", f"{commit.strip()}:refs/heads/{BRANCH}"], timeout=120)
    if rc:
        return False, "push: " + err.strip()[:200]
    return True, "ok"


# ---------------------------------------------------------------- commands
FLOW_RE = re.compile(r"^\d{1,2}(\.\d{1,2})?$")


def parse_flows(txt):
    toks = txt.replace(",", " ").split()
    if not toks or len(toks) > 12:
        raise ValueError("flows: give 1-12 values")
    for t in toks:
        if not FLOW_RE.match(t) or not (1 <= float(t) <= 80):
            raise ValueError(f"bad flow '{t}' (1-80 L/min)")
    return " ".join(toks)


def parse_body(body):
    kv = {}
    for line in (body or "").splitlines():
        m = re.match(r"^\s*([a-z\-]+)\s*:\s*(.*?)\s*$", line.strip("`* "))
        if m:
            kv[m.group(1)] = m.group(2)
    return kv


def fetch_issues():
    url = (f"https://api.github.com/repos/{OWNER}/{REPO_NAME}/issues"
           f"?state=all&creator={OWNER}&sort=created&direction=desc&per_page=30")
    req = urllib.request.Request(url, headers={"User-Agent": "lab-agent", "Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def handle_issue(s, iss, procs):
    n = iss["number"]
    kv = parse_body(iss.get("body"))
    act = kv.get("action", "").lower()
    if act in ("run", "campaign") and kv.get("pull", "").lower() in ("yes", "true", "1"):
        job = {"id": s["next_id"], "type": "pull", "issue": n, "queued": now_iso()}
        s["next_id"] += 1; s["queue"].append(job)
        cmdlog(s, n, f"queued job {job['id']}: git pull before the run")
    try:
        if act == "run":
            a = kv.get("assembly", "")
            if a not in ASSEMBLIES:
                raise ValueError(f"unknown assembly '{a}'")
            flows = parse_flows(kv.get("flows", ""))
            lev = kv.get("level", "medium") or "medium"
            if lev not in LEVELS:
                raise ValueError(f"bad level '{lev}'")
            job = {"id": s["next_id"], "type": "sweep", "assembly": a, "flows": flows, "level": lev, "issue": n, "queued": now_iso()}
            s["next_id"] += 1; s["queue"].append(job)
            cmdlog(s, n, f"queued job {job['id']}: {a} at {flows} L/min ({lev})")
        elif act == "campaign":
            flows = parse_flows(kv.get("flows", "30 6 8") or "30 6 8")
            al = kv.get("assemblies", "all").strip() or "all"
            if al != "all":
                names = al.replace(",", " ").split()
                for x in names:
                    if x not in ASSEMBLIES:
                        raise ValueError(f"unknown assembly '{x}'")
                al = " ".join(names)
            job = {"id": s["next_id"], "type": "campaign", "flows": flows, "assemblies": al, "issue": n, "queued": now_iso()}
            s["next_id"] += 1; s["queue"].append(job)
            cmdlog(s, n, f"queued job {job['id']}: campaign flows [{flows}] assemblies [{al}]")
        elif act == "pull":
            job = {"id": s["next_id"], "type": "pull", "issue": n, "queued": now_iso()}
            s["next_id"] += 1; s["queue"].append(job)
            cmdlog(s, n, f"queued job {job['id']}: git pull (runs when nothing is running)")
        elif act == "stop":
            tgt = kv.get("target", "current").lower()
            if tgt == "current":
                rc, _, _ = sh(["pkill", "-f", "mpirun -np"])
                cmdlog(s, n, "stopped the current case's mesher/solver (the sweep logs it as FAILED and moves on)" if rc == 0 else "stop current: nothing was meshing or solving")
            elif tgt == "all":
                for pat in ("run_campaign.sh", "run_sweep.sh", "mpirun -np", "snappyHexMesh", "foamRun"):
                    sh(["pkill", "-f", pat])
                s["launched"] = None
                cmdlog(s, n, "stopped everything (campaign, sweep, mesher, solver); queue kept")
            else:
                raise ValueError("target must be current or all")
        elif act == "cancel":
            jid = int(kv.get("job", "0"))
            before = len(s["queue"])
            s["queue"] = [j for j in s["queue"] if j["id"] != jid]
            if len(s["queue"]) == before:
                raise ValueError(f"no queued job {jid}")
            cmdlog(s, n, f"cancelled queued job {jid}")
        elif act == "clear-queue":
            s["queue"] = []
            cmdlog(s, n, "queue cleared")
        else:
            raise ValueError(f"unknown action '{act}'")
    except ValueError as e:
        cmdlog(s, n, str(e), ok=False)


def poll_issues(s, procs):
    try:
        issues = fetch_issues()
    except Exception as e:
        log(f"issue poll failed: {e}")
        return False
    since = s["since"]
    new = []
    for iss in issues:
        if "pull_request" in iss or iss["number"] in s["processed"]:
            continue
        if not iss.get("title", "").startswith("[lab]"):
            continue
        if (iss.get("user") or {}).get("login", "").lower() != OWNER.lower():
            continue
        if iss.get("author_association") not in ("OWNER", "COLLABORATOR", "MEMBER"):
            continue
        if datetime.datetime.fromisoformat(iss["created_at"].replace("Z", "+00:00")) < datetime.datetime.fromisoformat(since):
            s["processed"].append(iss["number"]); continue
        new.append(iss)
    for iss in sorted(new, key=lambda x: x["number"]):
        handle_issue(s, iss, procs)
        s["processed"].append(iss["number"])
    s["processed"] = s["processed"][-300:]
    return bool(new)


# ---------------------------------------------------------------- launching
def launch_window(jid, title, script):
    os.makedirs(JOBS_DIR, exist_ok=True)
    path = f"{JOBS_DIR}/job_{jid}.sh"
    open(path, "w").write(script)
    inner = (f"$Host.UI.RawUI.WindowTitle = ''{title}''; "
             f"wsl.exe -d {DISTRO} -- bash -l {path}")
    ps = f"Start-Process powershell.exe -ArgumentList '-NoExit','-Command','{inner}'"
    rc, out, err = sh(["powershell.exe", "-NoProfile", "-Command", ps], timeout=60)
    return rc == 0, err.strip()[:200]


def start_next(s, procs):
    if procs or not s["queue"]:
        return False
    if s.get("launched") and time.time() - s["launched"] < 90:
        return False                        # give the last launch time to show up in pgrep
    job = s["queue"].pop(0)
    jid = job["id"]
    if job["type"] == "pull":
        rc, out, err = sh(["git.exe", "-C", WIN_REPO, "pull", "--ff-only"], timeout=300)
        cmdlog(s, job["issue"], f"job {jid} git pull: " + ((out or err).strip().splitlines() or ["?"])[-1][:160], ok=(rc == 0))
        return True
    hdr = "cd /mnt/c/cfd/Venturi-Redesign/cfd/assembly3d\nexport PUSH_STATUS=1\n"
    if job["type"] == "sweep":
        lab = label(job["assembly"]).replace("'", "")
        script = hdr + f"export LABEL='{lab}'\nbash run_sweep.sh ../{job['assembly']} '{job['flows']}' {job['level']}\n"
        title = f"CFD JOB {jid} - {job['assembly']} {job['flows']} Lpm (do not close)"
    else:
        asm = "" if job["assemblies"] == "all" else f"export ASSEMBLIES='{job['assemblies']}'\n"
        script = hdr + asm + f"export FLOWS_ORDER='{job['flows']}'\nbash run_campaign.sh\n"
        title = f"CFD JOB {jid} - campaign {job['flows']} (do not close)"
    script += 'echo "=== job finished (exit $?) at $(date) - this window can be closed ==="\n'
    ok, err = launch_window(jid, title, script)
    s["launched"] = time.time() if ok else None
    cmdlog(s, job["issue"], f"started job {jid} in a new window" if ok else f"job {jid} failed to start: {err}", ok=ok)
    if not ok:
        s["queue"].insert(0, job)
    return True


# ---------------------------------------------------------------- main
def main():
    os.makedirs(STATE_DIR, exist_ok=True)
    s = load_state(); save_state(s)
    log(f"lab agent {VERSION} started; commands accepted from issues opened after {s['since']}")
    last_issue = last_push = 0
    while True:
        try:
            procs = processes()
            changed = False
            if time.time() - last_issue >= ISSUE_S:
                changed |= poll_issues(s, procs); last_issue = time.time()
            changed |= start_next(s, procs)
            if changed:
                save_state(s); procs = processes()
            st = build_status(s, procs)
            if changed or time.time() - last_push >= PUSH_S:
                ok, msg = push_status(st, s)
                s["last_push_ok"] = now_iso() if ok else s.get("last_push_ok")
                if not ok:
                    log("status push FAILED: " + msg)
                last_push = time.time(); save_state(s)
            cur = st.get("current") or {}
            log(f"{'RUNNING' if st['running'] else 'idle   '} {cur.get('assembly','')}/{cur.get('case','')} "
                f"{cur.get('stage','') or ''} it {cur.get('iteration','-')}/{cur.get('end','-')} | queue {len(s['queue'])} | last push {s.get('last_push_ok') or 'never'}")
        except Exception as e:
            log(f"loop error: {e!r}")
        time.sleep(LOOP_S)


if __name__ == "__main__":
    main()
