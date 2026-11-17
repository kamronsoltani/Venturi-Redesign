"""sweep_status.py <assembly dir> "<flows>" <level>  -> Markdown progress page on stdout.

Used by run_sweep.sh (pushed to the GitHub branch "sweep-status" when PUSH_STATUS=1).
Reads results.json of finished cases and the live work dir of the running case.
"""
import glob, json, os, sys, time

asm, flows, level = sys.argv[1], sys.argv[2].split(), sys.argv[3]
name = os.path.basename(asm.rstrip("/"))
suffix = "" if level == "medium" else "_" + level
wroot = os.path.join(os.path.expanduser("~/turitap_cfd"), asm.split("/cfd/", 1)[1])


def tail(fn, n=1):
    try:
        with open(fn) as f:
            return [l for l in f.read().splitlines() if l.strip() and not l.startswith("#")][-n:]
    except OSError:
        return []


def end_iter(work):
    for l in open(os.path.join(work, "system", "controlDict")):
        if l.strip().startswith("endTime"):
            return l.split()[1].rstrip(";")
    return "?"


rows = []
for q in flows:
    case = os.path.join(asm, f"Q{q}{suffix}")
    work = os.path.join(wroot, f"Q{q}{suffix}")
    rj = os.path.join(case, "results.json")
    if os.path.exists(rj):
        r = json.load(open(rj))
        st = "done"
        flag = " (inlet static NEGATIVE: would draw air)" if r.get("static_inlet_cm", 1) < 0 else ""
        res = r.get("final_residuals", {}).get("p")
        rows.append(f"| {q} | {st} | {r.get('final_iteration', '?')} | {res:.1e} | {r['level_inlet_cm']:.1f} | "
                    f"{r['static_inlet_cm']:.1f}{flag} | {r['port_wall_p_kPa']['mean']:.2f} |" if res is not None else
                    f"| {q} | {st} | {r.get('final_iteration', '?')} | - | {r['level_inlet_cm']:.1f} | {r['static_inlet_cm']:.1f}{flag} | {r['port_wall_p_kPa']['mean']:.2f} |")
        continue
    stf = os.path.join(work, "STATUS")
    if os.path.exists(stf):
        st = open(stf).read().strip()
        it, pres, lev = "", "", ""
        rf = sorted(glob.glob(os.path.join(work, "postProcessing", "residuals*", "*", "residuals.dat")))
        if rf and "solving" in st:
            last = tail(rf[-1])
            if last:
                v = last[0].split()
                it, pres = f"{v[0]} / {end_iter(work)}", v[1]
            af = sorted(glob.glob(os.path.join(work, "postProcessing", "inletAvg", "*", "*.dat")))
            if af and tail(af[-1]):
                lev = f"{float(tail(af[-1])[0].split()[1]) / 9.81 * 100:.1f} (static, live)"
        rows.append(f"| {q} | {st} | {it} | {pres} | | {lev} | |")
    else:
        rows.append(f"| {q} | waiting | | | | | |")

print(f"# Sweep status: {name}{' (' + level + ' mesh)' if suffix else ''}\n")
print(f"Updated {time.strftime('%Y-%m-%d %H:%M:%S')} (machine local time)\n")
print("| Flow (L/min) | Status | Iteration | p residual | Level above inlet centre, total (cm) | Static at inlet (cm) | Port wall p (kPa) |")
print("|---|---|---|---|---|---|---|")
print("\n".join(rows))
log = os.path.join(asm, "sweep.log")
if os.path.exists(log):
    print("\n## Sweep log (last 15 lines)\n\n```")
    print("\n".join(open(log).read().splitlines()[-15:]))
    print("```")
