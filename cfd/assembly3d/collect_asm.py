"""collect_asm.py <repo root>

Gathers every finished 3-D assembly case (cfd/asm*/Q*/results.json) into
cfd/results/assembly_3d_stage5b.csv and lists every viewer file in
docs/3d/cases/ in docs/3d/index.json (sorted by assembly, then flow).
"""
import csv, glob, json, os, re, sys

repo = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
UPSTREAM = {
    "asm34": "3/4in E-Sig (flat side down) + 3/4in threaded elbow",
    "asm1e": "1in E-Sig (flat side down) + 1in elbow",
    "asm1lr": "1in long-radius reducer",
    "asm34lr": "3/4in E-Sig (flat side down) + 3/4in long-radius elbow (no reduction)",
}


def case_key(path):
    m = re.search(r"Q([0-9.]+)(?:_(\w+))?$", path)
    return float(m.group(1)), m.group(2) or "medium"


rows = []
for rj in sorted(glob.glob(os.path.join(repo, "cfd", "asm*", "Q*", "results.json"))):
    case = os.path.dirname(rj)
    asm = os.path.basename(os.path.dirname(case))
    q, level = case_key(case)
    r = json.load(open(rj))
    L, P = r["loss_cm"], r["port_wall_p_kPa"]
    ms = open(os.path.join(case, "mesh_summary.txt")).read().split() if os.path.exists(os.path.join(case, "mesh_summary.txt")) else []
    res = r.get("final_residuals") or {}
    rows.append({
        "assembly": asm, "upstream": UPSTREAM.get(asm.split("_")[0], ""), "venturi": asm.split("_", 1)[-1],
        "Q_Lpm": f"{q:g}", "mesh_level": level,
        "cells": ms[1] if len(ms) > 1 else "", "max_nonortho_deg": ms[4] if len(ms) > 4 else "",
        "iterations": r.get("final_iteration", ""),
        "level_inlet_centre_cm": round(r["level_inlet_cm"], 2), "level_spigot_entry_cm": round(r["level_spigot_cm"], 2),
        "level_venturi_inlet_cm": round(r["level_venturi_in_cm"], 2),
        "static_inlet_cm": round(r["static_inlet_cm"], 2), "static_venturi_inlet_cm": round(r["static_venturi_in_cm"], 2),
        "inlet_static_negative_flag": "YES" if r["static_inlet_cm"] < 0 else "",
        "loss_inlet_ext_spigot_cm": round(L.get("inlet_extension+spigot", float("nan")), 2),
        "loss_esig_cm": round(L.get("esig_reducer", float("nan")), 2),
        "loss_elbow_cm": round(L.get("elbow_incl_thread_cavity", float("nan")), 2),
        "loss_venturi_cm": round(L.get("venturi_total", float("nan")), 2),
        "outlet_ke_head_cm": round(r["outlet_ke_head_cm"], 2), "outlet_alpha": round(r["outlet_alpha"], 3),
        "port_p_mean_kPa": round(P["mean"], 3), "port_p_min_kPa": round(P["min"], 3), "port_p_max_kPa": round(P["max"], 3),
        "port_p_esig_side_kPa": round(P["at_plus_y"], 3), "port_p_opposite_kPa": round(P["at_minus_y"], 3),
        "diffuser_recovery_frac": round(r["diffuser_recovery_frac"], 3),
        "diffuser_reversed_area_frac": round(r["diffuser_reversed_area_frac"], 3),
        "yplus_venturi_mean": r.get("yplus", {}).get("wall_venturi", {}).get("mean", ""),
        "mass_imbalance_pct": round(r["mass_imbalance_pct"], 3),
        "final_res_p": f"{res['p']:.1e}" if res.get("p") is not None else "",
        "inlet_level_drift_last500_cm": round(r.get("inlet_p_last500_range_cm", float("nan")), 2),
    })
rows.sort(key=lambda d: (d["assembly"], d["mesh_level"] != "medium", float(d["Q_Lpm"])))
if rows:
    out = os.path.join(repo, "cfd", "results", "assembly_3d_stage5b.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} cases -> {out}")

cases = []
for fn in glob.glob(os.path.join(repo, "docs", "3d", "cases", "*.json")):
    m = re.match(r"(.+)_Q([0-9.]+)(?:_(\w+))?\.json$", os.path.basename(fn))
    if not m:
        continue
    with open(fn) as f:
        head = f.read(4000)
    lab = re.search(r'"label"\s*:\s*"([^"]*)"', head)
    cases.append((m.group(1), m.group(3) is not None, float(m.group(2)),
                  {"file": "cases/" + os.path.basename(fn), "label": lab.group(1) if lab else os.path.basename(fn)}))
cases.sort(key=lambda c: c[:3])
json.dump({"cases": [c[3] for c in cases]}, open(os.path.join(repo, "docs", "3d", "index.json"), "w"))
print(f"{len(cases)} viewer cases -> docs/3d/index.json")
