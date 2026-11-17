#!/usr/bin/env python3
"""
build_gui_data.py -- collect every converged single-phase case into the data
files read by the CFD viewer (cfd/results/gui/).

Usage (WSL, from cfd/):  python3 build_gui_data.py
  * re-runs postprocess.py and viz_axisym.py on each case (so every case has
    the same metrics and a view3d.json)
  * writes results/gui/cases/<geom>__<case>.json (field data + residual and
    inlet-pressure histories) and results/gui/index.json (metrics per case)
VoF cases (current_vof) are skipped.
"""
import glob
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "gui")
RHO, G = 997.0, 9.81
DESC = {
    "current": "Current part", "A_1in": "Prototype A, 1in inlet", "A_3-4in": "Prototype A, 3/4in inlet",
    "B_1in": "Prototype B, 1in inlet", "B_3-4in": "Prototype B, 3/4in inlet",
    "V1": "V1 throat 14 / exit 26", "V2": "V2 throat 16 / exit 26", "V3": "V3 throat 18 / exit 26",
    "V4": "V4 throat 14 / exit 30, long",
}


def cases():
    out = []
    for g in DESC:
        for c in sorted(glob.glob(os.path.join(HERE, g, "Q*"))):
            if os.path.isfile(os.path.join(c, "results.json")):
                out.append(c)
    return out


def history(case, sub, fname, cols):
    rows = {}
    for f in sorted(glob.glob(os.path.join(case, "postProcessing", sub, "*", fname)),
                    key=lambda p: float(os.path.basename(os.path.dirname(p)))):
        for line in open(f):
            if line.startswith("#"):
                continue
            s = line.replace("(", " ").replace(")", " ").split()
            try:
                t = int(float(s[0]))
                rows[t] = [float(s[k]) if s[k] != "N/A" else None for k in cols]
            except (ValueError, IndexError):
                pass
    ts = sorted(rows)
    keep = [t for t in ts if t % 20 == 0 or t == ts[-1]] if ts else []
    return keep, [rows[t] for t in keep]


def one(case):
    rel = os.path.relpath(case, HERE)
    subprocess.run(["python3", "postprocess.py", rel], cwd=HERE, capture_output=True)
    subprocess.run(["python3", "viz_axisym.py", rel], cwd=HERE, capture_output=True)
    geom, cname = rel.split(os.sep)
    v = json.load(open(os.path.join(case, "viz", "view3d.json")))
    res = json.load(open(os.path.join(case, "results.json")))
    info = json.load(open(os.path.join(case, "mesh_info.json")))
    rdir = [d for d in os.listdir(os.path.join(case, "postProcessing")) if d.startswith("residuals")][0]
    t, r = history(case, rdir, "residuals.dat", [1, 2, 3, 5, 6])          # p Ux Uy k omega
    ti, pin = history(case, "inletAvg", "surfaceFieldValue.dat", [1])
    v["residuals"] = {"it": t, "p": [x[0] for x in r], "Ux": [x[1] for x in r], "Uy": [x[2] for x in r],
                      "k": [x[3] for x in r], "omega": [x[4] for x in r]}
    v["inlet_p_hist"] = {"it": ti, "level_cm": [None if x[0] is None else x[0] / G * 100 for x in pin]}
    key = f"{geom}__{cname}"
    json.dump(v, open(os.path.join(OUT, "cases", key + ".json"), "w"), separators=(",", ":"))
    st = info["stations"]
    drop = (st["outlet"] - st["inlet"]) / 10.0
    m = {k: val for k, val in res.items() if not k.startswith("_") and k not in ("residuals",)}
    m.update({"key": key, "geometry": geom, "case": cname, "desc": DESC[geom],
              "mesh_level": info["level"], "drop_cm": drop,
              "internal_loss_cm": res["level_above_inlet_cm"] + drop - res["outlet_KE_head_cm"],
              "inlet_D_mm": 2 * info["inlet_radius_mm"], "outlet_D_mm": 2 * info["outlet_radius_mm"],
              "throat_D_mm": 2 * min(info["R_nodes_mm"]), "length_mm": st["outlet"],
              "res_final": res.get("residuals", {})})
    return m


if __name__ == "__main__":
    os.makedirs(os.path.join(OUT, "cases"), exist_ok=True)
    with ThreadPoolExecutor(8) as ex:
        idx = list(ex.map(one, cases()))
    json.dump({"cases": idx}, open(os.path.join(OUT, "index.json"), "w"), indent=0)
    print(len(idx), "cases ->", OUT)
