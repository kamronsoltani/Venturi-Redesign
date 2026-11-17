#!/usr/bin/env python3
"""
postprocess_vof.py -- Stage 4 (VoF) summary of one case from its
postProcessing time histories. Writes results_vof.json in the case folder.

Averages are taken over the last AVG_WINDOW seconds of the run; the previous
window of the same length is used to judge whether the run is statistically
steady.

Usage: python3 postprocess_vof.py "<case dir>" [--window 0.5]
"""
import argparse
import glob
import json
import os
import re

RHO = 997.0
G = 9.81


def read_series(case, fo, fname):
    """Concatenate all time folders of postProcessing/<fo>/*/<fname>.
    Returns (times, rows, header comments)."""
    files = sorted(glob.glob(os.path.join(case, "postProcessing", fo, "*", fname)),
                   key=lambda p: float(os.path.basename(os.path.dirname(p))))
    data, comments = {}, []
    for f in files:
        for line in open(f):
            if line.startswith("#"):
                comments.append(line)
                continue
            line = re.sub(r"[()]", " ", line).split()
            if not line:
                continue
            data[float(line[0])] = [float(x) for x in line[1:]]
    t = sorted(data)
    return t, [data[k] for k in t], comments


def window(t, rows, t0, t1):
    return [r for tt, r in zip(t, rows) if t0 < tt <= t1]


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--window", type=float, default=0.5)
    a = ap.parse_args()
    case, W = a.case, a.window
    info = json.load(open(os.path.join(case, "mesh_info.json")))
    ci = dict(l.split(" ", 1) for l in open(os.path.join(case, "case_info.txt")).read().splitlines() if " " in l)
    zline = next(l for l in open(os.path.join(case, "case_info.txt")) if l.startswith("wall probe z_mm"))
    zs = [float(z) for z in zline.split()[3:]]

    tp, pw, _ = read_series(case, "wallProbes", "p")
    ta, aw, _ = read_series(case, "wallProbes", "alpha.water")
    tb, bw, _ = read_series(case, "boreWater", "volFieldValue.dat")
    ti, pin, _ = read_series(case, "inletAvg", "surfaceFieldValue.dat")
    tx, ax, _ = read_series(case, "axisProbes", "alpha.water")
    tend = max(tp[-1], ta[-1])
    t1, t0, tm = tend, tend - W, tend - 2 * W

    def stats(t, rows, col, lo, hi):
        v = [r[col] for r in window(t, rows, lo, hi)]
        return dict(mean=mean(v), min=min(v) if v else None, max=max(v) if v else None)

    r = dict(case=case, geometry=info["geometry"], Q_lpm=float(ci["Q_lpm"]), t_end_s=tend,
             avg_window_s=[t0, t1], n_cells=info["n_cells"])
    port = stats(tp, pw, 0, t0, t1)
    port_prev = stats(tp, pw, 0, tm, t0)
    r["port_wall_p_kPa_mean"] = port["mean"] / 1000
    r["port_wall_p_kPa_min"] = port["min"] / 1000
    r["port_wall_p_kPa_max"] = port["max"] / 1000
    r["port_wall_p_kPa_mean_prev_window"] = port_prev["mean"] / 1000
    r["port_wall_alpha_mean"] = stats(ta, aw, 0, t0, t1)["mean"]
    r["port_wall_alpha_min"] = stats(ta, aw, 0, t0, t1)["min"]
    prof = []
    for k, z in enumerate(zs):
        prof.append(dict(z_mm=z, alpha_mean=stats(ta, aw, k, t0, t1)["mean"],
                         alpha_min=stats(ta, aw, k, t0, t1)["min"],
                         p_kPa_mean=stats(tp, pw, k, t0, t1)["mean"] / 1000))
    r["wall_probe_profile"] = prof
    wet = [q for q in prof if q["alpha_mean"] < 0.5]
    r["air_on_wall_highest_z_mm"] = min(q["z_mm"] for q in wet) if wet else None
    r["bore_water_fraction_mean"] = stats(tb, bw, 0, t0, t1)["mean"]
    r["bore_water_fraction_prev_window"] = stats(tb, bw, 0, tm, t0)["mean"]
    r["bore_water_fraction_min"] = stats(tb, bw, 0, t0, t1)["min"]
    pin_m = stats(ti, pin, 0, t0, t1)["mean"]
    Rin = info["inlet_radius_mm"] * 1e-3
    Uin = r["Q_lpm"] / 60000 / (3.141592653589793 * Rin ** 2)
    r["inlet_static_Pa_mean"] = pin_m
    r["level_above_inlet_cm"] = (pin_m / RHO + 0.5 * Uin ** 2) / G * 100
    r["axis_alpha_outlet_face_mean"] = stats(tx, ax, 1, t0, t1)["mean"] if ax and len(ax[0]) > 1 else None
    dp = abs(port["mean"] - port_prev["mean"]) / max(abs(port["mean"]), 50.0)
    dw = abs(r["bore_water_fraction_mean"] - r["bore_water_fraction_prev_window"])
    r["steady_check"] = dict(port_p_change_rel=dp, bore_water_change_abs=dw,
                             steady=bool(dp < 0.05 and dw < 0.01))
    suction = r["port_wall_p_kPa_mean"] < 0 and r["port_wall_alpha_mean"] > 0.9
    r["dosing_verdict"] = ("suction at port (dosing possible)" if suction else
                           "no suction at port (dosing stops)")
    with open(os.path.join(case, "results_vof.json"), "w") as fh:
        json.dump(r, fh, indent=1)
    print(f"{info['geometry']} Q={r['Q_lpm']:g}: t_end={tend:.3f}s port p={r['port_wall_p_kPa_mean']:.3f} kPa "
          f"(prev {r['port_wall_p_kPa_mean_prev_window']:.3f}), port alpha={r['port_wall_alpha_mean']:.3f}, "
          f"bore water={r['bore_water_fraction_mean']:.4f}, air on wall up to z={r['air_on_wall_highest_z_mm']}, "
          f"level={r['level_above_inlet_cm']:.2f} cm, steady={r['steady_check']['steady']}, {r['dosing_verdict']}")


if __name__ == "__main__":
    main()
