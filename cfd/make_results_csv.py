#!/usr/bin/env python3
"""
make_results_csv.py -- collect results.json of finished cases into one CSV,
with the 1-D model predictions (analysis/turitap_venturi_analysis.m) alongside (current part only).

Usage: python3 make_results_csv.py <geom> <out.csv> "<case dir>" ["<case dir>" ...]
"""
import csv
import json
import os
import sys

# 1-D model predictions from analysis/turitap_venturi_analysis.m (current part, Venturi only, no elbow).
# None = "below inlet" (the model gives no positive level).
ONE_D = {
    "current": {6: (None, -0.90), 8: (None, -1.07), 20: (8, -3.2),
                25: (19, -4.7), 30: (32, -6.5)},
}

COLS = ["geometry", "level", "Q_lpm", "cells", "iterations", "residualControl_met",
        "converged_judgement", "inlet_p_change_last500_rel", "res_p_final", "res_Ux_final",
        "res_k_final", "res_omega_final", "mass_imbalance_rel",
        "level_above_inlet_cm", "inlet_total_Pa", "inlet_static_Pa", "inlet_static_negative_FLAG",
        "port_wall_p_kPa_orifice_avg", "port_wall_p_kPa_at_axis", "min_wall_p_kPa", "min_wall_p_z_mm",
        "diffuser_Cp", "diffuser_Cp_ideal", "diffuser_recovery_frac_of_ideal",
        "diffuser_separated", "diffuser_reversed_wall_fraction", "diffuser_wall_shear_min_Pa",
        "diffuser_reversed_wall_length_frac", "reversed_wall_flow_z_ranges_mm",
        "outlet_KE_head_cm", "outlet_KE_alpha", "outlet_drop_below_inlet_cm", "internal_loss_cm",
        "yplus_max", "yplus_mean", "checkMesh_ok", "checkMesh_max_nonortho_deg",
        "oneD_level_cm", "level_diff_pct", "oneD_port_kPa", "port_diff_pct", "flag_gt15pct"]


def judge(r):
    """Plain converged / not converged call used in the CSV."""
    res = r.get("residuals", {})
    p_ok = res.get("p", {}).get("final", 1) < 1e-4
    stable = r.get("inlet_p_change_last500_rel", 1) < 1e-3
    mass = abs(r.get("mass_imbalance_rel", 1)) < 1e-4
    if r.get("residualControl_met") and mass:
        return "converged (residual targets met)"
    if p_ok and stable and mass:
        return "converged (residuals flat, quantities stable <0.1%/500 it)"
    return "NOT converged"


def row(r, geom):
    out = {c: r.get(c) for c in COLS}
    res = r.get("residuals", {})
    out["res_p_final"] = res.get("p", {}).get("final")
    out["res_Ux_final"] = res.get("Ux", {}).get("final")
    out["res_k_final"] = res.get("k", {}).get("final")
    out["res_omega_final"] = res.get("omega", {}).get("final")
    out["Q_lpm"] = round(r["Q_lpm"], 3)
    out["converged_judgement"] = judge(r)
    # energy split: level + drop = outlet jet KE head + internal loss
    # (outlet static p = 0, the outlet face lies 'drop' below the inlet face)
    mi = r.get("_mesh_info")
    if mi:
        drop = (mi["stations"]["outlet"] - mi["stations"]["inlet"]) / 10.0
        out["outlet_drop_below_inlet_cm"] = drop
        out["internal_loss_cm"] = r["level_above_inlet_cm"] + drop - r["outlet_KE_head_cm"]
    q = round(r["Q_lpm"])
    lvl, port = ONE_D.get(geom, {}).get(q, (None, None))
    flags = []
    out["oneD_level_cm"] = "below inlet" if (q in ONE_D.get(geom, {}) and lvl is None) else lvl
    if lvl is not None:
        d = (r["level_above_inlet_cm"] - lvl) / abs(lvl) * 100
        out["level_diff_pct"] = round(d, 1)
        if abs(d) > 15:
            flags.append("level")
    elif q in ONE_D.get(geom, {}) and r["level_above_inlet_cm"] > 0:
        flags.append("level (CFD above inlet, 1-D below)")
    out["oneD_port_kPa"] = port
    if port is not None:
        d = (r["port_wall_p_kPa_orifice_avg"] - port) / abs(port) * 100
        out["port_diff_pct"] = round(d, 1)
        if abs(d) > 15:
            flags.append("port")
    out["flag_gt15pct"] = "; ".join(flags)
    for k, v in out.items():
        if isinstance(v, float):
            out[k] = f"{v:.5g}"
    return out


def main():
    geom, path, cases = sys.argv[1], sys.argv[2], sys.argv[3:]
    rows = []
    for c in cases:
        r = json.load(open(os.path.join(c, "results.json")))
        r["_mesh_info"] = json.load(open(os.path.join(c, "mesh_info.json")))
        rows.append(row(r, geom))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {path}")


if __name__ == "__main__":
    main()
