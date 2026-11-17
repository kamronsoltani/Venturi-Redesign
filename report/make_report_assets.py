#!/usr/bin/env python3
"""
make_report_assets.py -- figures (PDF) and LaTeX tables for the CFD report,
generated from the same data as the viewer (cfd/results/gui/).

Run from report/ (WSL):  python3 make_report_assets.py
Writes figures/*.pdf and tables/*.tex. Re-run after cfd/build_gui_data.py.
"""
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors

HERE = os.path.dirname(os.path.abspath(__file__))
GUI = os.path.join(HERE, "..", "cfd", "results", "gui")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

RHO, G, NU = 997.0, 9.81, 8.9e-7
IDX = json.load(open(os.path.join(GUI, "index.json")))["cases"]
BY = {c["key"]: c for c in IDX}
_cache = {}


def case(key):
    if key not in _cache:
        _cache[key] = json.load(open(os.path.join(GUI, "cases", key + ".json")))
    return _cache[key]


def get(geom, q, level="medium"):
    for c in IDX:
        if c["geometry"] == geom and abs(c["Q_lpm"] - q) < 0.01 and c["mesh_level"] == level:
            return c
    return None


def flows(geom):
    return sorted({c["Q_lpm"] for c in IDX if c["geometry"] == geom and c["mesh_level"] == "medium"})


VSET = ["V1", "V2", "V3", "V4"]
PROTO = ["A_3-4in", "A_1in", "B_3-4in", "B_1in"]
ALLG = ["current"] + PROTO + VSET
TEXNAME = {"current": "Current", "A_1in": "A\\_1in", "A_3-4in": "A\\_3-4in", "B_1in": "B\\_1in", "B_3-4in": "B\\_3-4in",
           "V1": "V1", "V2": "V2", "V3": "V3", "V4": "V4"}
# report palette (matches the template: Navy, Blue, plus distinguishable accents)
COL = {"current": "#5A626E", "A_1in": "#8FA8C8", "A_3-4in": "#B9C6D6", "B_1in": "#D6780A", "B_3-4in": "#F0BE28",
       "V1": "#0B2545", "V2": "#1C5898", "V3": "#5AAA73", "V4": "#A0408C"}
DASH = {"B_1in": (0, (4, 2)), "B_3-4in": (0, (1, 1.5)), "current": (0, (6, 2, 1, 2))}

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Latin Modern Roman", "CMU Serif", "DejaVu Serif"], "mathtext.fontset": "cm",
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 1.2, "axes.grid": True, "grid.color": "#dde2e8", "grid.linewidth": 0.5,
    "axes.edgecolor": "#5A626E", "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
HALF, FULL = 3.45, 7.1


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    plt.close(fig)


def f(x, d=2):
    return "--" if x is None else f"{x:.{d}f}"


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def write(name, s):
    open(os.path.join(TAB, name + ".tex"), "w", newline="\n").write(s)


def geom_info(g):
    c = next(x for x in IDX if x["geometry"] == g)
    d = case(c["key"])
    zn, Rn = np.array(d["z_nodes"]), np.array(d["R_nodes"])
    z0, z1 = d["diffuser"]
    r0, r1 = np.interp(z0, zn, Rn), np.interp(z1, zn, Rn)
    theta = math.degrees(math.atan((r1 - r0) / (z1 - z0)))
    col = zn[-1] - d["port_z"]
    return dict(Din=2 * Rn[0], Dt=2 * Rn.min(), Dex=2 * Rn[-1], L=zn[-1], port=d["port_z"], theta=theta,
                col=col, zd=(z0, z1), zn=zn, Rn=Rn)


GI = {g: geom_info(g) for g in ALLG}


def dyn_port(c):
    """port wall pressure with the hydrostatic column below the port removed (kPa)"""
    return c["port_wall_p_kPa_orifice_avg"] + RHO * G * GI[c["geometry"]]["col"] * 1e-3 / 1000


# ======================================================================== tables
# geometry
s = ""
for g in ALLG:
    i = GI[g]
    s += (f"{TEXNAME[g]} & {i['Din']:.1f} & {i['Dt']:.1f} & {i['Dex']:.1f} & {(i['Dt']/i['Din'])**2:.3f} & "
          f"{(i['Dex']/i['Dt'])**2:.2f} & {i['theta']:.1f} & {i['L']:.1f} & {i['port']:.2f} & {i['col']:.1f}\\\\ \\hline\n")
write("geometry", s)

# Reynolds numbers at the throat and inlet
s = ""
for g in ALLG:
    i = GI[g]
    re = lambda q, D: 4 * q / 60000 / (math.pi * D * 1e-3 * NU)
    s += f"{TEXNAME[g]} & " + " & ".join(f"{re(q, i['Dt'])/1000:.1f}" for q in (6, 8, 20, 30, 50)) + f" & {re(25, i['Din'])/1000:.1f}\\\\ \\hline\n"
write("reynolds", s)

# key results: level, V-set + B_1in
QV = [20, 25, 30, 40, 50]
s = ""
for g in VSET + ["B_1in"]:
    row = []
    lv = [get(g, q)["level_above_inlet_cm"] for q in QV]
    best = None
    for q, v in zip(QV, lv):
        mn = min(get(h, q)["level_above_inlet_cm"] for h in VSET)
        row.append(f"\\textbf{{{v:.1f}}}" if (g in VSET and abs(v - mn) < 1e-9) else f"{v:.1f}")
    s += f"{TEXNAME[g]} & " + " & ".join(row) + "\\\\ \\hline\n"
write("key_level", s)

# zero-crossing flow (linear interpolation)
s = []
for g in VSET + ["B_1in"]:
    qs = flows(g)
    lv = [get(g, q)["level_above_inlet_cm"] for q in qs]
    qz = None
    for k in range(len(qs) - 1):
        if lv[k] < 0 <= lv[k + 1]:
            qz = qs[k] + (qs[k + 1] - qs[k]) * (-lv[k]) / (lv[k + 1] - lv[k])
    s.append(f"{qz:.1f}" if qz else "--")
write("zero_cross", " / ".join(s))

# mesh table + study
lv_names = [("coarse", 54, 0.300, 0.060, 4.5), ("medium", 80, 0.200, 0.040, 3.0), ("fine", 120, 0.1333, 0.0267, 2.0)]
s = ""
for name, nr, hc, hk, h1 in lv_names:
    c = get("current", 25, name)
    s += f"{name.capitalize()} & {nr} & {hc:.3f} & {hk:.4f} & {h1:.1f} & {c['cells']:,} & {c['yplus_max']:.2f} & {c['yplus_mean']:.2f} & {c['checkMesh_max_nonortho_deg']:.1f} & {c['checkMesh_max_skewness']:.2f}\\\\ \\hline\n"
write("mesh_levels", s.replace(",", "{,}"))

QTY = [("level_above_inlet_cm", "Level above inlet (cm)", 2), ("inlet_static_Pa", "Inlet static pressure (Pa)", 0),
       ("port_wall_p_kPa_orifice_avg", "Port wall pressure (kPa)", 3), ("min_wall_p_kPa", "Minimum wall pressure (kPa)", 3),
       ("diffuser_recovery_frac_of_ideal", "Diffuser recovery (--)", 3), ("outlet_KE_head_cm", "Outlet jet head (cm)", 2),
       ("outlet_KE_alpha", "Outlet $\\alpha$ (--)", 3)]
s = ""
r = 1.5
for k, lab, d in QTY:
    f3, f2, f1 = (get("current", 25, n)[k] for n in ("coarse", "medium", "fine"))
    e21, e32 = f1 - f2, f2 - f3
    if e21 * e32 > 0 and abs(e21) > 0:
        p = math.log(abs(e32 / e21)) / math.log(r)
        gci = 1.25 * abs(e21 / f1) / (r ** p - 1) * 100
        ps, gs = f"{p:.2f}", f"{gci:.2f}"
    else:
        ps, gs = "osc.", f"{1.25 * max(abs(e21), abs(e32)) / abs(f1) * 100:.2f}$^\\dagger$"
    s += f"{lab} & {f3:.{d}f} & {f2:.{d}f} & {f1:.{d}f} & {100*(f2-f3)/abs(f2):+.2f} & {100*(f1-f2)/abs(f2):+.2f} & {ps} & {gs}\\\\ \\hline\n"
write("mesh_study", s)

# 1-D comparison (current part)
ONED = {6: (None, -0.90), 8: (None, -1.07), 20: (8, -3.2), 25: (19, -4.7), 30: (32, -6.5)}
s = ""
for q, (l1, p1) in ONED.items():
    c = get("current", q)
    lv, pp = c["level_above_inlet_cm"], c["port_wall_p_kPa_orifice_avg"]
    dl = "--" if l1 is None else f"{100*(lv-l1)/l1:+.0f}"
    dp = 100 * (pp - p1) / abs(p1)
    flag = "\\textbf{yes}" if (l1 is not None and abs(lv - l1) / l1 > 0.15) or abs(dp) > 15 else "no"
    s += f"{q} & {'below inlet' if l1 is None else l1} & {lv:.1f} & {dl} & {p1:.2f} & {pp:.2f} & {dp:+.0f} & {flag}\\\\ \\hline\n"
write("oned", s)

# full result rows for a geometry
def result_rows(geoms, qs=None):
    s = ""
    for g in geoms:
        for q in (qs or flows(g)):
            c = get(g, q)
            if c is None:
                continue
            fl = "\\cellcolor{StProg}" if c["inlet_static_negative_FLAG"] else ""
            s += (f"{TEXNAME[g]} & {q:g} & {c['level_above_inlet_cm']:.2f} & {fl}{c['inlet_static_Pa']:.0f} & {c['port_wall_p_kPa_orifice_avg']:.3f} & "
                  f"{c['diffuser_recovery_frac_of_ideal']:.3f} & {100*c['diffuser_reversed_wall_length_frac']:.0f} & {c['outlet_KE_head_cm']:.2f} & "
                  f"{c['outlet_KE_alpha']:.2f} & {c['internal_loss_cm']:.2f}\\\\ \\hline\n")
    return s


write("res_current", result_rows(["current"]))
write("res_proto", result_rows(PROTO, [20, 25, 30]))
write("res_vset", result_rows(VSET + ["B_1in"], QV))

# head budget V-set
s = ""
for g in VSET + ["B_1in"]:
    for q in QV:
        c = get(g, q)
        tot = c["level_above_inlet_cm"] + c["drop_cm"]
        s += (f"{TEXNAME[g]} & {q} & {c['drop_cm']:.1f} & {c['level_above_inlet_cm']:+.2f} & {tot:.2f} & {c['outlet_KE_head_cm']:.2f} & "
              f"{100*c['outlet_KE_head_cm']/tot:.0f} & {c['internal_loss_cm']:.2f} & {100*c['internal_loss_cm']/tot:.0f} & {c['outlet_KE_alpha']:.2f}\\\\ \\hline\n")
write("head_budget", s)

# diffuser behaviour
s = ""
for g in VSET + ["B_1in"]:
    for q in QV:
        c = get(g, q)
        rr = c["reversed_wall_flow_z_ranges_mm"]
        rs = ", ".join(f"{a:.1f}--{b:.1f}" for a, b in rr) if rr else "none"
        s += (f"{TEXNAME[g]} & {q} & {c['diffuser_Cp']:.3f} & {c['diffuser_Cp_ideal']:.3f} & {c['diffuser_recovery_frac_of_ideal']:.3f} & "
              f"{c['diffuser_wall_shear_min_Pa']:.2f} & {100*c['diffuser_reversed_wall_length_frac']:.0f} & {rs}\\\\ \\hline\n")
write("diffuser", s)

# port pressure split
s = ""
for g in VSET + ["B_1in"]:
    hs = -RHO * G * GI[g]["col"] * 1e-3 / 1000
    vals = []
    for q in QV:
        c = get(g, q)
        vals.append(f"{c['port_wall_p_kPa_orifice_avg']:.2f} & {dyn_port(c):.2f}")
    s += f"{TEXNAME[g]} & {GI[g]['col']:.1f} & {hs:.2f} & " + " & ".join(vals) + "\\\\ \\hline\n"
write("port_split", s)

# low-flow port pressure for designs run at 6-8 L/min
s = ""
for g in ["current"] + PROTO:
    hs = -RHO * G * GI[g]["col"] * 1e-3 / 1000
    c6, c8 = get(g, 6), get(g, 8)
    s += (f"{TEXNAME[g]} & {GI[g]['col']:.1f} & {hs:.2f} & {c6['port_wall_p_kPa_orifice_avg']:.2f} & {dyn_port(c6):.2f} & "
          f"{c8['port_wall_p_kPa_orifice_avg']:.2f} & {dyn_port(c8):.2f}\\\\ \\hline\n")
write("port_lowflow", s)

# V1 vs B_1in
s = ""
for q in (25, 30):
    a, b = get("V1", q), get("B_1in", q)
    for k, lab, d in [("level_above_inlet_cm", "Level above inlet (cm)", 2), ("port_wall_p_kPa_orifice_avg", "Port wall pressure (kPa)", 3),
                      ("diffuser_recovery_frac_of_ideal", "Diffuser recovery (--)", 3), ("diffuser_reversed_wall_length_frac", "Reversed diffuser wall (fraction)", 2),
                      ("outlet_KE_head_cm", "Outlet jet head (cm)", 2), ("outlet_KE_alpha", "Outlet $\\alpha$", 2), ("internal_loss_cm", "Internal loss (cm)", 2),
                      ("drop_cm", "Drop inlet to outlet (cm)", 1)]:
        s += f"{q} & {lab} & {b[k]:.{d}f} & {a[k]:.{d}f} & {a[k]-b[k]:+.{d}f}\\\\ \\hline\n"
write("v1_vs_b", s)

# convergence summary per geometry
s = ""
for g in ALLG:
    cs = [c for c in IDX if c["geometry"] == g]
    pr = max(c["res_final"]["p"]["final"] for c in cs)
    mi = max(abs(c["mass_imbalance_rel"]) for c in cs)
    dp = max(abs(c["inlet_p_change_last500_rel"]) for c in cs)
    yp = max(c["yplus_max"] for c in cs)
    s += f"{TEXNAME[g]} & {len(cs)} & {sci(pr)} & {sci(dp)} & {sci(mi)} & {yp:.2f}\\\\ \\hline\n"
write("convergence", s)

# appendix: every case
s = ""
for g in ALLG:
    for c in sorted([c for c in IDX if c["geometry"] == g], key=lambda c: (c["Q_lpm"], c["mesh_level"] != "medium", c["mesh_level"])):
        tag = "" if c["mesh_level"] == "medium" else f" ({c['mesh_level']})"
        fl = "\\cellcolor{StProg}" if c["inlet_static_negative_FLAG"] else ""
        s += (f"{TEXNAME[g]}{tag} & {c['Q_lpm']:g} & {c['cells']} & {c['level_above_inlet_cm']:.2f} & {c['inlet_total_Pa']:.0f} & {fl}{c['inlet_static_Pa']:.0f} & "
              f"{c['port_wall_p_kPa_orifice_avg']:.3f} & {c['min_wall_p_kPa']:.3f} & {c['diffuser_Cp']:.3f} & {c['diffuser_recovery_frac_of_ideal']:.3f} & "
              f"{100*c['diffuser_reversed_wall_length_frac']:.0f} & {c['outlet_KE_head_cm']:.2f} & {c['outlet_KE_alpha']:.2f} & {c['internal_loss_cm']:.2f} & "
              f"{c['yplus_max']:.2f} & {c['res_final']['p']['final']:.1e} & {c['mass_imbalance_rel']:.0e}\\\\ \\hline\n")
write("appendix_all", s)

# numbers quoted in the text
vals = {}
for g in VSET + ["B_1in"]:
    for q in QV:
        c = get(g, q)
        vals[f"{g}_{q}"] = c
open(os.path.join(TAB, "numbers.txt"), "w").write(json.dumps({k: {kk: v[kk] for kk in ("level_above_inlet_cm", "port_wall_p_kPa_orifice_avg", "outlet_KE_alpha", "internal_loss_cm", "outlet_KE_head_cm", "diffuser_reversed_wall_length_frac")} for k, v in vals.items()}, indent=0))

# ======================================================================== figures
# bore profiles
fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.5), constrained_layout=True)
for g in ["current", "A_1in", "B_1in"]:
    i = GI[g]; axs[0].plot(i["zn"], i["Rn"] * 2, color=COL[g], ls=DASH.get(g, "-"), label=TEXNAME[g].replace("\\", ""))
for g in VSET:
    i = GI[g]; axs[1].plot(i["zn"], i["Rn"] * 2, color=COL[g], label=g)
for a, t in zip(axs, ["Current part and Stage 3 prototypes", "Printed set V1--V4"]):
    a.set_xlabel("$z$ below inlet face (mm)"); a.set_ylabel("bore diameter $D$ (mm)"); a.set_title(t.replace("--", "–"), fontsize=8); a.legend(frameon=False)
save(fig, "bore_profiles")

# level and port vs flow
def sweep(key, ylabel, name, geoms, zero=True, transform=None):
    fig, ax = plt.subplots(figsize=(HALF, 2.5), constrained_layout=True)
    for g in geoms:
        qs = flows(g)
        ys = [get(g, q)[key] if transform is None else transform(get(g, q)) for q in qs]
        ax.plot(qs, ys, marker="o", ms=2.8, color=COL[g], ls=DASH.get(g, "-"), label=TEXNAME[g].replace("\\", ""))
    if zero:
        ax.axhline(0, color="#5A626E", lw=0.7)
    ax.set_xlabel("flow rate $Q$ (L/min)"); ax.set_ylabel(ylabel); ax.legend(frameon=False, ncol=2)
    save(fig, name)


sweep("level_above_inlet_cm", "level above inlet face $h$ (cm)", "level_vs_q", VSET + ["B_1in", "current"])
sweep("port_wall_p_kPa_orifice_avg", "port wall pressure $p_{\\mathrm{port}}$ (kPa)", "port_vs_q", VSET + ["B_1in", "current"])
sweep(None, "Venturi part $p_{\\mathrm{port}}+\\rho g h_{\\mathrm{col}}$ (kPa)", "port_dyn_vs_q", VSET + ["B_1in", "current"], transform=dyn_port)
sweep("outlet_KE_alpha", "outlet energy factor $\\alpha$ (--)", "alpha_vs_q", VSET + ["B_1in"], zero=False)

# wall shear and wall pressure along the wall at 30 L/min
for key, ylabel, name, ylim in [("wall_tau_Pa", "axial wall shear $\\tau_{w,x}$ (Pa)", "wall_shear_30", (-2, 20)),
                                ("wall_p_kPa", "wall static pressure (kPa)", "wall_p_30", None)]:
    fig, ax = plt.subplots(figsize=(FULL, 2.4), constrained_layout=True)
    for g in VSET + ["B_1in"]:
        d = case(get(g, 30)["key"])
        zc = 0.5 * (np.array(d["z_nodes"][1:]) + np.array(d["z_nodes"][:-1]))
        ax.plot(zc, d[key], color=COL[g], ls=DASH.get(g, "-"), label=TEXNAME[g].replace("\\", ""))
        ax.axvline(d["port_z"], color=COL[g], lw=0.5, ls=":")
    ax.axhline(0, color="#5A626E", lw=0.7)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel("$z$ below inlet face (mm)"); ax.set_ylabel(ylabel); ax.legend(frameon=False, ncol=5)
    save(fig, name)

# outlet velocity profiles at 30 L/min
fig, ax = plt.subplots(figsize=(HALF, 2.5), constrained_layout=True)
for g in VSET + ["B_1in"]:
    c = get(g, 30); d = case(c["key"])
    e = np.array(d["cut_eta"]); u = np.array(d["cut_Ux"][-1])
    um = c["Q_lpm"] / 60000 / (math.pi * (d["R_nodes"][-1] * 1e-3) ** 2)
    ax.plot(np.r_[-e[::-1], e], np.r_[u[::-1], u] / um, color=COL[g], ls=DASH.get(g, "-"), label=TEXNAME[g].replace("\\", ""))
ax.axhline(0, color="#5A626E", lw=0.7); ax.set_xlim(-1, 1)
ax.set_xlabel("$r/R$ at the outlet"); ax.set_ylabel("$u_x/\\bar{V}$"); ax.legend(frameon=False, ncol=2)
save(fig, "outlet_profiles_30")

# head budget bars at 30 and 50
fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.3), constrained_layout=True, sharey=True)
gs = VSET + ["B_1in"]
for ax, q in zip(axs, (30, 50)):
    jet = [get(g, q)["outlet_KE_head_cm"] for g in gs]; loss = [get(g, q)["internal_loss_cm"] for g in gs]
    drop = [get(g, q)["drop_cm"] for g in gs]
    x = np.arange(len(gs))
    ax.bar(x, jet, 0.6, color="#1C5898", label="outlet jet $\\alpha \\bar V_o^2/2g$")
    ax.bar(x, loss, 0.6, bottom=jet, color="#F0BE28", label="internal loss $h_L$")
    ax.scatter(x, drop, marker="_", s=260, color="#0B2545", zorder=3, label="drop $\\Delta z$ (level $=$ bar top $-$ tick)")
    ax.set_xticks(x); ax.set_xticklabels([g.replace("_", "\\_") if False else g for g in gs]); ax.set_title(f"{q} L/min", fontsize=8)
    ax.set_ylabel("head (cm)")
axs[0].legend(frameon=False, fontsize=6.5, loc="upper left")
save(fig, "head_budget")

# speed contours at 30 L/min, five designs
fig, axs = plt.subplots(1, 5, figsize=(FULL, 3.6), constrained_layout=True, sharey=True)
smax = max(max(map(max, case(get(g, 30)["key"])["cut_speed"])) for g in gs)
norm = colors.Normalize(0, smax)
for ax, g in zip(axs, gs):
    d = case(get(g, 30)["key"])
    zc = np.array(d["cut_z"]); Rc = np.array(d["cut_R"]); e = np.array(d["cut_eta"])
    Z = np.repeat(zc[:, None], len(e), 1); Rr = Rc[:, None] * e[None, :]
    S = np.array(d["cut_speed"]); U = np.array(d["cut_Ux"])
    Zm = np.concatenate([Z[:, ::-1], Z], 1); Rm = np.concatenate([-Rr[:, ::-1], Rr], 1)
    cf = ax.contourf(Rm, Zm, np.concatenate([S[:, ::-1], S], 1), levels=24, cmap="viridis", norm=norm)
    if (U < 0).any():
        ax.contourf(Rm, Zm, np.concatenate([U[:, ::-1], U], 1), levels=[U.min() - 1, 0], colors=["#E0503C"], alpha=0.85)
    for L in d["streamlines_rz"]:
        L = np.array(L); ax.plot(L[:, 0], L[:, 1], color="w", lw=0.35); ax.plot(-L[:, 0], L[:, 1], color="w", lw=0.35)
    zn, Rn = np.array(d["z_nodes"]), np.array(d["R_nodes"])
    ax.plot(Rn, zn, color="#0B2545", lw=0.8); ax.plot(-Rn, zn, color="#0B2545", lw=0.8)
    ax.plot([-Rn[np.searchsorted(zn, d['port_z'])] - 3, -Rn[np.searchsorted(zn, d['port_z'])]], [d["port_z"]] * 2, color="#D6780A", lw=1.6)
    ax.set_aspect("equal"); ax.set_title(g, fontsize=8); ax.set_xlabel("$r$ (mm)"); ax.grid(False)
axs[0].set_ylabel("$z$ below inlet face (mm)"); axs[0].invert_yaxis()
cb = fig.colorbar(plt.cm.ScalarMappable(norm, "viridis"), ax=axs, shrink=0.7, pad=0.01); cb.set_label("$|\\mathbf{U}|$ (m/s)")
save(fig, "speed_cuts_30")

# mesh study: wall pressure and level on three meshes
fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.3), constrained_layout=True)
for n, ls in (("coarse", ":"), ("medium", "-"), ("fine", "--")):
    d = case(get("current", 25, n)["key"])
    zc = 0.5 * (np.array(d["z_nodes"][1:]) + np.array(d["z_nodes"][:-1]))
    axs[0].plot(zc, d["wall_p_kPa"], ls=ls, color="#0B2545", label=n)
    axs[1].plot(zc, d["wall_tau_Pa"], ls=ls, color="#1C5898", label=n)
axs[0].set_ylabel("wall pressure (kPa)"); axs[1].set_ylabel("axial wall shear (Pa)")
for a in axs:
    a.set_xlabel("$z$ below inlet face (mm)"); a.legend(frameon=False)
save(fig, "mesh_study")

# convergence history, V1 at 30 L/min
d = case(get("V1", 30)["key"])
fig, axs = plt.subplots(1, 2, figsize=(FULL, 2.2), constrained_layout=True)
for k, c in zip(("p", "Ux", "Uy", "k", "omega"), ("#0B2545", "#1C5898", "#8FA8C8", "#D6780A", "#5AAA73")):
    axs[0].semilogy(d["residuals"]["it"], d["residuals"][k], color=c, label={"omega": "$\\omega$", "Ux": "$U_x$", "Uy": "$U_r$"}.get(k, k))
axs[0].set_xlabel("iteration"); axs[0].set_ylabel("initial residual"); axs[0].legend(frameon=False, ncol=5)
axs[1].plot(d["inlet_p_hist"]["it"], d["inlet_p_hist"]["level_cm"], color="#0B2545")
axs[1].set_ylim(-5, 10); axs[1].set_xlabel("iteration"); axs[1].set_ylabel("inlet static $\\bar p_1/\\rho g$ (cm)")
save(fig, "convergence_V1_30")
print("figures:", sorted(os.listdir(FIG)))
print("tables:", sorted(os.listdir(TAB)))
