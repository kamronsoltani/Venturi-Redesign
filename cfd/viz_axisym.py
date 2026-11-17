#!/usr/bin/env python3
"""
viz_axisym.py -- pictures of one converged 2-D axisymmetric (wedge) case.

Usage: python3 viz_axisym.py "<case dir>"   -> writes <case>/viz/*.png

Reads the latest time directory, mesh_info.json and results.json. The wedge
solution is axisymmetric, so the "centre cut" is the wedge mirrored about the
axis and the 3-D view is the wedge revolved. Streamlines are contours of the
Stokes stream function psi(z, r) = integral_0^r Ux 2 pi r dr (exact for
axisymmetric flow); reversed flow is outlined by the Ux = 0 contour.
Flow is drawn downward, as installed.
"""
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

RHO = 997.0
case = sys.argv[1]
info = json.load(open(os.path.join(case, "mesh_info.json")))
res = json.load(open(os.path.join(case, "results.json")))
NX, NR = info["NX"], info["NR"]
zn = np.array(info["z_nodes_mm"])
Rn = np.array(info["R_nodes_mm"])
eta = np.array(info["eta"])
st = info["stations"]
tdir = os.path.join(case, str(max(int(d) for d in os.listdir(case) if d.isdigit())))
out = os.path.join(case, "viz")
os.makedirs(out, exist_ok=True)
geom = info["geometry"]
q = res.get("Q_lpm") or float(open(os.path.join(case, "case_info.txt")).read().split("Q_lpm")[1].split()[0])
title = f"{geom}, {q:g} L/min"


def internal(name):
    s = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    n = int(m.group(2)); b = s[m.end():]
    if m.group(1) == "vector":
        return np.array(re.findall(r"\(([^()]+)\)", b)[:n], dtype=object)
    return np.array(b.split(")")[0].split()[:n], dtype=float)


def wall_vec(name):
    s = open(os.path.join(tdir, name)).read()
    w = s[s.index("\n    wall"):]
    m = re.search(r"value\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", w)
    n = int(m.group(1))
    return np.array([list(map(float, x.split())) for x in re.findall(r"\(([^()]+)\)", w[m.end():])[:n]])


U = np.array([list(map(float, x.split())) for x in internal("U")]).reshape(NX, NR, 3)
Ux, Ur = U[..., 0], U[..., 1]
p = internal("p").reshape(NX, NR) * RHO / 1000.0          # kPa gauge
# OpenFOAM's wallShearStress is the stress the wall puts on the fluid (opposes the
# near-wall flow); flip it so positive = forward (downward) wall flow, as in postprocess.py
tau = -wall_vec("wallShearStress")[:, 0] * RHO            # Pa, axial component
speed = np.hypot(Ux, Ur)

# node grid (axial x radial)
Zg = np.repeat(zn[:, None], NR + 1, 1)
Rg = Rn[:, None] * eta[None, :]
zc = 0.5 * (zn[1:] + zn[:-1])
Rm = 0.5 * (Rn[1:] + Rn[:-1])
rc = Rm[:, None] * 0.5 * (eta[1:] + eta[:-1])[None, :]

# stream function on nodes: psi_ij = sum_{k<j} Ux * pi (r_{k+1}^2 - r_k^2), cells -> node columns
ring = np.pi * ((Rm[:, None] * eta[None, 1:]) ** 2 - (Rm[:, None] * eta[None, :-1]) ** 2) * 1e-6
psi_c = np.concatenate([np.zeros((NX, 1)), np.cumsum(Ux * ring, 1)], 1)       # (NX, NR+1) at cell-column centres
Qm3 = psi_c[:, -1].mean()
psiN = psi_c / Qm3
Zp = np.repeat(zc[:, None], NR + 1, 1)
Rp = Rm[:, None] * eta[None, :]


def mirror(Z, R, F):
    return np.concatenate([Z[:, ::-1], Z], 1), np.concatenate([-R[:, ::-1], R], 1), np.concatenate([F[:, ::-1], F], 1)


def outline(ax):
    wz, wr = zn, Rn
    ax.plot(wr, wz, "k", lw=1.2); ax.plot(-wr, wz, "k", lw=1.2)
    ax.axhline(st["port"], color="0.4", lw=0.6, ls=":")
    ax.annotate("port", (-Rn.max() * 1.02, st["port"]), ha="right", va="center", fontsize=8, color="0.3")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xlabel("r (mm)"); ax.set_ylabel("distance below inlet face (mm)")


# ---------------------------------------------------------------- 1. speed + streamlines + reversal
fig, axs = plt.subplots(1, 2, figsize=(10, 9), constrained_layout=True)
ax = axs[0]
# pcolormesh with mirrored node grid (2*NR+1 nodes across)
Xn = np.concatenate([-Rg[:, :0:-1], Rg], 1)
Yn = np.concatenate([Zg[:, :0:-1], Zg], 1)
Sc = np.concatenate([speed[:, ::-1], speed], 1)
pc = ax.pcolormesh(Xn, Yn, Sc, shading="flat", cmap="viridis", rasterized=True)
fig.colorbar(pc, ax=ax, shrink=0.6, label="speed (m/s)")
Xp = np.concatenate([-Rp[:, :0:-1], Rp], 1); Yp = np.concatenate([Zp[:, :0:-1], Zp], 1)
Pp = np.concatenate([psiN[:, :0:-1], psiN], 1)
ax.contour(Xp, Yp, Pp, levels=np.linspace(0.05, 0.95, 10), colors="w", linewidths=0.6, alpha=0.8)
Xc = np.concatenate([-rc[:, ::-1], rc], 1); Yc = np.repeat(zc[:, None], 2 * NR, 1)
Uc = np.concatenate([Ux[:, ::-1], Ux], 1)
if (Ux < 0).any():
    ax.contour(Xc, Yc, Uc, levels=[0], colors="r", linewidths=1.2)
    ax.contourf(Xc, Yc, Uc, levels=[Uc.min() - 1, 0], colors="r", alpha=0.35)
outline(ax); ax.set_title("speed, streamlines (white)\nreversed flow (red)" + ("" if (Ux < 0).any() else ": none"))

ax = axs[1]
Pc = np.concatenate([p[:, ::-1], p], 1)
lim = np.abs(Pc).max()
pc = ax.pcolormesh(Xn, Yn, Pc, shading="flat", cmap="RdBu_r", norm=colors.TwoSlopeNorm(0, -lim, max(Pc.max(), 1e-3)), rasterized=True)
fig.colorbar(pc, ax=ax, shrink=0.6, label="static pressure (kPa gauge)")
ax.contour(Xc, Yc, np.concatenate([p[:, ::-1], p], 1), levels=[0], colors="k", linewidths=0.6, linestyles="--")
outline(ax); ax.set_title("static pressure\n(dashed: atmospheric)")
fig.suptitle(f"{title}: centre cut (flow downward)", fontsize=12)
fig.savefig(os.path.join(out, "1_centre_cut.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 2. wall pressure + axial shear along the wall
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True)
ax[0].plot(zc, p[:, -1], color="C0")
ax[0].axvline(st["port"], color="C3", ls=":"); ax[0].text(st["port"], p[:, -1].min(), " port", color="C3", va="bottom")
ax[0].axhline(0, color="0.5", lw=0.6)
ax[0].set_ylabel("wall pressure (kPa gauge)")
ax[1].plot(zc, tau, color="C2"); ax[1].axhline(0, color="0.5", lw=0.6)
ax[1].fill_between(zc, tau, 0, where=tau < 0, color="r", alpha=0.4, label="reversed (tau < 0)")
ax[1].set_ylabel("axial wall shear (Pa)"); ax[1].set_xlabel("distance below inlet face (mm)")
for a in ax:
    a.axvspan(st["diffuser_start"], st["diffuser_end"], color="0.9", zorder=0)
ax[0].text(0.5 * (st["diffuser_start"] + st["diffuser_end"]), ax[0].get_ylim()[1], "diffuser", ha="center", va="top", color="0.4")
if (tau < 0).any(): ax[1].legend()
fig.suptitle(f"{title}: along the wall", fontsize=12)
fig.savefig(os.path.join(out, "2_wall_profiles.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 3. velocity profiles at stations
fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
mid = 0.5 * (st["diffuser_start"] + st["diffuser_end"])
for zs, lab in [(1.0, "inlet (z=1)"), (st["port"], f"throat / port (z={st['port']:.0f})"),
                (mid, f"mid-diffuser (z={mid:.0f})"), (zn[-1] - 1.0, f"outlet (z={zn[-1]-1:.0f})")]:
    i = int(np.argmin(np.abs(zc - zs)))
    um = Qm3 / (np.pi * (Rm[i] * 1e-3) ** 2)
    rr = 0.5 * (eta[1:] + eta[:-1])
    ax.plot(np.r_[-rr[::-1], rr], np.r_[Ux[i, ::-1], Ux[i]] / um, label=lab)
ax.axhline(0, color="0.5", lw=0.6)
ax.set_xlabel("r / R (local)"); ax.set_ylabel("axial velocity / mean")
ax.legend(fontsize=8); ax.set_title(f"{title}: velocity profiles")
fig.savefig(os.path.join(out, "3_profiles.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 4. head split
lvl = res["level_above_inlet_cm"]; ke = res["outlet_KE_head_cm"]
drop = (st["outlet"] - st["inlet"]) / 10.0
loss = lvl + drop - ke
fig, ax = plt.subplots(figsize=(8, 3.2), constrained_layout=True)
ax.barh(["available", "spent"], [drop, 0], color="0.7", label=f"drop to outlet {drop:.1f}")
ax.barh(["available"], [lvl], left=[drop], color="C3" if lvl > 0 else "C2", label=f"water level {lvl:+.1f}")
ax.barh(["spent"], [ke], color="C0", label=f"outlet jet {ke:.1f} (alpha {res['outlet_KE_alpha']:.2f})")
ax.barh(["spent"], [loss], left=[ke], color="C1", label=f"internal loss {loss:.1f}")
ax.set_xlabel("head (cm of water)"); ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
ax.set_title(f"{title}: where the head goes")
fig.savefig(os.path.join(out, "4_head_split.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 5. data for the interactive 3-D viewer
# subsampled meridional field (axial x radial) + wall pressure + streamline paths
ia = np.unique(np.r_[np.linspace(0, NX - 1, 160).round().astype(int)])
jr = np.unique(np.r_[np.linspace(0, NR - 1, 24).round().astype(int)])
lines = []
cs = plt.figure().add_subplot().contour(Rp, Zp, psiN, levels=np.linspace(0.1, 0.9, 9))
for seg in cs.allsegs:
    for s in seg:
        if len(s) > 4:
            k = np.linspace(0, len(s) - 1, min(len(s), 120)).round().astype(int)
            lines.append([[round(float(a), 2), round(float(b), 2)] for a, b in s[k]])
plt.close("all")
view = {
    "title": title, "geom": geom, "Q_lpm": q,
    "z_nodes": [round(float(x), 2) for x in zn], "R_nodes": [round(float(x), 3) for x in Rn],
    "eta_nodes": [round(float(x), 5) for x in eta],
    "port_z": st["port"], "diffuser": [st["diffuser_start"], st["diffuser_end"]],
    "wall_p_kPa": [round(float(x), 3) for x in p[:, -1]],
    "wall_tau_Pa": [round(float(x), 3) for x in tau],
    "cut_z": [round(float(zc[i]), 2) for i in ia],
    "cut_eta": [round(float(0.5 * (eta[j] + eta[j + 1])), 4) for j in jr],
    "cut_R": [round(float(Rm[i]), 3) for i in ia],
    "cut_speed": [[round(float(speed[i, j]), 3) for j in jr] for i in ia],
    "cut_p_kPa": [[round(float(p[i, j]), 3) for j in jr] for i in ia],
    "cut_Ux": [[round(float(Ux[i, j]), 3) for j in jr] for i in ia],
    "cut_Ur": [[round(float(Ur[i, j]), 4) for j in jr] for i in ia],
    "streamlines_rz": lines,
    "results": {k: res.get(k) for k in ("level_above_inlet_cm", "inlet_static_Pa", "port_wall_p_kPa_orifice_avg",
                                         "diffuser_recovery_frac_of_ideal", "outlet_KE_head_cm", "outlet_KE_alpha",
                                         "diffuser_reversed_wall_length_frac", "yplus_max")},
    "drop_cm": drop, "internal_loss_cm": loss,
}
json.dump(view, open(os.path.join(out, "view3d.json"), "w"), separators=(",", ":"))
print("wrote", out)
