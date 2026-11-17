"""review_assembly.py <fluid dir> <out.png> [title]

Review picture of a water volume from build_assembly.py / build_fluid.py, in the
installed orientation (up is up, the pump is to the right):
  left : the bend (symmetry) plane, water coloured by part, with flow and gravity arrows
  right: 3-D view of the walls (half domain mirrored), coloured by part
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

d, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(os.path.dirname(os.path.abspath(d)))
info = json.load(open(os.path.join(d, "assembly_info.json")))
sol, cur = {}, None
for l in open(os.path.join(d, "fluid_mm.stl")):
    l = l.strip()
    if l.startswith("solid"):
        cur = l.split()[1]; sol.setdefault(cur, [])
    elif l.startswith("vertex"):
        sol[cur].append([float(v) for v in l.split()[1:]])
tri = {k: np.array(v).reshape(-1, 3, 3) for k, v in sol.items() if v}

COL = {"wall_esig": "#4c78a8", "wall_elbow": "#f2a03d", "wall_venturi": "#8e6bbf", "inlet": "#54a24b", "outlet": "#e45756"}
up = info.get("upstream", "esig34")
NAME = {"wall_esig": "long-radius reducer, straight part" if up == "lr1" else "E-Sig reducer",
        "wall_elbow": "reducing bend" if up == "lr1" else "elbow", "wall_venturi": "Venturi"}
vx, split = info["venturi_inlet_x"], info["esig_end_face_y"]


def part_of(c):
    if c[0] < vx - 0.05:
        return "wall_venturi"
    return "wall_esig" if c[1] > split + 0.05 else "wall_elbow"


fig = plt.figure(figsize=(22, 8.5))
ax = fig.add_subplot(1, 3, 1)
S = tri["symmetry"]
cen = S.mean(1)
for k in ("wall_esig", "wall_elbow", "wall_venturi"):
    m = np.array([part_of(c) == k for c in cen])
    ax.add_collection(PolyCollection(S[m][:, :, [1, 0]], facecolor=COL[k], edgecolor="none", alpha=0.55, label=NAME[k]))
for k in ("inlet", "outlet"):
    T = tri[k]
    ax.plot(T[:, :, 1].ravel(), T[:, :, 0].ravel(), ".", ms=1, color=COL[k])
ax.autoscale(); ax.set_aspect("equal"); ax.grid(alpha=0.3)
ax.set_xlabel("horizontal, toward the pump (mm)"); ax.set_ylabel("height (mm)")
ax.set_title("bend plane (water), installed orientation")
xi, yi = info["inlet_centre_x"], info["inlet_y"]
ax.annotate("from pump", xy=(yi - 25, xi), xytext=(yi + 15, xi), color=COL["inlet"], va="center",
            arrowprops=dict(arrowstyle="->", color=COL["inlet"], lw=2))
ax.annotate("to outlet\n(open air)", xy=(0, info["outlet_x"] - 22), xytext=(28, info["outlet_x"] - 18),
            color=COL["outlet"], ha="left", va="center",
            arrowprops=dict(arrowstyle="<-", color=COL["outlet"], lw=2))
ax.annotate("", xy=(0, info["outlet_x"] - 30), xytext=(0, info["outlet_x"] - 2),
            arrowprops=dict(arrowstyle="->", color=COL["outlet"], lw=2))
x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
ax.annotate("up", xy=(x0 + 8, y1 - 8), xytext=(x0 + 8, y1 - 40), ha="center",
            arrowprops=dict(arrowstyle="->", lw=1.5))
ax.set_ylim(y0 - 35, y1)
ax.legend(loc="lower right", fontsize=9)
txt = [f"Venturi inlet face at height {vx:.1f} mm", f"outlet at {info['outlet_x']:.1f} mm",
       f"thread engagement {info.get('engage_mm', 14.0):.1f} mm", f"inlet centreline at {xi:.1f} mm"]
ax.text(0.02, 0.02, "\n".join(txt), transform=ax.transAxes, fontsize=8, va="bottom",
        bbox=dict(fc="white", ec="0.7"))

# close-up of the elbow / Venturi inlet joint
az = fig.add_subplot(1, 3, 2)
for k in ("wall_esig", "wall_elbow", "wall_venturi"):
    m = np.array([part_of(c) == k for c in cen])
    az.add_collection(PolyCollection(S[m][:, :, [1, 0]], facecolor=COL[k], edgecolor="k", linewidth=0.05, alpha=0.6))
for k in ("wall_esig", "wall_elbow", "wall_venturi"):
    T = tri[k]; e = T[np.abs(T[:, :, 2]).max(1) < 0.4]            # wall triangles touching the bend plane
    az.plot(e[:, :, 1].ravel(), e[:, :, 0].ravel(), ".", ms=1.5, color=COL[k])
cx0 = info["esig_end_face_y"]
az.set_xlim(-30, max(cx0 + 25, 60)); az.set_ylim(vx - 45, max(info["inlet_centre_x"], 0) + 30)
az.set_aspect("equal"); az.grid(alpha=0.3)
az.axhline(vx, color=COL["wall_venturi"], lw=0.8, ls="--"); az.text(-28, vx + 1, "Venturi inlet face", fontsize=8)
az.axvline(cx0, color=COL["wall_esig"], lw=0.8, ls="--")
az.text(cx0 + 1, vx - 40, "reducer straight / bend" if up == "lr1" else "E-Sig end face", fontsize=8, rotation=90)
az.set_title("close-up: elbow and Venturi inlet (bend plane)")
az.set_xlabel("horizontal, toward the pump (mm)")
ax3 = fig.add_subplot(1, 3, 3, projection="3d")
for k in ("wall_esig", "wall_elbow", "wall_venturi", "inlet", "outlet"):
    if k not in tri:
        continue
    T = tri[k]
    for sgn in (1, -1):
        Q = T.copy(); Q[:, :, 2] *= sgn
        P = Q[:, :, [1, 2, 0]]          # plot (toward pump, across, up)
        ax3.add_collection3d(Poly3DCollection(P[::2], facecolor=COL[k], edgecolor="none", alpha=0.35 if k.startswith("wall") else 0.9))
allp = np.vstack([t.reshape(-1, 3) for t in tri.values()])
rng = np.ptp(allp, axis=0).max() / 2
mid = [(allp[:, 1].min() + allp[:, 1].max()) / 2, 0, (allp[:, 0].min() + allp[:, 0].max()) / 2]
ax3.set_xlim(mid[0] - rng, mid[0] + rng); ax3.set_ylim(-rng, rng); ax3.set_zlim(mid[2] - rng, mid[2] + rng)
ax3.set_xlabel("toward pump"); ax3.set_ylabel("across"); ax3.set_zlabel("up")
ax3.view_init(elev=18, azim=-62)
ax3.set_title("3-D (half domain mirrored); green inlet, red outlet")
fig.suptitle(title, fontsize=13)
fig.text(0.5, 0.01, "Assumptions: " + "; ".join(info.get("assumptions", [])), ha="center", fontsize=7, wrap=True)
plt.tight_layout(rect=(0, 0.04, 1, 0.96))
plt.savefig(out, dpi=90)
