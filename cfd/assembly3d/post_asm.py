#!/usr/bin/env python3
"""
post_asm.py <case dir> -- Stage 1 quantities for a 3-D assembly case from the
VTK surfaces written by post_asm.sh. Writes <case>/results.json.

Pressures: OpenFOAM p is kinematic (m2/s2) and includes the hydrostatic part
(gravity enters as a body force), gauge relative to the atmospheric outlet.
Levels (cm of water) = p / g * 100 at a cross-section, area-averaged.

Reported:
  * level_inlet_cm       area-averaged TOTAL pressure on the inlet patch, as water
                         level above the inlet centre (the overflow number, for this
                         assembly with its inlet extension)
  * level_spigot_cm      the same at the E-Sig spigot entry plane
  * level_venturi_in_cm  the same at the Venturi inlet face (compare with 2-D)
  * static levels and a negative-static flag
  * part losses: drop in mechanical energy head (p/g + V^2/2g + z) between
    sections, with flux-weighted values (energy flux / volume flux)
  * port wall pressure (kPa gauge): mean, min and max around the circumference
    at the port station; the plane of the bend is y (+y = side of the incoming
    E-Sig leg, -y = the opposite side)
  * diffuser recovery vs ideal, wall shear sign along the Venturi wall
  * outlet KE flux as head and KE factor alpha
  * mass imbalance, y+ statistics, final residuals
"""
import glob, json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtkio import read_vtk, face_geom, latest_dir  # noqa: E402

G, RHO = 9.81, 997.0
case = sys.argv[1]
info = json.load(open(os.path.join(case, "assembly_info.json")))
pp = os.path.join(case, "postProcessing")
latest = lambda func: latest_dir(pp, func)


def section(name, keep=None, func="sections"):
    fn = glob.glob(os.path.join(latest(func), f"{name}.vtk"))[0]
    pts, polys, cd, _ = read_vtk(fn)
    A, C = face_geom(pts, polys)
    m = np.ones(len(polys), bool) if keep is None else keep(C)
    return dict(A=A[m], C=C[m], p=cd["p"][m], U=cd["U"][m], area=np.linalg.norm(A[m], axis=1))


def stats(s, normal, z_up_axis=0):
    """Area-averaged static/total head, flux-weighted mechanical energy head (m)."""
    a = s["area"]; n = np.asarray(normal, float)
    un = s["U"] @ n
    U2 = (s["U"] ** 2).sum(1)
    Q = (un * a).sum()
    z = s["C"][:, z_up_axis]                          # elevation (m), +x is up
    p_avg = (s["p"] * a).sum() / a.sum()
    pt_avg = ((s["p"] + 0.5 * U2) * a).sum() / a.sum()
    E = ((s["p"] + 0.5 * U2 + G * z) * un * a).sum() / Q      # flux-weighted energy (m2/s2)
    KE = (0.5 * U2 * un * a).sum() / Q
    Vm = Q / a.sum()
    return dict(Q_m3s=2 * Q, area_m2=2 * a.sum(), V_mean=Vm, p_head_m=p_avg / G, pt_head_m=pt_avg / G,
                energy_head_m=E / G, ke_head_m=KE / G, alpha=KE / (0.5 * Vm ** 2) if Vm else float("nan"),
                z_centroid_m=float((z * a).sum() / a.sum()))


mm = 1e-3
ven_y = lambda C: np.abs(C[:, 1]) < 12 * mm          # keep the Venturi/elbow-leg region of x-normal planes
S = {
    "inlet": stats(section("inlet"), (0, -1, 0)),
    "spigot": stats(section("spigot"), (0, -1, 0)),
    "esig_end": stats(section("esigEnd"), (0, -1, 0)),
    "venturi_in": stats(section("ventIn", ven_y), (-1, 0, 0)),
    "port": stats(section("port", ven_y), (-1, 0, 0)),
    "diffuser_end": stats(section("diffEnd", ven_y), (-1, 0, 0)),
    "outlet": stats(section("outlet"), (-1, 0, 0)),
}
# levels: total head referenced to each section's own centroid elevation
lv = lambda k: S[k]["pt_head_m"] * 100
res = dict(
    level_inlet_cm=lv("inlet"), level_spigot_cm=lv("spigot"), level_venturi_in_cm=lv("venturi_in"),
    static_inlet_cm=S["inlet"]["p_head_m"] * 100, static_venturi_in_cm=S["venturi_in"]["p_head_m"] * 100,
    static_inlet_negative_FLAG=bool(S["inlet"]["p_head_m"] < 0),
    static_venturi_in_negative_FLAG=bool(S["venturi_in"]["p_head_m"] < 0),
)
E = {k: S[k]["energy_head_m"] for k in S}
res["loss_cm"] = {
    "inlet_extension+spigot": (E["inlet"] - E["spigot"]) * 100,
    "esig_reducer": (E["spigot"] - E["esig_end"]) * 100,
    "elbow_incl_thread_cavity": (E["esig_end"] - E["venturi_in"]) * 100,
    "venturi_to_port": (E["venturi_in"] - E["port"]) * 100,
    "port_to_diffuser_end": (E["port"] - E["diffuser_end"]) * 100,
    "diffuser_end_to_outlet": (E["diffuser_end"] - E["outlet"]) * 100,
    "venturi_total": (E["venturi_in"] - E["outlet"]) * 100,
}
res["outlet_ke_head_cm"] = S["outlet"]["ke_head_m"] * 100
res["outlet_alpha"] = S["outlet"]["alpha"]
res["venturi_inlet_alpha"] = S["venturi_in"]["alpha"]
res["flow_in_Lpm"] = -S["inlet"]["Q_m3s"] * 6e4 if S["inlet"]["Q_m3s"] < 0 else S["inlet"]["Q_m3s"] * 6e4
res["flow_out_Lpm"] = S["outlet"]["Q_m3s"] * 6e4
res["mass_imbalance_pct"] = (res["flow_in_Lpm"] - res["flow_out_Lpm"]) / res["flow_in_Lpm"] * 100

# velocity non-uniformity at Venturi inlet and port section: (max-min)/mean of axial velocity
for k, nm in (("venturi_in", "ventIn"), ("port", "port")):
    s = section(nm, ven_y)
    ua = -s["U"][:, 0]
    res[f"{k}_ux_max_over_mean"] = float(ua.max() / S[k]["V_mean"])
    res[f"{k}_ux_cov"] = float(np.sqrt(((ua - ua.mean()) ** 2 * s["area"]).sum() / s["area"].sum()) / S[k]["V_mean"])

# walls
W = {}
for nm in ("wall_esig", "wall_elbow", "wall_venturi"):
    fn = glob.glob(os.path.join(latest("walls"), f"{nm}.vtk"))[0]
    pts, polys, cd, _ = read_vtk(fn)
    A, C = face_geom(pts, polys)
    W[nm] = dict(C=C, area=np.linalg.norm(A, axis=1), p=cd["p"], tau=cd["wallShearStress"], yp=cd["yPlus"])
res["yplus"] = {nm: dict(mean=float((w["yp"] * w["area"]).sum() / w["area"].sum()), min=float(w["yp"].min()),
                         max=float(w["yp"].max()), p95=float(np.percentile(w["yp"], 95))) for nm, w in W.items()}
wv = W["wall_venturi"]
x_in = info["venturi_inlet_x"] * mm
zb = (x_in - wv["C"][:, 0]) / mm                     # mm below the Venturi inlet face
port_zb = info["venturi_inlet_x"] - info["port_x"]
band = np.abs(zb - port_zb) < 0.75
theta = np.degrees(np.arctan2(wv["C"][band, 2], wv["C"][band, 1]))  # 0 = +y (E-Sig side), 180 = -y
pp_kPa = wv["p"][band] * RHO / 1000
res["port_wall_p_kPa"] = dict(mean=float((pp_kPa * wv["area"][band]).sum() / wv["area"][band].sum()),
                              min=float(pp_kPa.min()), max=float(pp_kPa.max()),
                              at_plus_y=float(pp_kPa[np.argmin(np.abs(theta))]),
                              at_minus_y=float(pp_kPa[np.argmin(np.abs(np.abs(theta) - 180))]),
                              at_plus_z=float(pp_kPa[np.argmin(np.abs(theta - 90))]),
                              n_faces=int(band.sum()))
# axial wall shear along the Venturi: OpenFOAM wallShearStress is the force of the wall on
# the fluid per area (kinematic); forward-flow friction points upstream (+x here), so the
# forward-flow shear is -(-tau_x) ... flow is -x; wall acts on fluid in +x for attached flow.
tau_fwd = wv["tau"][:, 0] * RHO                      # >0 attached (wall drags fluid back, +x), <0 reversed
wall = info["venturi_wall"]
dif0 = [z for z, r in wall if z <= 60.36][-1]
dif1 = 79.4
bins = np.arange(0, wall[-1][0] + 2.5, 2.5)
prof = []
for a, b in zip(bins, bins[1:]):
    m = (zb >= a) & (zb < b)
    if m.sum():
        prof.append(dict(z_mm=float(0.5 * (a + b)), tau_fwd_Pa=float((tau_fwd[m] * wv["area"][m]).sum() / wv["area"][m].sum()),
                         frac_reversed=float(wv["area"][m][tau_fwd[m] < 0].sum() / wv["area"][m].sum()),
                         p_kPa=float((wv["p"][m] * wv["area"][m]).sum() / wv["area"][m].sum() * RHO / 1000)))
res["venturi_wall_profile"] = prof
dm = (zb >= dif0) & (zb <= wall[-1][0])
res["diffuser_reversed_area_frac"] = float(wv["area"][dm][tau_fwd[dm] < 0].sum() / wv["area"][dm].sum())
rev = [q["z_mm"] for q in prof if q["frac_reversed"] > 0.05 and q["z_mm"] >= dif0]
res["diffuser_reversal_starts_mm_below_inlet"] = rev[0] if rev else None
# elbow: reversed area fraction (separation behind the bend / steps)
we = W["wall_elbow"]
res["elbow_reversed_area_frac_note"] = "see viz; elbow wall shear direction varies with the bend, not tabulated"
# diffuser recovery vs ideal (port section -> diffuser end), area-averaged static
dp_real = (S["diffuser_end"]["p_head_m"] - S["port"]["p_head_m"]) - (S["port"]["z_centroid_m"] - S["diffuser_end"]["z_centroid_m"])
A1, A2 = S["port"]["area_m2"], S["diffuser_end"]["area_m2"]
dp_ideal = S["port"]["V_mean"] ** 2 / (2 * G) * (1 - (A1 / A2) ** 2)
res["diffuser_recovery_frac"] = dp_real / dp_ideal
res["sections"] = S
# residuals
rf = sorted(glob.glob(os.path.join(pp, "residuals*", "*", "residuals.dat")))
if rf:
    last = [l for l in open(rf[-1]) if l.strip() and not l.startswith("#")][-1].split()
    res["final_iteration"] = int(float(last[0]))
    res["final_residuals"] = dict(zip(["p", "Ux", "Uy", "Uz", "k", "omega"],
                                      [float(v) if v != "N/A" else None for v in last[1:7]]))
# inlet static-pressure drift over the last 500 iterations (convergence check)
af = sorted(glob.glob(os.path.join(pp, "inletAvg", "*", "*.dat")))
if af:
    rows = [l.split() for l in open(af[-1]) if l.strip() and not l.startswith("#")]
    pin = [(float(r[0]), float(r[1])) for r in rows]
    tail = [p for t, p in pin if t > pin[-1][0] - 500]
    res["inlet_p_last500_range_cm"] = (max(tail) - min(tail)) / G * 100   # kinematic p -> cm
json.dump(res, open(os.path.join(case, "results.json"), "w"), indent=1)
keys = ["level_inlet_cm", "level_spigot_cm", "level_venturi_in_cm", "static_inlet_cm", "static_venturi_in_cm",
        "outlet_ke_head_cm", "outlet_alpha", "diffuser_recovery_frac", "diffuser_reversed_area_frac",
        "flow_in_Lpm", "flow_out_Lpm", "mass_imbalance_pct"]
for k in keys:
    print(f"{k:32s} {res[k]:.4g}" if isinstance(res[k], float) else f"{k:32s} {res[k]}")
print("losses cm:", {k: round(v, 2) for k, v in res["loss_cm"].items()})
print("port wall p kPa:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["port_wall_p_kPa"].items()})
print("y+:", {k: {a: round(b, 1) for a, b in v.items()} for k, v in res["yplus"].items()})
print("final:", res.get("final_iteration"), res.get("final_residuals"))
