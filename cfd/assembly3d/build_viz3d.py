#!/usr/bin/env python3
"""
build_viz3d.py <case dir> <out json> [label]
Packs a 3-D assembly result for the web viewer (docs/3d/): the wall surface with
pressure and wall-shear sign, the centre (symmetry) plane with speed and
pressure, and inlet-seeded streamlines with a travel-time coordinate so the page
can animate particles. Half-domain data; the page mirrors it about z = 0.
Coordinates are written in mm in a display frame: X = horizontal (toward the
pump), Y = up, Z = across.
"""
import glob, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtkio import read_vtk, latest_dir  # noqa: E402

case, out = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(os.path.dirname(os.path.abspath(case)))
RHO = 997.0
info = json.load(open(os.path.join(case, "assembly_info.json")))
latest = lambda func: latest_dir(os.path.join(case, "postProcessing"), func)
res = json.load(open(os.path.join(case, "results.json")))


def to_disp(P):
    """elbow frame (m): x up, y toward pump, z across  ->  display mm (X toward pump, Y up, Z across)."""
    return np.stack([P[:, 1], P[:, 0], P[:, 2]], 1) * 1000.0


def tris(polys):
    t = []
    for f in polys:
        for j in range(1, len(f) - 1):
            t.append((f[0], f[j], f[j + 1]))
    return np.array(t, int)


def surf(func, name, fields):
    fn = glob.glob(os.path.join(latest(func), f"{name}.vtk"))[0]
    pts, polys, cd, pd = read_vtk(fn)
    T = tris(polys)
    used = np.unique(T)
    remap = -np.ones(len(pts), int); remap[used] = np.arange(len(used))
    d = dict(pos=np.round(to_disp(pts[used]), 1).ravel().tolist(), tri=remap[T].ravel().tolist())
    for k, f in fields.items():
        d[k] = np.round(f(pd)[used], 3).tolist()
    return d


ex = (0, 1, 0)
wall = surf("viz", "wall", {
    "p": lambda pd: pd["p"] * RHO / 1000.0,                                   # kPa gauge
    "tau": lambda pd: np.linalg.norm(pd["wallShearStress"], axis=1) * RHO,     # Pa
})
centre = surf("viz", "centre", {
    "p": lambda pd: pd["p"] * RHO / 1000.0,
    "spd": lambda pd: np.linalg.norm(pd["U"], axis=1),
})

# streamlines
fn = sorted(glob.glob(os.path.join(latest("streamlines"), "*.vtk")))[0]
pts, lines, cd, pd = read_vtk(fn)
U = pd["U"]; spd = np.linalg.norm(U, axis=1)
P = to_disp(pts)
SL = []
for L in lines:
    L = np.array(L)
    if len(L) < 5:
        continue
    q = P[L]; s = spd[L]
    seg = np.linalg.norm(np.diff(q, axis=0), axis=1)
    keep = [0]; acc = 0.0
    for i, d in enumerate(seg, 1):          # resample to ~1.6 mm spacing
        acc += d
        if acc >= 1.6 or i == len(seg):
            keep.append(i); acc = 0.0
    q, s = q[keep], s[keep]
    seg = np.linalg.norm(np.diff(q, axis=0), axis=1) / 1000.0
    vm = 0.5 * (s[1:] + s[:-1])
    t = np.concatenate([[0.0], np.cumsum(seg / np.maximum(vm, 1e-4))])
    if t[-1] > 20:   # stuck in a recirculation; keep only the first 20 s
        n = int(np.searchsorted(t, 20))
        q, s, t = q[:n], s[:n], t[:n]
    SL.append(dict(p=np.round(q, 1).ravel().tolist(), s=np.round(s, 2).tolist(), t=np.round(t, 4).tolist()))

o = dict(
    label=label, Q_lpm=res.get("flow_in_Lpm"), frame="mm; X toward pump, Y up, Z across; half domain z>=0",
    wall=wall, centre=centre, streamlines=SL,
    stations=dict(venturi_inlet_Y=info["venturi_inlet_x"], port_Y=info["port_x"], outlet_Y=info["outlet_x"],
                  esig_end_X=info["esig_end_face_y"], spigot_X=info["spigot_entry_y"]),
    results={k: res[k] for k in ("level_inlet_cm", "level_spigot_cm", "level_venturi_in_cm", "static_inlet_cm",
                                 "static_venturi_in_cm", "outlet_ke_head_cm", "outlet_alpha", "diffuser_recovery_frac",
                                 "diffuser_reversed_area_frac", "mass_imbalance_pct", "flow_in_Lpm")},
    losses=res["loss_cm"], port=res["port_wall_p_kPa"], yplus=res["yplus"],
    wall_profile=res["venturi_wall_profile"], final_residuals=res.get("final_residuals"),
    final_iteration=res.get("final_iteration"), assumptions=info["assumptions"],
)
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
json.dump(o, open(out, "w"), separators=(",", ":"))
print(out, round(os.path.getsize(out) / 1e6, 2), "MB;", len(wall["tri"]) // 3, "wall tris,",
      len(centre["tri"]) // 3, "centre tris,", len(SL), "streamlines")
