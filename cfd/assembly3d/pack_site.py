#!/usr/bin/env python3
"""
pack_site.py <viewer json dir> <results csv> <out docs dir>
Packs the 3-D viewer files (docs/3d/cases/*.json from build_viz3d.py) into compact binary files for the
combined site:
  data/geom/<asm>.bin      wall + centre-plane geometry (identical for every flow of an assembly)
  data/flow/<asm>_Q<q>.bin fields, streamlines and the case's result tables
  data/site.json           catalogue: assemblies, flows, per-case summary numbers from the results CSV
Binary layout: uint32 header length, JSON header, then arrays (each 4-byte aligned, little endian).
uint16 arrays are quantised: value = lo + q * (hi - lo) / 65535 (per component for xyz arrays).
"""
import csv, glob, json, os, re, struct, sys
import numpy as np

src, csv_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(os.path.join(out, "data", "geom"), exist_ok=True)
os.makedirs(os.path.join(out, "data", "flow"), exist_ok=True)

ASM = {  # id: (upstream group, venturi name, short label)
    "asm34_current": ("e34", "Current", "3/4in E-Sig + 3/4in elbow + current Venturi"),
    "asm34_A":       ("e34", "A",       "3/4in E-Sig + 3/4in elbow + Venturi A (14 mm throat, 19 mm exit)"),
    "asm34_Bvane":   ("e34", "B + vane", "3/4in E-Sig + 3/4in elbow + Venturi B (14 mm throat, 26 mm exit) + cross-vane"),
    "asm1e_V1":  ("e1", "V1", "1in E-Sig + McMaster 1in elbow + V1 (14 mm throat, 26 mm exit)"),
    "asm1e_V2":  ("e1", "V2", "1in E-Sig + McMaster 1in elbow + V2 (16 mm throat, 26 mm exit)"),
    "asm1e_V3":  ("e1", "V3", "1in E-Sig + McMaster 1in elbow + V3 (18 mm throat, 26 mm exit)"),
    "asm1e_V4":  ("e1", "V4", "1in E-Sig + McMaster 1in elbow + V4 (14 mm throat, 30 mm exit)"),
    "asm1lr_V1": ("lr", "V1", "1in long-radius reducing elbow + V1 (14 mm throat, 26 mm exit)"),
    "asm1lr_V2": ("lr", "V2", "1in long-radius reducing elbow + V2 (16 mm throat, 26 mm exit)"),
    "asm1lr_V3": ("lr", "V3", "1in long-radius reducing elbow + V3 (18 mm throat, 26 mm exit)"),
    "asm1lr_V4": ("lr", "V4", "1in long-radius reducing elbow + V4 (14 mm throat, 30 mm exit)"),
    # aligned-mesh pair (symmetry plane exactly on the centreline): compare these two with each other
    "asm34_current_zfix": ("p34", "Threaded elbow", "3/4in E-Sig + 3/4in threaded elbow + current Venturi (aligned mesh)"),
    "asm34lr_current":    ("p34", "Long-radius 42 mm", "3/4in E-Sig + 3/4in long-radius elbow (R 42 mm, no reduction) + current Venturi (aligned mesh)"),
}
GROUPS = {"e34": "3/4in E-Sig + 3/4in elbow", "e1": "1in E-Sig + McMaster elbow", "lr": "1in long-radius elbow",
          "p34": "3/4in elbow retrofit, aligned mesh"}


def write_bin(path, meta, arrays):
    hdr = dict(meta=meta, arrays=[])
    blobs, off = [], 0
    for name, a, kind in arrays:
        a = np.asarray(a)
        ent = dict(name=name, n=int(a.shape[0]), comp=int(a.shape[1]) if a.ndim == 2 else 1)
        if kind == "q16":
            a2 = a.reshape(len(a), -1).astype(np.float64)
            lo, hi = a2.min(0), a2.max(0)
            span = np.where(hi > lo, hi - lo, 1.0)
            q = np.round((a2 - lo) / span * 65535).astype("<u2")
            ent.update(type="u16", lo=lo.tolist(), hi=(lo + span).tolist()); b = q.tobytes()
        elif kind == "u32":
            b = a.astype("<u4").tobytes(); ent["type"] = "u32"
        elif kind == "f32":
            b = a.astype("<f4").tobytes(); ent["type"] = "f32"
        ent["offset"] = off
        pad = (-len(b)) % 4
        blobs.append(b + b"\0" * pad); off += len(b) + pad
        hdr["arrays"].append(ent)
    h = json.dumps(hdr, separators=(",", ":")).encode()
    h += b" " * ((-len(h)) % 4)
    with open(path, "wb") as f:
        f.write(struct.pack("<I", len(h))); f.write(h)
        for b in blobs:
            f.write(b)
    return os.path.getsize(path)


rows = {}
for r in csv.DictReader(open(csv_path)):
    if r["mesh_level"] != "medium" or r["assembly"] not in ASM:
        continue
    rows[(r["assembly"], int(float(r["Q_Lpm"])))] = r


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


site = dict(groups=GROUPS, assemblies=[], cases={})
total = 0
for aid, (grp, vname, label) in ASM.items():
    files = sorted(glob.glob(os.path.join(src, f"{aid}_Q*.json")))
    files = [f for f in files if re.search(rf"{aid}_Q(\d+)\.json$", f)]
    flows = sorted(int(re.search(r"_Q(\d+)\.json$", f).group(1)) for f in files)
    geo_done = False
    for q in flows:
        d = json.load(open(os.path.join(src, f"{aid}_Q{q}.json")))
        if not geo_done:
            W, C = d["wall"], d["centre"]
            total += write_bin(os.path.join(out, "data", "geom", f"{aid}.bin"),
                               dict(asm=aid, frame=d["frame"], stations=d["stations"]),
                               [("wall_pos", np.reshape(W["pos"], (-1, 3)), "q16"), ("wall_tri", W["tri"], "u32"),
                                ("centre_pos", np.reshape(C["pos"], (-1, 3)), "q16"), ("centre_tri", C["tri"], "u32")])
            geo_done = True
        SL = d["streamlines"]
        lens = [len(s["s"]) for s in SL]
        P = np.concatenate([np.reshape(s["p"], (-1, 3)) for s in SL])
        S = np.concatenate([s["s"] for s in SL]); T = np.concatenate([s["t"] for s in SL])
        meta = {k: d[k] for k in ("label", "Q_lpm", "results", "losses", "port", "yplus", "final_residuals",
                                  "final_iteration", "assumptions")}
        meta["label"] = f"{label}, {q} L/min"
        total += write_bin(os.path.join(out, "data", "flow", f"{aid}_Q{q}.bin"), meta,
                           [("wall_p", d["wall"]["p"], "q16"), ("wall_tau", d["wall"]["tau"], "q16"),
                            ("centre_p", d["centre"]["p"], "q16"), ("centre_spd", d["centre"]["spd"], "q16"),
                            ("sl_len", lens, "u32"), ("sl_pos", P, "q16"), ("sl_spd", S, "q16"), ("sl_t", T, "f32")])
        r = rows.get((aid, q), {})
        site["cases"][f"{aid}_Q{q}"] = dict(
            level=fnum(r.get("level_inlet_centre_cm")), drift=fnum(r.get("inlet_level_drift_last500_cm")),
            port=fnum(r.get("port_p_mean_kPa")), jet=fnum(r.get("outlet_ke_head_cm")),
            recovery=fnum(r.get("diffuser_recovery_frac")), imbalance=fnum(r.get("mass_imbalance_pct")))
    site["assemblies"].append(dict(id=aid, group=grp, venturi=vname, label=label, flows=flows))
    print(aid, flows)

json.dump(site, open(os.path.join(out, "data", "site.json"), "w"), indent=1)
print("total binary MB", round(total / 1e6, 1), "cases", len(site["cases"]))
