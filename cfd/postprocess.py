#!/usr/bin/env python3
"""
postprocess.py -- Stage 1 tracked quantities for one finished case.

Reads mesh_info.json, the latest time directory (p, U, wallShearStress,
yPlus), postProcessing/ and log.foamRun / log.checkMesh. No OpenFOAM needed.
Writes results.json in the case folder and prints a short summary.

Pressures: OpenFOAM's p is kinematic (p/rho); converted with rho = 997.
p includes gravity (buoyancyForce, g along +x), outlet static = 0 gauge.

Usage: python3 postprocess.py "<case dir>"
"""
import glob
import json
import math
import os
import re
import sys

RHO = 997.0
G = 9.81
ORIFICE_D_MM = 1.67


# ------------------------------------------------------------ file readers
def _parse_values(s):
    s = s.strip()
    if s.startswith("uniform"):
        v = s[len("uniform"):].split(";", 1)[0].strip()
        return ("uniform", _parse_one(v))
    m = re.match(r"nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    if not m:
        raise ValueError(s[:80])
    n = int(m.group(2))
    body = s[m.end():]
    if m.group(1) == "scalar":
        vals = [float(x) for x in body.split(None, n)[:n]]
    else:
        vals = [tuple(map(float, t.split()))
                for t in re.findall(r"\(([^()]*)\)", body[: body.find("\n)") + 2])][:n]
    if len(vals) != n:
        raise ValueError("short list")
    return ("list", vals)


def _parse_one(v):
    v = v.strip()
    if v.startswith("("):
        return tuple(map(float, v.strip("()").split()))
    return float(v)


def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    internal = _parse_values(txt[i + len("internalField"):])
    patches = {}
    bf = txt[txt.index("boundaryField"):]
    for m in re.finditer(r"\n    (\w+)\n    \{(.*?)\n    \}", bf, re.S):
        blk = m.group(2)
        vm = re.search(r"\n        value\s+", blk)    # patch-level entry only
        patches[m.group(1)] = _parse_values(blk[vm.end():]) if vm else None
    return internal, patches


def expand(v, n):
    kind, val = v
    return [val] * n if kind == "uniform" else val


def last_line(path):
    rows = [l for l in open(path) if l.strip() and not l.startswith("#")]
    return rows[-1].split() if rows else None


def pp_file(case, name):
    f = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "*.dat")))
    return f[-1] if f else None


def pp_series(path):
    hdr, rows = None, []
    for l in open(path):
        if l.startswith("#"):
            hdr = l[1:].split()
        elif l.strip():
            rows.append(l.split())
    return hdr, rows


# --------------------------------------------------------------- analysis
def analyse(case):
    info = json.load(open(os.path.join(case, "mesh_info.json")))
    NX, NR = info["NX"], info["NR"]
    zn = [z * 1e-3 for z in info["z_nodes_mm"]]
    Rn = [r * 1e-3 for r in info["R_nodes_mm"]]
    eta = info["eta"]
    st = info["stations"]
    times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"]
    if not times:
        raise SystemExit("no solution time directory")
    tdir = os.path.join(case, max(times, key=float))
    it_final = int(float(os.path.basename(tdir)))

    p_int, p_bc = read_field(os.path.join(tdir, "p"))
    U_int, U_bc = read_field(os.path.join(tdir, "U"))
    p_c = expand(p_int, NX * NR)
    U_c = expand(U_int, NX * NR)
    cell = lambda i, j: i * NR + j

    def band_w(i):     # cross-section area weights of radial bands at node i
        return [(eta[j + 1] * Rn[i]) ** 2 - (eta[j] * Rn[i]) ** 2 for j in range(NR)]

    out = dict(case=os.path.relpath(case), geometry=info["geometry"], level=info["level"],
               cells=info["n_cells"], iterations=it_final)

    # flow rate (from the inlet BC, full 360 deg)
    U_in = expand(U_bc["inlet"], NR)
    Rin, Rout = Rn[0], Rn[-1]
    w_in = band_w(0)
    Ubar_in = sum(u[0] * w for u, w in zip(U_in, w_in)) / sum(w_in)
    Q = Ubar_in * math.pi * Rin ** 2
    out["Q_lpm"] = Q * 60000

    # 1, 2: inlet total and static pressure (area averages)
    p_in = expand(p_bc["inlet"], NR)
    ps = sum(p * w for p, w in zip(p_in, w_in)) / sum(w_in)
    pt = sum((p + 0.5 * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2)) * w
             for p, u, w in zip(p_in, U_in, w_in)) / sum(w_in)
    out["inlet_total_Pa"] = pt * RHO
    out["level_above_inlet_cm"] = pt / G * 100
    out["inlet_static_Pa"] = ps * RHO
    out["inlet_static_negative_FLAG"] = ps < 0

    # 3: wall pressure at the chlorine port
    p_w = expand(p_bc["wall"], NX)
    zf = [0.5 * (zn[i] + zn[i + 1]) for i in range(NX)]
    zp = st["port"] * 1e-3
    a, b = zp - 0.5 * ORIFICE_D_MM * 1e-3, zp + 0.5 * ORIFICE_D_MM * 1e-3
    num = den = 0.0
    for i in range(NX):
        ov = max(0.0, min(zn[i + 1], b) - max(zn[i], a))
        num += p_w[i] * ov
        den += ov

    def interp(xs, ys, x):
        for k in range(len(xs) - 1):
            if xs[k] <= x <= xs[k + 1]:
                return ys[k] + (ys[k + 1] - ys[k]) * (x - xs[k]) / (xs[k + 1] - xs[k])
        return ys[0] if x < xs[0] else ys[-1]

    out["port_z_mm"] = st["port"]
    out["port_wall_p_kPa_orifice_avg"] = num / den * RHO / 1000
    out["port_wall_p_kPa_at_axis"] = interp(zf, p_w, zp) * RHO / 1000
    imin = min(range(NX), key=lambda i: p_w[i])
    out["min_wall_p_kPa"] = p_w[imin] * RHO / 1000
    out["min_wall_p_z_mm"] = zf[imin] * 1000

    # 4: diffuser pressure recovery (section-averaged piezometric p = p - g x)
    def section_p(z):
        i = min(range(len(zn)), key=lambda k: abs(zn[k] - z))
        cols = [c for c in (i - 1, i) if 0 <= c < NX]
        vals, xs = [], []
        for c in cols:
            w = [((eta[j + 1] ** 2 - eta[j] ** 2)) for j in range(NR)]
            pc = sum(p_c[cell(c, j)] * w[j] for j in range(NR)) / sum(w)
            xc = 0.5 * (zn[c] + zn[c + 1])
            vals.append(pc - G * xc)
            xs.append(xc)
        if len(vals) == 1:
            return vals[0], i
        return vals[0] + (vals[1] - vals[0]) * (zn[i] - xs[0]) / (xs[1] - xs[0]), i

    if st["diffuser_start"] is not None and st["diffuser_end"] is not None:
        p1, i1 = section_p(st["diffuser_start"] * 1e-3)
        p2, i2 = section_p(st["diffuser_end"] * 1e-3)
        V1 = Q / (math.pi * Rn[i1] ** 2)
        cp = (p2 - p1) / (0.5 * V1 ** 2)
        cpi = 1 - (Rn[i1] / Rn[i2]) ** 4
        out["diffuser_Cp"] = cp
        out["diffuser_Cp_ideal"] = cpi
        out["diffuser_recovery_frac_of_ideal"] = cp / cpi

    # 5: axial wall shear sign along the diffuser (and anywhere on the wall)
    ut = []
    for i in range(NX):
        tz, tr = zn[i + 1] - zn[i], Rn[i + 1] - Rn[i]
        n = math.hypot(tz, tr)
        u = U_c[cell(i, NR - 1)]
        ut.append((u[0] * tz + u[1] * tr) / n)
    tau = None
    wss_path = os.path.join(tdir, "wallShearStress")
    if os.path.exists(wss_path):
        _, wbc = read_field(wss_path)
        tw = expand(wbc["wall"], NX)
        tau = []
        for i in range(NX):
            tz, tr = zn[i + 1] - zn[i], Rn[i + 1] - Rn[i]
            n = math.hypot(tz, tr)
            tau.append((tw[i][0] * tz + tw[i][1] * tr) / n * RHO)
        # orient so that positive = same direction as the near-wall flow in
        # the most attached cell (OpenFOAM's sign convention is fluid-on-wall
        # vs wall-on-fluid depending on version)
        iref = max(range(NX), key=lambda i: abs(ut[i]))
        if tau[iref] * ut[iref] < 0:
            tau = [-t for t in tau]
    d0, d1 = st["diffuser_start"] * 1e-3, st["diffuser_end"] * 1e-3
    idiff = [i for i in range(NX) if d0 <= zf[i] <= d1]
    rev_d = [i for i in idiff if ut[i] < 0]
    rev_all = [i for i in range(NX) if ut[i] < 0]
    out["diffuser_wall_shear_min_Pa"] = min(tau[i] for i in idiff) if tau else None
    out["diffuser_wall_shear_max_Pa"] = max(tau[i] for i in idiff) if tau else None
    out["diffuser_separated"] = bool(rev_d)
    out["diffuser_reversed_wall_fraction"] = len(rev_d) / len(idiff)
    dzf = [zn[i + 1] - zn[i] for i in range(NX)]
    out["diffuser_reversed_wall_length_frac"] = sum(dzf[i] for i in rev_d) / sum(dzf[i] for i in idiff)

    def ranges(idx):
        rs, start, prev = [], None, None
        for i in idx:
            if start is None:
                start = prev = i
            elif i == prev + 1:
                prev = i
            else:
                rs.append((start, prev)); start = prev = i
        if start is not None:
            rs.append((start, prev))
        return [(round(zn[a] * 1000, 2), round(zn[b + 1] * 1000, 2)) for a, b in rs]

    out["reversed_wall_flow_z_ranges_mm"] = ranges(rev_all)

    # 6: kinetic-energy flux leaving the outlet, as head
    U_out = expand(U_bc["outlet"], NR)
    w_out = band_w(NX)
    q_out = sum(u[0] * w for u, w in zip(U_out, w_out))
    ke = sum(0.5 * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) * u[0] * w for u, w in zip(U_out, w_out))
    Vout = Q / (math.pi * Rout ** 2)
    out["outlet_KE_head_cm"] = ke / q_out / G * 100
    out["outlet_KE_alpha"] = (ke / q_out) / (0.5 * Vout ** 2)
    out["outlet_backflow_faces"] = sum(1 for u in U_out if u[0] < 0)

    # 7: mass imbalance and residuals
    fi, fo = pp_file(case, "inletFlow"), pp_file(case, "outletFlow")
    if fi and fo:
        qi, qo = float(last_line(fi)[1]), float(last_line(fo)[1])
        out["mass_imbalance_rel"] = (qi + qo) / abs(qi)
    rf = pp_file(case, "residuals*")
    if rf:
        hdr, rows = pp_series(rf)
        rows = [r for r in rows if r[1] != "N/A"]
        cols = hdr[1:]
        res = {}
        for k, c in enumerate(cols, start=1):
            vals = [float(r[k]) for r in rows if r[k] != "N/A"]
            if vals:
                res[c] = dict(first=vals[0], final=vals[-1], max_last500=max(vals[-500:]),
                              min_last500=min(vals[-500:]))
        out["residuals"] = res
        # Uz on a wedge is ~0, so its normalised residual is noise (not used)
    log = open(os.path.join(case, "log.foamRun")).read() if os.path.exists(
        os.path.join(case, "log.foamRun")) else ""
    m = re.search(r"SIMPLE solution converged in (\S+) iterations", log)
    out["residualControl_met"] = bool(m)
    out["ended_cleanly"] = "\nEnd" in log
    ia = pp_file(case, "inletAvg")
    if ia:
        hdr, rows = pp_series(ia)
        k = hdr.index("areaAverage(p)")
        series = [float(r[k]) for r in rows]
        n = min(500, len(series) - 1)
        out["inlet_p_change_last500_rel"] = abs(series[-1] - series[-1 - n]) / max(abs(series[-1]), 1e-12)

    # y+ on the wall
    yp_path = os.path.join(tdir, "yPlus")
    if os.path.exists(yp_path):
        _, ybc = read_field(yp_path)
        yp = expand(ybc["wall"], NX)
        out["yplus_max"] = max(yp)
        out["yplus_mean"] = sum(yp) / len(yp)
        out["yplus_max_z_mm"] = zf[max(range(NX), key=lambda i: yp[i])] * 1000

    # checkMesh summary
    cm = os.path.join(case, "log.checkMesh")
    if os.path.exists(cm):
        t = open(cm).read()
        m = re.search(r"Mesh non-orthogonality Max: ([0-9.eE+-]+)", t)
        out["checkMesh_ok"] = "Mesh OK" in t
        out["checkMesh_max_nonortho_deg"] = float(m.group(1)) if m else None
        m = re.search(r"Max skewness = ([0-9.eE+-]+)", t)
        out["checkMesh_max_skewness"] = float(m.group(1)) if m else None

    # wall profiles for plotting / later checks
    out["_wall_profile"] = dict(z_mm=[z * 1000 for z in zf],
                                p_kPa=[p * RHO / 1000 for p in p_w],
                                u_t_nearwall=ut, tau_t_Pa=tau)
    return out


def main():
    case = sys.argv[1]
    out = analyse(case)
    with open(os.path.join(case, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k, v in out.items():
        if k.startswith("_"):
            continue
        if isinstance(v, float):
            v = f"{v:.5g}"
        print(f"{k:36s} {v}")


if __name__ == "__main__":
    main()
