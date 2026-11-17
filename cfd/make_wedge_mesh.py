#!/usr/bin/env python3
"""
make_wedge_mesh.py -- 2-D axisymmetric wedge mesh of a Venturi bore.

Builds an OpenFOAM polyMesh (points/faces/owner/neighbour/boundary) directly
from one geometry in geometry/bore_stations.csv. No blockMesh needed.

Geometry
  * Stations (z down from the inlet face, bore diameter) are joined by
    straight lines (cones / cylinders).
  * A station whose label contains "fillet R<n>" gets a circular fillet of
    radius <n> mm, tangent to both neighbouring lines (the station is the
    sharp corner the fillet replaces).
  * Flow axis = OpenFOAM x (flow in +x, downward). Radius in the x-y plane.

Mesh
  * Structured H-grid: radial lines at constant x, nodes on every profile
    vertex (so the wall is reproduced exactly), axial clustering at sharp
    corners, smooth growth elsewhere.
  * Radial tanh clustering to the wall. The first-cell height equals h1 at
    the narrowest section and scales with local radius elsewhere (y+ ~ 1).
  * Wedge of 5 deg total (+/-2.5 deg about the x-y plane). Cells next to the
    axis are prisms, there is no axis patch.
  * Patches: inlet (patch), outlet (patch), wall (wall), front/back (wedge).

Also writes mesh_info.json in the case folder: node positions, key stations
(port, diffuser start/end) and the wall-adjacent cell-centre distance of
every wall face (used for the wall omega value).

Usage
  python3 make_wedge_mesh.py --geom current --case "<case dir>" [--level medium]
                             [--csv "<path to bore_stations.csv>"]
"""
import argparse
import csv
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "..", "geometry", "bore_stations.csv")

# Mesh levels. Lengths in mm. Refinement ratio 1.5 in every direction.
LEVELS = {
    "coarse": dict(nr=54,  h_core=0.300, h_corner=0.060, h1=4.5e-3),
    "medium": dict(nr=80,  h_core=0.200, h_corner=0.040, h1=3.0e-3),
    "fine":   dict(nr=120, h_core=0.1333, h_corner=0.0267, h1=2.0e-3),
}
AXIAL_GROWTH = 1.10          # max axial cell growth ratio away from corners
WEDGE_HALF_DEG = 2.5
SHARP_TURN_DEG = 0.1         # wall turn above this (no fillet) = sharp corner


# --------------------------------------------------------------- geometry
def load_stations(csv_path, geom):
    rows = []
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["geometry"].strip() == geom:
                rows.append(dict(z=float(row["z_mm_below_inlet_face"]),
                                 r=0.5 * float(row["bore_diameter_mm"]),
                                 label=row["station"].strip()))
    if len(rows) < 2:
        raise SystemExit(f"geometry '{geom}' not found (or <2 stations) in {csv_path}")
    rows.sort(key=lambda s: s["z"])
    return rows


def _unit(dx, dy):
    n = math.hypot(dx, dy)
    return dx / n, dy / n


def build_wall(st, h_arc):
    """Return wall polyline [(z, r, sharp)] in mm, fillets discretised."""
    pts = [(st[0]["z"], st[0]["r"], False)]
    for k in range(1, len(st) - 1):
        P = (st[k]["z"], st[k]["r"])
        A = (st[k - 1]["z"], st[k - 1]["r"])
        B = (st[k + 1]["z"], st[k + 1]["r"])
        u1 = _unit(P[0] - A[0], P[1] - A[1])
        u2 = _unit(B[0] - P[0], B[1] - P[1])
        cross = u1[0] * u2[1] - u1[1] * u2[0]
        phi = math.atan2(abs(cross), u1[0] * u2[0] + u1[1] * u2[1])  # turn angle
        m = re.search(r"fillet\s*R\s*([0-9.]+)", st[k]["label"], re.I)
        if m and phi > 1e-9:
            Rf = float(m.group(1))
            t = Rf * math.tan(phi / 2)
            la = math.hypot(P[0] - A[0], P[1] - A[1])
            lb = math.hypot(B[0] - P[0], B[1] - P[1])
            if t > 0.999 * min(la, lb):
                raise SystemExit(f"fillet R{Rf} at z={P[0]} does not fit the adjacent segments")
            s = 1.0 if cross > 0 else -1.0
            T1 = (P[0] - u1[0] * t, P[1] - u1[1] * t)
            n1 = (-s * u1[1], s * u1[0])                 # towards the centre
            C = (T1[0] + n1[0] * Rf, T1[1] + n1[1] * Rf)
            n = max(3, math.ceil(Rf * phi / h_arc))
            for q in range(n + 1):
                a = s * phi * q / n
                u = (u1[0] * math.cos(a) - u1[1] * math.sin(a),
                     u1[0] * math.sin(a) + u1[1] * math.cos(a))
                pts.append((C[0] + s * u[1] * Rf, C[1] - s * u[0] * Rf, False))
        else:
            pts.append((P[0], P[1], math.degrees(phi) > SHARP_TURN_DEG))
    pts.append((st[-1]["z"], st[-1]["r"], False))
    # sanity: strictly increasing z
    for a, b in zip(pts, pts[1:]):
        if b[0] <= a[0] + 1e-9:
            raise SystemExit("wall profile is not monotonic in z")
    return pts


def key_stations(st):
    """Port and diffuser stations from the CSV labels (z in mm)."""
    lab = [s["label"].lower() for s in st]
    z = [s["z"] for s in st]
    port = next((z[i] for i, l in enumerate(lab) if "port" in l), None)
    d0 = next((z[i] for i, l in enumerate(lab)
               if "start of diffuser" in l or "throat end" in l), None)
    d1 = next((z[i] for i, l in enumerate(lab)
               if ("diffuser" in l and ("end of" in l or "diffuser end" in l))), None)
    return dict(port=port, diffuser_start=d0, diffuser_end=d1,
                inlet=z[0], outlet=z[-1])


# ------------------------------------------------------------ distributions
def axial_nodes(wall, h_core, h_corner):
    corners = [p[0] for p in wall if p[2]]

    def h(z):
        if not corners:
            return h_core
        d = min(abs(z - c) for c in corners)
        return min(h_core, h_corner + (AXIAL_GROWTH - 1.0) * d)

    nodes = [wall[0][0]]
    for (za, _, _), (zb, _, _) in zip(wall, wall[1:]):
        m = 400
        dz = (zb - za) / m
        cum = [0.0]
        for q in range(m):
            zm = za + (q + 0.5) * dz
            cum.append(cum[-1] + dz / h(zm))
        n = max(1, round(cum[-1]))
        for k in range(1, n):
            target = cum[-1] * k / n
            q = next(i for i in range(1, m + 1) if cum[i] >= target)
            f = (target - cum[q - 1]) / (cum[q] - cum[q - 1])
            nodes.append(za + (q - 1 + f) * dz)
        nodes.append(zb)
    return nodes


def wall_radius(wall, z):
    for (za, ra, _), (zb, rb, _) in zip(wall, wall[1:]):
        if za - 1e-12 <= z <= zb + 1e-12:
            return ra + (rb - ra) * (z - za) / (zb - za)
    raise ValueError(z)


def radial_eta(nr, first_frac):
    """eta_j = tanh(d j/nr)/tanh(d), d chosen so 1-eta_{nr-1} = first_frac."""
    def first(d):
        return 1.0 - math.tanh(d * (1.0 - 1.0 / nr)) / math.tanh(d)
    if first(1e-6) < first_frac:
        return [j / nr for j in range(nr + 1)]
    lo, hi = 1e-6, 30.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if first(mid) > first_frac:
            lo = mid
        else:
            hi = mid
    d = 0.5 * (lo + hi)
    return [math.tanh(d * j / nr) / math.tanh(d) for j in range(nr + 1)]


# --------------------------------------------------------------- writing
def foam_header(cls, obj, note=None):
    s = ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
         "  Generated by cfd/make_wedge_mesh.py\n"
         "\\*---------------------------------------------------------------------------*/\n"
         "FoamFile\n{\n    format      ascii;\n"
         f"    class       {cls};\n")
    if note:
        s += f'    note        "{note}";\n'
    s += f'    location    "constant/polyMesh";\n    object      {obj};\n}}\n\n'
    return s


def newell(P):
    nx = ny = nz = 0.0
    for k in range(len(P)):
        x1, y1, z1 = P[k]
        x2, y2, z2 = P[(k + 1) % len(P)]
        nx += (y1 - y2) * (z1 + z2)
        ny += (z1 - z2) * (x1 + x2)
        nz += (x1 - x2) * (y1 + y2)
    return nx, ny, nz


def build_mesh(geom, case, level="medium", csv_path=DEFAULT_CSV):
    L = LEVELS[level]
    st = load_stations(csv_path, geom)
    wall = build_wall(st, L["h_core"])
    keys = key_stations(st)
    zn = axial_nodes(wall, L["h_core"], L["h_corner"])
    Rn = [wall_radius(wall, z) for z in zn]
    Rmin = min(Rn)
    eta = radial_eta(L["nr"], L["h1"] / Rmin)

    NI = len(zn)          # axial nodes
    NX = NI - 1           # axial cells
    NR = L["nr"]          # radial cells
    th = math.radians(WEDGE_HALF_DEG)
    c, s = math.cos(th), math.sin(th)
    mm = 1e-3

    # points: axis (NI), front (NI*NR), back (NI*NR)
    pts = [(zn[i] * mm, 0.0, 0.0) for i in range(NI)]
    for side in (-1.0, 1.0):
        for i in range(NI):
            for j in range(1, NR + 1):
                r = eta[j] * Rn[i] * mm
                pts.append((zn[i] * mm, r * c, side * r * s))

    def P(i, j, side):          # side 0 = front (-z), 1 = back (+z)
        if j == 0:
            return i
        return NI + side * NI * NR + i * NR + (j - 1)

    def cell(i, j):
        return i * NR + j

    def orient(f, want):
        """Reverse face f unless its normal has positive dot with 'want'."""
        n = newell([pts[k] for k in f])
        if n[0] * want[0] + n[1] * want[1] + n[2] * want[2] < 0:
            f = f[::-1]
        return f

    def axial_face(i, j):       # at node column i, cell row j
        if j == 0:
            return [P(i, 0, 0), P(i, 1, 0), P(i, 1, 1)]
        return [P(i, j, 0), P(i, j + 1, 0), P(i, j + 1, 1), P(i, j, 1)]

    def radial_face(i, j):      # at node row j, cell column i
        return [P(i, j, 0), P(i + 1, j, 0), P(i + 1, j, 1), P(i, j, 1)]

    def wedge_face(i, j, side):
        return [P(i, j, side), P(i + 1, j, side), P(i + 1, j + 1, side), P(i, j + 1, side)]

    faces, owner, neigh = [], [], []
    X, Y, Zm, Zp = (1, 0, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)
    # internal faces, upper-triangular order
    for i in range(NX):
        for j in range(NR):
            co = cell(i, j)
            if j < NR - 1:
                faces.append(orient(radial_face(i, j + 1), Y))
                owner.append(co); neigh.append(cell(i, j + 1))
            if i < NX - 1:
                faces.append(orient(axial_face(i + 1, j), X))
                owner.append(co); neigh.append(cell(i + 1, j))
    n_internal = len(faces)
    patches = []

    def add_patch(name, ptype, flist):
        start = len(faces)
        for f, o in flist:
            faces.append(f); owner.append(o)
        patches.append((name, ptype, len(flist), start))

    add_patch("inlet", "patch",
              [(orient(axial_face(0, j), (-1, 0, 0)), cell(0, j)) for j in range(NR)])
    add_patch("outlet", "patch",
              [(orient(axial_face(NX, j), X), cell(NX - 1, j)) for j in range(NR)])
    add_patch("wall", "wall",
              [(orient(radial_face(i, NR), Y), cell(i, NR - 1)) for i in range(NX)])
    add_patch("front", "wedge",
              [(orient(wedge_face(i, j, 0), Zm), cell(i, j)) for i in range(NX) for j in range(NR)])
    add_patch("back", "wedge",
              [(orient(wedge_face(i, j, 1), Zp), cell(i, j)) for i in range(NX) for j in range(NR)])

    n_cells = NX * NR
    pm = os.path.join(case, "constant", "polyMesh")
    os.makedirs(pm, exist_ok=True)
    note = f"nPoints:{len(pts)}  nCells:{n_cells}  nFaces:{len(faces)}  nInternalFaces:{n_internal}"
    with open(os.path.join(pm, "points"), "w", newline="\n") as fh:
        fh.write(foam_header("vectorField", "points"))
        fh.write(f"{len(pts)}\n(\n")
        fh.write("".join(f"({x:.16g} {y:.16g} {z:.16g})\n" for x, y, z in pts))
        fh.write(")\n")
    with open(os.path.join(pm, "faces"), "w", newline="\n") as fh:
        fh.write(foam_header("faceList", "faces"))
        fh.write(f"{len(faces)}\n(\n")
        fh.write("".join(f"{len(f)}({' '.join(map(str, f))})\n" for f in faces))
        fh.write(")\n")
    for name, data in (("owner", owner), ("neighbour", neigh)):
        with open(os.path.join(pm, name), "w", newline="\n") as fh:
            fh.write(foam_header("labelList", name, note))
            fh.write(f"{len(data)}\n(\n")
            fh.write("\n".join(map(str, data)))
            fh.write("\n)\n")
    with open(os.path.join(pm, "boundary"), "w", newline="\n") as fh:
        fh.write(foam_header("polyBoundaryMesh", "boundary"))
        fh.write(f"{len(patches)}\n(\n")
        for name, ptype, n, start in patches:
            fh.write(f"    {name}\n    {{\n        type            {ptype};\n"
                     f"        nFaces          {n};\n        startFace       {start};\n    }}\n")
        fh.write(")\n")

    # wall-adjacent cell-centre distance (2-D quad centroid to wall segment), m
    y1 = []
    for i in range(NX):
        quad = [(zn[i], eta[NR - 1] * Rn[i]), (zn[i + 1], eta[NR - 1] * Rn[i + 1]),
                (zn[i + 1], Rn[i + 1]), (zn[i], Rn[i])]
        cz = sum(q[0] for q in quad) / 4
        cr = sum(q[1] for q in quad) / 4
        tz, tr = _unit(zn[i + 1] - zn[i], Rn[i + 1] - Rn[i])
        y1.append(abs((cz - zn[i]) * tr - (cr - Rn[i]) * tz) * mm)

    info = dict(geometry=geom, level=level, params=L, axial_growth=AXIAL_GROWTH,
                wedge_half_deg=WEDGE_HALF_DEG, NX=NX, NR=NR, n_cells=n_cells,
                z_nodes_mm=zn, R_nodes_mm=Rn, eta=eta, stations=keys,
                wall_profile_mm=[(p[0], p[1]) for p in wall],
                wall_y1_m=y1, inlet_radius_mm=Rn[0], outlet_radius_mm=Rn[-1])
    with open(os.path.join(case, "mesh_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geom", required=True, help="geometry name in bore_stations.csv")
    ap.add_argument("--case", required=True, help="case directory to write into")
    ap.add_argument("--level", default="medium", choices=sorted(LEVELS))
    ap.add_argument("--csv", default=DEFAULT_CSV)
    a = ap.parse_args()
    info = build_mesh(a.geom, a.case, a.level, a.csv)
    print(f"{a.geom} [{a.level}]: {info['NX']} x {info['NR']} = {info['n_cells']} cells, "
          f"first cell {info['params']['h1']*1000:.1f} um at narrowest section, "
          f"wall y1 {min(info['wall_y1_m'])*1e6:.2f}-{max(info['wall_y1_m'])*1e6:.2f} um")


if __name__ == "__main__":
    main()
