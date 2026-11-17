#!/usr/bin/env python3
"""
make_vof_mesh.py -- 2-D axisymmetric wedge mesh of a Venturi bore PLUS an air
region below the outlet, for Stage 4 two-phase (VoF) runs.

The bore part is identical to make_wedge_mesh.py (same stations, fillets,
axial/radial distributions, wall-resolved first cell), so cell numbering of
the bore is unchanged: bore cell (i, j) = i*NR + j for i < NXd.

Below the outlet face an air box is added:
  * core  : continuation of the bore cells at the outlet radius R_out
  * annulus: R_out .. BOX_R_FACTOR*R_out, NA cells, first cell matched to the
    core's outermost cell in every column (geometric growth)
  * length BOX_L_FACTOR*D_out, axial cells growing from the last bore cell up
    to PLUME_H_FACTOR*h_core
  * the bore's wall clustering is relaxed over the first BLEND_FACTOR*D_out of
    the box so the air box does not inherit 3-um radial cells everywhere.

Patches: inlet (patch), wall (wall, bore only), atmosphere (patch: box side,
box bottom and the annular face level with the outlet), front/back (wedge).
cellZone "bore" = all bore cells (for water-volume monitoring).

Usage
  python3 make_vof_mesh.py --geom current --case "<case dir>" [--level medium]
"""
import argparse
import json
import math
import os

import make_wedge_mesh as mwm

BOX_R_FACTOR = 3.0      # box radius / outlet radius
BOX_L_FACTOR = 4.0      # box length / outlet diameter
BLEND_FACTOR = 0.6      # radial-clustering relaxation length / outlet diameter
PLUME_H_FACTOR = 2.5    # max axial cell in the box / h_core
PLUME_GROWTH = 1.05
NA_CELLS = 60           # annulus radial cells
SOFT_FIRST = 0.5        # relaxed core edge cell / h_core


def geometric_sizes(a, L, n):
    """n sizes a*r^k summing to L (r found by bisection, r may be < 1)."""
    def tot(r):
        return a * n if abs(r - 1) < 1e-12 else a * (r ** n - 1) / (r - 1)
    lo, hi = 0.5, 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tot(mid) < L:
            lo = mid
        else:
            hi = mid
    r = 0.5 * (lo + hi)
    s = [a * r ** k for k in range(n)]
    f = L / sum(s)
    return [x * f for x in s], r


# Stage 4 levels. The Stage 1 levels put the first cell at 2-4.5 um (sized for
# y+ ~ 1 at 30 L/min); in VoF an interface in such cells drives parasitic
# air velocities of ~1 m/s and time steps of ~3e-6 s. At 4-10 L/min a 25 um
# first cell still gives y+ <~ 2 (Stage 1: y+max 0.30 at 6 L/min with 3 um),
# so the wall stays resolved without wall functions.
VOF_LEVELS = {
    "vof":      dict(nr=60, h_core=0.250, h_corner=0.050, h1=25e-3),
    "vof_fine": dict(nr=90, h_core=0.1667, h_corner=0.0333, h1=16.7e-3),
}


def build_mesh(geom, case, level="vof", csv_path=mwm.DEFAULT_CSV):
    L = VOF_LEVELS[level] if level in VOF_LEVELS else mwm.LEVELS[level]
    st = mwm.load_stations(csv_path, geom)
    wall = mwm.build_wall(st, L["h_core"])
    keys = mwm.key_stations(st)
    zd = mwm.axial_nodes(wall, L["h_core"], L["h_corner"])
    Rd = [mwm.wall_radius(wall, z) for z in zd]
    NR = L["nr"]
    eta_w = mwm.radial_eta(NR, L["h1"] / min(Rd))
    Rout = Rd[-1]
    Dout = 2 * Rout
    eta_s = mwm.radial_eta(NR, SOFT_FIRST * L["h_core"] / Rout)

    # ---- box axial nodes
    zp = [zd[-1]]
    h = zd[-1] - zd[-2]
    hmax = PLUME_H_FACTOR * L["h_core"]
    zend = zd[-1] + BOX_L_FACTOR * Dout
    while zp[-1] + h < zend - 0.3 * h:
        zp.append(zp[-1] + h)
        h = min(hmax, h * PLUME_GROWTH)
    zp.append(zend)
    zn = zd + zp[1:]
    NXd = len(zd) - 1
    NI = len(zn)
    NX = NI - 1
    NA = NA_CELLS
    Lr = (BOX_R_FACTOR - 1.0) * Rout

    # ---- radial node positions per column (mm); core j=0..NR, annulus k=1..NA
    rcore, rann = [], []
    for i, z in enumerate(zn):
        if i <= NXd:
            eta = eta_w
            R = Rd[i]
        else:
            w = min(1.0, (z - zd[-1]) / (BLEND_FACTOR * Dout))
            w = 0.5 - 0.5 * math.cos(math.pi * w)
            eta = [(1 - w) * a + w * b for a, b in zip(eta_w, eta_s)]
            R = Rout
        rcore.append([e * R for e in eta])
        if i >= NXd:
            a = R - rcore[-1][NR - 1]
            sz, _ = geometric_sizes(a, Lr, NA)
            r, acc = [], Rout
            for x in sz:
                acc += x
                r.append(acc)
            r[-1] = BOX_R_FACTOR * Rout
            rann.append(r)
        else:
            rann.append(None)

    th = math.radians(mwm.WEDGE_HALF_DEG)
    c, s = math.cos(th), math.sin(th)
    mm = 1e-3

    # ---- points: axis, then (i, j>=1, side) for existing nodes
    pid = {}
    pts = []
    for i in range(NI):
        pid[(i, 0, 0)] = pid[(i, 0, 1)] = len(pts)
        pts.append((zn[i] * mm, 0.0, 0.0))
    for side in (0, 1):
        sg = -1.0 if side == 0 else 1.0
        for i in range(NI):
            jmax = NR + (NA if i >= NXd else 0)
            for j in range(1, jmax + 1):
                r = (rcore[i][j] if j <= NR else rann[i][j - NR - 1]) * mm
                pid[(i, j, side)] = len(pts)
                pts.append((zn[i] * mm, r * c, sg * r * s))

    def P(i, j, side):
        return pid[(i, j, side)]

    # ---- cells
    NJ = NR + NA
    cid = {}
    for i in range(NX):
        for j in range(NJ):
            if j < NR or i >= NXd:
                cid[(i, j)] = len(cid)
    n_cells = len(cid)

    def active(i, j):
        return (i, j) in cid

    def orient(f, want):
        n = mwm.newell([pts[k] for k in f])
        return f if n[0] * want[0] + n[1] * want[1] + n[2] * want[2] >= 0 else f[::-1]

    def axial_face(i, j):
        if j == 0:
            return [P(i, 0, 0), P(i, 1, 0), P(i, 1, 1)]
        return [P(i, j, 0), P(i, j + 1, 0), P(i, j + 1, 1), P(i, j, 1)]

    def radial_face(i, j):
        return [P(i, j, 0), P(i + 1, j, 0), P(i + 1, j, 1), P(i, j, 1)]

    def wedge_face(i, j, side):
        return [P(i, j, side), P(i + 1, j, side), P(i + 1, j + 1, side), P(i, j + 1, side)]

    X, Y = (1, 0, 0), (0, 1, 0)
    faces, owner, neigh = [], [], []
    for i in range(NX):
        for j in range(NJ):
            if not active(i, j):
                continue
            co = cid[(i, j)]
            if active(i, j + 1):
                faces.append(orient(radial_face(i, j + 1), Y)); owner.append(co); neigh.append(cid[(i, j + 1)])
            if i + 1 < NX and active(i + 1, j):
                faces.append(orient(axial_face(i + 1, j), X)); owner.append(co); neigh.append(cid[(i + 1, j)])
    n_internal = len(faces)

    inlet = [(orient(axial_face(0, j), (-1, 0, 0)), cid[(0, j)]) for j in range(NR)]
    wallp = [(orient(radial_face(i, NR), Y), cid[(i, NR - 1)]) for i in range(NXd)]
    atm = [(orient(axial_face(NXd, j), (-1, 0, 0)), cid[(NXd, j)]) for j in range(NR, NJ)]   # annular top
    atm += [(orient(radial_face(i, NJ), Y), cid[(i, NJ - 1)]) for i in range(NXd, NX)]       # side
    atm += [(orient(axial_face(NX, j), X), cid[(NX - 1, j)]) for j in range(NJ)]             # bottom
    order = sorted(cid, key=lambda k: cid[k])
    front = [(orient(wedge_face(i, j, 0), (0, 0, -1)), cid[(i, j)]) for i, j in order]
    back = [(orient(wedge_face(i, j, 1), (0, 0, 1)), cid[(i, j)]) for i, j in order]

    patches = []
    for name, ptype, fl in (("inlet", "patch", inlet), ("wall", "wall", wallp),
                            ("atmosphere", "patch", atm), ("front", "wedge", front),
                            ("back", "wedge", back)):
        patches.append((name, ptype, len(fl), len(faces)))
        for f, o in fl:
            faces.append(f); owner.append(o)

    pm = os.path.join(case, "constant", "polyMesh")
    os.makedirs(pm, exist_ok=True)
    note = f"nPoints:{len(pts)}  nCells:{n_cells}  nFaces:{len(faces)}  nInternalFaces:{n_internal}"
    with open(os.path.join(pm, "points"), "w", newline="\n") as fh:
        fh.write(mwm.foam_header("vectorField", "points"))
        fh.write(f"{len(pts)}\n(\n" + "".join(f"({x:.16g} {y:.16g} {z:.16g})\n" for x, y, z in pts) + ")\n")
    with open(os.path.join(pm, "faces"), "w", newline="\n") as fh:
        fh.write(mwm.foam_header("faceList", "faces"))
        fh.write(f"{len(faces)}\n(\n" + "".join(f"{len(f)}({' '.join(map(str, f))})\n" for f in faces) + ")\n")
    for name, data in (("owner", owner), ("neighbour", neigh)):
        with open(os.path.join(pm, name), "w", newline="\n") as fh:
            fh.write(mwm.foam_header("labelList", name, note))
            fh.write(f"{len(data)}\n(\n" + "\n".join(map(str, data)) + "\n)\n")
    with open(os.path.join(pm, "boundary"), "w", newline="\n") as fh:
        fh.write(mwm.foam_header("polyBoundaryMesh", "boundary"))
        fh.write(f"{len(patches)}\n(\n")
        for name, ptype, n, start in patches:
            fh.write(f"    {name}\n    {{\n        type            {ptype};\n"
                     f"        nFaces          {n};\n        startFace       {start};\n    }}\n")
        fh.write(")\n")
    n_bore = NXd * NR
    with open(os.path.join(pm, "cellZones"), "w", newline="\n") as fh:
        fh.write(mwm.foam_header("regIOobject", "cellZones"))
        fh.write(f"1\n(\nbore\n{{\n    type cellZone;\ncellLabels      List<label> {n_bore}\n(\n")
        fh.write("\n".join(str(k) for k in range(n_bore)))
        fh.write("\n)\n;\n}\n)\n")

    y1 = []
    for i in range(NXd):
        quad = [(zd[i], eta_w[NR - 1] * Rd[i]), (zd[i + 1], eta_w[NR - 1] * Rd[i + 1]),
                (zd[i + 1], Rd[i + 1]), (zd[i], Rd[i])]
        cz = sum(q[0] for q in quad) / 4
        cr = sum(q[1] for q in quad) / 4
        tz, tr = mwm._unit(zd[i + 1] - zd[i], Rd[i + 1] - Rd[i])
        y1.append(abs((cz - zd[i]) * tr - (cr - Rd[i]) * tz) * mm)

    info = dict(geometry=geom, level=level, params=L, kind="vof",
                wedge_half_deg=mwm.WEDGE_HALF_DEG, NX=NX, NXd=NXd, NR=NR, NA=NA,
                n_cells=n_cells, n_bore_cells=n_bore,
                z_nodes_mm=zn, R_nodes_mm=Rd + [Rout] * (NI - len(Rd)), eta=eta_w,
                stations=keys, wall_profile_mm=[(p[0], p[1]) for p in wall],
                wall_y1_m=y1, inlet_radius_mm=Rd[0], outlet_radius_mm=Rout,
                box_radius_mm=BOX_R_FACTOR * Rout, box_end_mm=zend)
    with open(os.path.join(case, "mesh_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geom", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--level", default="vof")
    ap.add_argument("--csv", default=mwm.DEFAULT_CSV)
    a = ap.parse_args()
    info = build_mesh(a.geom, a.case, a.level, a.csv)
    print(f"{a.geom} [{a.level}] VoF mesh: {info['n_cells']} cells "
          f"({info['n_bore_cells']} bore), box R {info['box_radius_mm']:.1f} mm to z {info['box_end_mm']:.1f} mm")


if __name__ == "__main__":
    main()
