"""Minimal readers for OpenFOAM legacy ASCII VTK surface/line output."""
import glob, os
import numpy as np


def latest_dir(pp, func):
    d = sorted(glob.glob(os.path.join(pp, func, "*")), key=lambda p: float(os.path.basename(p)))
    return d[-1]


def read_vtk(fn):
    """Legacy ASCII VTK polydata -> points (n,3), cells (list of index lists:
    polygons or lines), cell-data dict, point-data dict. Token-stream parser, so
    line wrapping does not matter."""
    tok = open(fn).read().split()
    n = len(tok)
    i = 0
    pts, cells, cdata, pdata, cur = None, [], {}, {}, None

    def take(k):
        nonlocal i
        a = tok[i:i + k]; i += k
        return a

    while i < n:
        t = tok[i]
        if t == "POINTS":
            npts = int(tok[i + 1]); i += 3
            pts = np.array(take(3 * npts), float).reshape(-1, 3)
        elif t in ("POLYGONS", "LINES", "VERTICES", "TRIANGLE_STRIPS"):
            nc, size = int(tok[i + 1]), int(tok[i + 2]); i += 3
            flat = np.array(take(size), int)
            j = 0
            for _ in range(nc):
                m = flat[j]; cells.append(flat[j + 1:j + 1 + m].tolist()); j += 1 + m
        elif t == "CELL_DATA":
            cur = cdata; i += 2
        elif t == "POINT_DATA":
            cur = pdata; i += 2
        elif t == "FIELD":
            nf = int(tok[i + 2]); i += 3
            for _ in range(nf):
                name, nc, nt = tok[i], int(tok[i + 1]), int(tok[i + 2]); i += 4   # name nComp nTuples type
                a = np.array(take(nc * nt), float)
                cur[name] = a.reshape(nt, nc) if nc > 1 else a
        elif t in ("SCALARS", "VECTORS"):
            name = tok[i + 1]; isvec = t == "VECTORS"
            i += 3
            if t == "SCALARS" and tok[i].isdigit():
                i += 1
            if tok[i] == "LOOKUP_TABLE":
                i += 2
            cnt = (len(pts) if cur is pdata else len(cells)) * (3 if isvec else 1)
            a = np.array(take(cnt), float)
            cur[name] = a.reshape(-1, 3) if isvec else a
        else:
            i += 1
    return pts, cells, cdata, pdata


def face_geom(pts, polys):
    """Area vectors and centres of polygons (triangle fan about the centroid)."""
    A = np.zeros((len(polys), 3)); C = np.zeros((len(polys), 3))
    for k, f in enumerate(polys):
        P = pts[f]; c = P.mean(0); a = np.zeros(3)
        for j in range(len(f)):
            a += 0.5 * np.cross(P[j] - c, P[(j + 1) % len(f)] - c)
        A[k] = a; C[k] = c
    return A, C
