#!/usr/bin/env python3
"""
build_assembly.py -- 3-D water volume (half domain) for a complete assembly:
pump-spout coupler (straight inlet extension) -> reducer -> elbow -> Venturi.

Usage (WSL, gmsh >= 4.11):
    python3 build_assembly.py <upstream> <venturi tag> <out dir> [--vane]

  upstream:  esig34  3/4in E-Sig (flat side down) + 3/4in threaded elbow
             esig1   1in E-Sig (flat side down) + McMaster-Carr 4880K464 1in elbow
             lr1     1in long-radius reducing bend (one part), bend turning down
  venturi:   a tag in geometry/bore_stations.csv (current, A_3-4in, B_3-4in, V1..V4)
  --vane:    cut the outlet cross-vane (blades 1.2 mm thick, 12 mm long; 4 blades,
             or 6 for V4), one blade in the bend (symmetry) plane

Writes fluid.stl (metres, one solid per patch: inlet, outlet, symmetry, wall_esig,
wall_elbow, wall_venturi), fluid_mm.stl and assembly_info.json (stations used by
setup_asm_case.py, post_asm.sh/.py and build_viz3d.py). For lr1 "wall_esig" is the
straight 44 mm inlet part of the reducer and "wall_elbow" its reducing bend.

Frame (mm): the Venturi axis is x with flow in -x (down; gravity -x); the leg toward
the pump is +y with flow in -y; the bend plane is z = 0 (symmetry), half domain z >= 0.

(build_fluid.py built the first assembly, esig34 + current; this script reproduces it
with the same dimensions and generalises it.)
"""
import json, math, os, shutil, sys, tempfile
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
GEOM = os.path.join(REPO, "geometry")
sys.path.insert(0, os.path.join(HERE, ".."))
import make_wedge_mesh as mwm  # noqa: E402

args = [a for a in sys.argv[1:] if not a.startswith("--")]
UP, TAG, out_dir = args
VANE = "--vane" in sys.argv
os.makedirs(out_dir, exist_ok=True)

INLET_EXT = 130.0                       # straight inlet extension at the spigot bore
NPT_ENGAGE = {"3/4": 8.6 + 3 * 25.4 / 14, "1": 10.16 + 3 * 25.4 / 11.5}   # L1 hand-tight + 3 wrench turns
VANE_T, VANE_L = 1.2, 12.0

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", int(os.environ.get("GMSH_VERBOSE", "0")))
gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
occ = gmsh.model.occ
TMP = tempfile.mkdtemp()


def step(rel):
    """Copy a STEP file to a space-free temp path before import (OCC/gmsh are safer that way)."""
    dst = os.path.join(TMP, os.path.basename(rel).replace(" ", "_"))
    shutil.copy(os.path.join(GEOM, rel), dst)
    return dst


assumptions = []
parts = []          # (3, tag) solids, full (both halves); halved and fused at the end

# =================================================================== E-Sig (both sizes)
ES_X0, ES_XJ, ES_XS, ES_XE, ES_W = -43.82, 14.0, 89.8, 197.4, 2.54
ES_SPIGOT_R, ES_SPIGOT_ZC = 21.965, -24.505


def build_esig(step_rel, small_r, end_y, x_first=ES_XJ):
    """Inner water volume of an E-Sig (CAD frame: axis x, flow +x, flat side at z = 0),
    installed flat side DOWN, its small-bore axis on the elbow leg axis (x = 0, z = 0),
    its end face at elbow y = end_y. Returns (parts, inlet y, inlet centre x, spigot entry y)."""
    small_zc = -ES_W - small_r

    def to_elbow(xE, zE):          # flat side down: turned 180 deg about the small-bore axis
        return -(zE - small_zc), end_y + (ES_XE - xE)

    vE = occ.importShapes(step(step_rel))
    occ.synchronize()

    def outer_bottom(X):
        r = occ.addRectangle(-200, -200, 0, 400, 400)
        occ.rotate([(2, r)], 0, 0, 0, 0, 1, 0, -math.pi / 2)
        occ.translate([(2, r)], X, 0, 0)
        o, _ = occ.intersect(vE, [(2, r)], removeObject=False, removeTool=True)
        occ.synchronize()
        zb = min(gmsh.model.getBoundingBox(2, t)[2] for d, t in o)
        occ.remove(o, recursive=True)
        occ.synchronize()
        return zb

    xs = [x_first + i * (ES_XS - x_first) / 25 for i in range(26)]
    inner_bot = [(0.0, -46.47), (4.0, -46.45), (8.0, -46.32), (11.0, -46.1)]   # spigot/collar (same in both sizes)
    inner_bot += [(x, outer_bottom(x) + ES_W) for x in xs]
    occ.remove(vE, recursive=True)
    occ.synchronize()
    loops = []
    for xE, zb in inner_bot:
        top = -ES_W
        r = 0.5 * (top - zb)
        c = occ.addCircle(0, 0, 0, r)
        occ.rotate([(1, c)], 0, 0, 0, 0, 0, 1, math.pi / 2)      # seam off the symmetry plane
        occ.rotate([(1, c)], 0, 0, 0, 1, 0, 0, -math.pi / 2)     # normal -> +y
        X, Y = to_elbow(xE, top - r)
        occ.translate([(1, c)], X, Y, 0)
        loops.append(occ.addCurveLoop([c]))
    loft = occ.addThruSections(loops, makeSolid=True, makeRuled=False)
    X0, Y0 = to_elbow(ES_X0 - INLET_EXT, ES_SPIGOT_ZC)
    spig = occ.addCylinder(X0, Y0, 0, 0, -(-(ES_X0 - INLET_EXT)), 0, ES_SPIGOT_R)
    Xs, Ys = to_elbow(ES_XS - 0.5, small_zc)
    small = occ.addCylinder(Xs, Ys, 0, 0, -(ES_XE - ES_XS + 0.5), 0, small_r)
    return [(3, spig), (3, small)] + loft, Y0, X0, to_elbow(ES_X0, 0)[1]


# =================================================================== elbows
def build_elbow(shoulder, bore_r, bend_r, minor_r, esig_end_y, vent_in_x, vent_threaded):
    """Elbow water volume: female-thread cavities (minor radius) from the E-Sig end face and
    from the Venturi inlet face to each socket shoulder, plus the quarter-torus bore.
    Bend centre at (-bend_r, bend_r); the bore leaves the bend tangent at x = 0 / y = 0."""
    p = []
    p.append(occ.addCylinder(0, esig_end_y, 0, 0, shoulder - esig_end_y, 0, minor_r))
    tor = occ.addTorus(-bend_r, bend_r, 0, bend_r, bore_r, angle=math.pi / 2)
    occ.rotate([(3, tor)], -bend_r, bend_r, 0, 0, 0, 1, -math.pi / 2)
    p.append(tor)
    if shoulder > bend_r:          # straight bore between the shoulders and the bend
        p.append(occ.addCylinder(0, shoulder + 0.5, 0, 0, bend_r - shoulder - 0.5, 0, bore_r))
        p.append(occ.addCylinder(-shoulder - 0.5, 0, 0, -(bend_r - shoulder - 0.5), 0, 0, bore_r))
    else:
        p.append(occ.addCylinder(0, bend_r + 0.5, 0, 0, -0.5, 0, bore_r))
        p.append(occ.addCylinder(-bend_r + 0.5, 0, 0, -0.5, 0, 0, bore_r))
    if vent_threaded:                  # threaded Venturi: cavity below the shoulder down to its inlet face
        p.append(occ.addCylinder(-shoulder, 0, 0, vent_in_x - (-shoulder), 0, 0, minor_r))
    elif vent_in_x < -bend_r + 0.5:    # socket below the bend: leg bore from the bend down to the Venturi face
        p.append(occ.addCylinder(-bend_r + 0.5, 0, 0, vent_in_x - (-bend_r + 0.5), 0, 0, bore_r))
    # else: the Venturi face (socket shoulder) is inside the bend; the torus alone reaches it
    return [(3, t) for t in p]


# =================================================================== Venturi
st = mwm.load_stations(os.path.join(GEOM, "bore_stations.csv"), TAG)
wall = mwm.build_wall(st, 0.3)
ks = mwm.key_stations(st)
three_q = TAG in ("current",) or TAG.endswith("3-4in")
if UP == "esig34" and not three_q:
    raise SystemExit(f"{TAG} has a 1in inlet; esig34 is the 3/4in train")
if UP in ("esig1", "lr1") and three_q:
    raise SystemExit(f"{TAG} has a 3/4in inlet; {UP} is a 1in train")

if UP == "esig34":
    # 3/4in threaded elbow (from its STEP, inch units): female leg face y 35.56, shoulder
    # 14.73, bore 20.93, bend radius 14.73; Venturi socket shoulder at x = -15.49.
    FACE, SHOULDER, BORE_R, BEND_R, MINOR_R = 35.56, 14.73, 10.465, 14.732, 11.9
    ENG = 14.0
    ES_END_Y = FACE - ENG
    VENT_IN_X = -15.49
    es, Y_IN, X_IN, SPIG_Y = build_esig("upstream/ESig_reducer_3-4in_NPT_male.step", 10.46, ES_END_Y)
    parts += es
    parts += build_elbow(SHOULDER, BORE_R, BEND_R, MINOR_R, ES_END_Y, VENT_IN_X, False)
    LIP_R = BORE_R
    FINE = [None, -16.0, -1.0, 14.0, ES_END_Y + 12.0, 25.0]
    INSIDE = (0.0, 100.0, 5.0)
    assumptions += [
        f"E-Sig 3/4in NPT male engaged {ENG:.1f} mm into the elbow female leg (L1 8.6 mm + 3 wrench turns)",
        f"{TAG} Venturi spigot fully seated in the elbow socket; inlet face on the socket shoulder (same as the current part)",
        "Female-thread cavity between E-Sig end and elbow shoulder modelled as a 23.8 mm cylinder",
        "E-Sig flat side down (spigot end high, S-curve drops to the elbow)"]
    WALL_SPLIT_Y = ES_END_Y
elif UP == "esig1":
    # McMaster-Carr 4880K464 (from its STEP, inch units): 1in NPT female both legs, faces
    # 38.89 mm from the leg-axis intersection, socket shoulders at 15.88, bore 26.4,
    # bend centreline radius 19.81, thread minor (crest) radius about 14.65.
    FACE, SHOULDER, BORE_R, BEND_R, MINOR_R = 38.89, 15.88, 13.2, 19.81, 14.65
    ENG = NPT_ENGAGE["1"]
    ES_END_Y = FACE - ENG
    # V1-V4 have a plain 33.4 mm (1in pipe OD) spigot, a press fit: the Venturi leg is treated
    # as a 1in slip socket with the Venturi seated on the socket shoulder (a slip x FNPT elbow
    # with the McMaster bore and bend). The 26.4 mm bore meets the 31.2 mm mouth: a real step.
    VENT_IN_X = -SHOULDER
    es, Y_IN, X_IN, SPIG_Y = build_esig("upstream/ESig_reducer_1in_NPT_male.step", 13.32, ES_END_Y, x_first=ES_XJ + 0.5)   # x = 14 is the collar face on the 1in part
    parts += es
    parts += build_elbow(SHOULDER, BORE_R, BEND_R, MINOR_R, ES_END_Y, VENT_IN_X, False)
    LIP_R = 12.5
    FINE = [None, -18.0, -1.0, 16.0, ES_END_Y + 12.0, 25.0]
    INSIDE = (0.0, 100.0, 5.0)
    assumptions += [
        f"1in E-Sig NPT male engaged {ENG:.1f} mm into the McMaster 4880K464 elbow (L1 10.2 mm + 3 wrench turns); thread cavity to the shoulder modelled as a 29.3 mm cylinder",
        f"{TAG} plain 33.4 mm spigot press-fitted and seated on the Venturi-leg socket shoulder (leg treated as a 1in slip socket); 26.4 mm elbow bore steps to the {2 * wall[0][1]:.1f} mm Venturi mouth",
        "1in E-Sig flat side down (spigot end high, S-curve drops to the elbow)"]
    WALL_SPLIT_Y = ES_END_Y
elif UP == "lr1":
    # Long-radius reducing bend (one part; CAD frame: spigot along +x, bend turning to +z,
    # symmetric about y = 868.7). OCC cannot run booleans on this STEP solid, but plane
    # sections work, so the bore is rebuilt from sections: a straight 43.95 mm bore, then
    # N_SEC sections in a fan of planes about the line (x 614.43, z 1334.36) // y, each fitted
    # with a circle (centre and radius in its plane) and lofted; then the bore at the port
    # shoulder and the gap up to the Venturi face. Installed with the port DOWN.
    import numpy as np
    PX, PY = 690.52, 868.70                 # port axis (x, y)
    PORT_FACE, PORT_R, SHOULDER_Z, SHOULDER_R = 1352.36, 16.90, 1334.36, 15.58
    SPIG_END, SPIG_R, SPIG_ZC = 526.16, 21.975, 1262.29
    BEND_START = 614.43
    C = np.array([BEND_START, PY, SHOULDER_Z])
    # The Venturi's plain 33.4 mm spigot is a press fit in the plain 33.8 mm port and is
    # pushed in to the port shoulder (spigot 18 mm long, port 18 mm deep): its 31.2 mm mouth
    # meets the reducer's 31.16 mm bore end flush.
    VZ = SHOULDER_Z                         # Venturi inlet face (CAD z)
    ENG = PORT_FACE - SHOULDER_Z
    vL = occ.importShapes(step("upstream/LongRadius_reducer_1in_NPT_female.step"))
    occ.synchronize()
    ey = np.array([0.0, 1.0, 0.0])

    def frame_matrix(a, b, c, o):
        """affineTransform list mapping local (x, y, z) -> o + x a + y b + z c."""
        return [a[0], b[0], c[0], o[0], a[1], b[1], c[1], o[1], a[2], b[2], c[2], o[2]]

    def fit_section(th):
        d = np.array([math.sin(th), 0.0, -math.cos(th)])       # from C toward the bore
        n = np.cross(d, ey)                                      # plane normal (flow direction)
        r = occ.addRectangle(20.0, -60.0, 0.0, 120.0, 120.0)     # local u along d, v along y
        occ.affineTransform([(2, r)], frame_matrix(d, ey, n, C))
        o, _ = occ.intersect(vL, [(2, r)], removeObject=False, removeTool=True)
        occ.synchronize()
        P = []
        for _, cv in gmsh.model.getBoundary(o, oriented=False, combined=False):
            cv = abs(cv)
            lo, hi = gmsh.model.getParametrizationBounds(1, cv)
            P.append(np.array(gmsh.model.getValue(1, cv, list(np.linspace(lo[0], hi[0], 40)))).reshape(-1, 3))
        occ.remove(o, recursive=True)
        occ.synchronize()
        P = np.vstack(P) - C
        uv = np.c_[P @ d, P @ ey]
        rr = np.linalg.norm(uv - uv.mean(0), axis=1)
        inner = uv[rr < 0.5 * (rr.min() + rr.max())]               # bore loop, not the outer wall
        A = np.c_[2 * inner, np.ones(len(inner))]
        cu, cv_, k = np.linalg.lstsq(A, (inner ** 2).sum(1), rcond=None)[0]
        rad = math.sqrt(k + cu ** 2 + cv_ ** 2)
        return C + cu * d + cv_ * ey, n, rad

    N_SEC = 31
    secs = [fit_section(math.radians(89.0 * i / (N_SEC - 1))) for i in range(N_SEC)]
    occ.remove(vL, recursive=True)
    occ.synchronize()
    loops = []
    for c0, n, rad in secs:
        a = ey                                       # circle seam off the bend (symmetry) plane
        b2 = np.cross(n, a)
        cc = occ.addCircle(0, 0, 0, rad)
        occ.affineTransform([(1, cc)], frame_matrix(a, b2, n, c0))
        loops.append(occ.addCurveLoop([cc]))
    lr_parts = occ.addThruSections(loops, makeSolid=True, makeRuled=False)
    lr_parts.append((3, occ.addCylinder(SPIG_END + 0.2, PY, SPIG_ZC, BEND_START + 0.3 - SPIG_END - 0.2, 0, 0, SPIG_R)))
    zl = secs[-1][0][2] - 1.0
    lr_parts.append((3, occ.addCylinder(PX, PY, zl, 0, 0, SHOULDER_Z - zl, SHOULDER_R)))
    if os.environ.get("DBG"):
        for c0, n, rad in secs[::5]:
            print("section", [round(v, 2) for v in c0], round(rad, 2))
    # CAD -> assembly frame: x' = -(z - VZ), y' = -(x - PX), z' = y - PY  (proper rotation)
    occ.affineTransform(lr_parts, [0, 0, -1, VZ, -1, 0, 0, PX, 0, 1, 0, -PY])
    parts += lr_parts
    SPIG_Y = PX - (SPIG_END + 0.2)
    X_IN = VZ - SPIG_ZC
    Y_IN = SPIG_Y + INLET_EXT
    parts.append((3, occ.addCylinder(X_IN, Y_IN, 0, 0, -INLET_EXT - 0.5, 0, SPIG_R)))
    ES_END_Y = PX - BEND_START              # "reducer straight part" / "bend" split
    VENT_IN_X = 0.0
    LIP_R = min(PORT_R, wall[0][1]) - 0.2
    # level-2 box over the Venturi and the last 25 mm of the bend above it (as for the elbows);
    # a box over the whole reducing bend gave ~3x the cells of the other assemblies
    FINE = [None, -20.0, -1.0, 25.0, 40.0, 25.0]
    INSIDE = (X_IN, SPIG_Y - 20.0, 5.0)
    WALL_SPLIT_Y = ES_END_Y
    assumptions += [
        "Long-radius reducing bend installed with the 1in port DOWN (spigot horizontal toward the pump)",
        f"{TAG} plain spigot press-fitted into the reducer's 1in port and seated on the port shoulder ({ENG:.1f} mm deep); bore steps 31.16 -> {2 * wall[0][1]:.2f} mm, flush",
        f"Reducer bore rebuilt from {N_SEC} circle-fitted sections of its STEP (blank port, thread not modelled)"]
else:
    raise SystemExit("upstream must be esig34, esig1 or lr1")

assumptions += [f"Inlet extension of {INLET_EXT:.0f} mm at the spigot bore stands in for the rubber coupler and pump spout",
                "Half domain: symmetric about the bend plane (the chlorine port is not modelled)"]

# Venturi bore, revolved about the x axis below VENT_IN_X
pts = [occ.addPoint(VENT_IN_X - z, r, 0) for z, r, _ in wall]
z_out = VENT_IN_X - wall[-1][0]
pa = [occ.addPoint(z_out, 0, 0), occ.addPoint(VENT_IN_X, 0, 0)]
allp = pts + pa
lines = [occ.addLine(allp[i], allp[(i + 1) % len(allp)]) for i in range(len(allp))]
prof = occ.addPlaneSurface([occ.addCurveLoop(lines)])
vent = [dt for dt in occ.revolve([(2, prof)], 0, 0, 0, 1, 0, 0, 2 * math.pi) if dt[0] == 3]
lip = occ.addCylinder(VENT_IN_X + 0.3, 0, 0, -0.6, 0, 0, LIP_R)
if VANE:
    nbl = 6 if TAG == "V4" else 4
    dex = 2 * wall[-1][1]
    cutters = []
    for k in range(nbl // 2):
        bx = occ.addBox(z_out - 1.0, -(dex + 2) / 2, -VANE_T / 2, VANE_L + 1.0, dex + 2, VANE_T)
        occ.rotate([(3, bx)], 0, 0, 0, 1, 0, 0, k * 2 * math.pi / nbl)
        cutters.append((3, bx))
    vent, _ = occ.cut(vent, cutters)
    assumptions.append(f"Outlet cross-vane: {nbl} blades, {VANE_T} mm thick, {VANE_L:.0f} mm long, one blade in the bend plane")
parts += vent + [(3, lip)]

# ------------------------------------------------------------------ halve and fuse
ZCUT = 0.0
halves = []
for dt in parts:
    bx = occ.addBox(-2000, -2000, -2000, 4000, 4000, 2000 + ZCUT)
    m0 = occ.getMass(*dt)
    h, _ = occ.cut([dt], [(3, bx)])
    if os.environ.get("DBG"):
        print("half", dt, round(m0), "->", [round(occ.getMass(d, t)) for d, t in h])
    halves += h
fluid, _ = occ.fuse(halves[:1], halves[1:])
occ.synchronize()
vols = gmsh.model.getEntities(3)
vol_mm3 = [round(occ.getMass(3, t)) for _, t in vols]
print("volumes", vols, "volume mm3", vol_mm3)
if len(vols) != 1:
    raise SystemExit("water volume is not a single connected solid")
bb = gmsh.model.getBoundingBox(3, vols[0][1])

# ------------------------------------------------------------------ patches
patches = {"inlet": [], "outlet": [], "symmetry": [], "wall_esig": [], "wall_elbow": [], "wall_venturi": []}
for _, s in gmsh.model.getBoundary(vols, oriented=False):
    s = abs(s)
    b = gmsh.model.getBoundingBox(2, s)
    cx, cy, cz = occ.getCenterOfMass(2, s)
    if abs(b[2] - ZCUT) < 1e-3 and abs(b[5] - ZCUT) < 1e-3:
        patches["symmetry"].append(s)
    elif abs(b[1] - Y_IN) < 1e-3 and abs(b[4] - Y_IN) < 1e-3:
        patches["inlet"].append(s)
    elif abs(b[0] - z_out) < 1e-3 and abs(b[3] - z_out) < 1e-3:
        patches["outlet"].append(s)
    elif cx < VENT_IN_X - 0.05:
        patches["wall_venturi"].append(s)
    elif cy > WALL_SPLIT_Y + 0.05:
        patches["wall_esig"].append(s)
    else:
        patches["wall_elbow"].append(s)
print({k: len(v) for k, v in patches.items()})
for k, v in patches.items():
    if not v:
        raise SystemExit(f"patch {k} is empty")
    gmsh.model.addPhysicalGroup(2, v, name=k)

gmsh.option.setNumber("Mesh.Algorithm", 1)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.write(os.path.join(out_dir, "fluid.stl"))
gmsh.option.setNumber("Mesh.ScalingFactor", 1.0)
gmsh.write(os.path.join(out_dir, "fluid_mm.stl"))

FINE[0] = z_out - 1.0
info = dict(
    frame="mm; Venturi axis = x, flow -x (down, gravity -x); leg toward the pump = +y, flow -y; z>=0 half",
    upstream=UP, venturi=TAG, vane=VANE, engage_mm=ENG,
    esig_end_face_y=ES_END_Y, venturi_inlet_x=VENT_IN_X, outlet_x=z_out,
    inlet_y=Y_IN, inlet_centre_x=X_IN, inlet_radius=ES_SPIGOT_R if UP != "lr1" else SPIG_R,
    inlet_extension_mm=INLET_EXT, spigot_entry_y=SPIG_Y,
    port_x=VENT_IN_X - ks["port"], port_side="not modelled; wall pressure sampled all round",
    diffuser_start_mm=ks["diffuser_start"], diffuser_end_mm=ks["diffuser_end"],
    venturi_wall=[(z, r) for z, r, _ in wall], fluid_volume_mm3_half=vol_mm3[0],
    bbox_mm=[bb[0], bb[3], bb[1], bb[4], bb[5]], inside_point_mm=list(INSIDE), fine_box_mm=FINE,
    assumptions=assumptions)
json.dump(info, open(os.path.join(out_dir, "assembly_info.json"), "w"), indent=1)
gmsh.finalize()
shutil.rmtree(TMP, ignore_errors=True)
