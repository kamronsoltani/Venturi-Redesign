#!/usr/bin/env python3
"""
build_fluid.py -- 3-D water volume for: 3/4in E-Sig reducer -> 3/4in threaded
elbow -> current Venturi (group 1, current part). Half domain (z >= 0, the bend
plane is the symmetry plane). Writes STL patches (metres) for snappyHexMesh and
a JSON with key stations.

Frame (mm, "elbow frame"): the Venturi axis is the x axis with flow in -x
(downward; gravity is -x). The threaded elbow leg is the y axis, flow in -y,
coming from the pump. +x is up.

Assembly assumptions (see assembly_info.json):
  * E-Sig 3/4in NPT male end engaged ENGAGE_34 mm into the elbow's female leg
    (hand-tight L1 = 8.6 mm + 3 wrench turns = 5.4 mm).
  * The current Venturi's 26.7 mm OD inlet spigot is fully seated in the elbow's
    26.6 mm socket, so its inlet face sits on the socket shoulder.
  * The female-thread cavity between the E-Sig end face and the elbow shoulder is
    modelled as a cylinder at the thread minor diameter (threads not resolved).
  * E-Sig installed with its flat side DOWN: the 44 mm spigot end sits high and
    the eccentric S-curve drops to the small pipe (ESIG_FLAT_SIDE).

Usage (WSL, gmsh >= 4.11):  python3 build_fluid.py <out dir> <E-Sig 3/4 STEP>
"""
import json, math, os, sys
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import make_wedge_mesh as mwm  # noqa: E402

out_dir, esig_step = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# ---------------------------------------------------------------- dimensions (mm)
ENGAGE_34 = 14.0
# 3/4in threaded elbow (from its STEP, inch units converted)
ELB_FACE_Y = 35.56          # female leg entrance face
ELB_SHOULDER_Y = 14.73      # end of female thread / start of bend
ELB_BORE_R = 10.465         # 20.93 mm bore through the bend
ELB_BEND_R = 14.732         # bend centreline radius
ELB_SOCKET_SHOULDER_X = -15.49
THREAD_MINOR_R = 11.9       # 3/4in NPT female minor radius over the cavity (approx.)
# E-Sig 3/4in (E-Sig frame: axis x, flow +x, flat top at z = 0)
ES_X0, ES_XJ, ES_XS, ES_XE = -43.82, 14.0, 89.8, 197.4
ES_SPIGOT_R, ES_SPIGOT_ZC = 21.965, -24.505
ES_W = 2.54
ES_SMALL_R, ES_SMALL_ZC = 10.46, -13.0
INLET_EXT = 130.0           # straight inlet extension upstream of the spigot (~3 D)
ESIG_FLAT_SIDE = "down"     # as installed: spigot end high, S-curve drops toward the elbow

ES_END_Y = ELB_FACE_Y - ENGAGE_34          # E-Sig end face in the elbow frame
VENT_IN_X = ELB_SOCKET_SHOULDER_X

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", int(os.environ.get("GMSH_VERBOSE", "0")))
gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
occ = gmsh.model.occ


def esig_to_elbow(xE, zE):
    """E-Sig (x, z) -> elbow (x, y); E-Sig y -> elbow z. The E-Sig CAD has its flat
    side at z = 0 (top); with the flat side down it is turned 180 deg about its axis."""
    up = -1.0 if ESIG_FLAT_SIDE == "down" else 1.0
    return up * (zE - ES_SMALL_ZC), ES_END_Y + (ES_XE - xE)


# ------------------------------------------------ E-Sig inner S-curve from its STEP
vE = occ.importShapes(esig_step)
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


xs = [ES_XJ + i * (ES_XS - ES_XJ) / 25 for i in range(26)]
inner_bot = [(0.0, -46.47), (4.0, -46.45), (8.0, -46.32), (11.0, -46.1)]
inner_bot += [(x, outer_bottom(x) + ES_W) for x in xs]
occ.remove(vE, recursive=True)
occ.synchronize()


def ring_x(xE, r, zc):
    """Circle in the E-Sig cross-section plane at xE, placed in the elbow frame."""
    c = occ.addCircle(0, 0, 0, r)                      # in xy plane, normal z
    occ.rotate([(1, c)], 0, 0, 0, 0, 0, 1, math.pi / 2)    # move the seam off the symmetry plane
    occ.rotate([(1, c)], 0, 0, 0, 1, 0, 0, -math.pi / 2)   # normal -> +y
    X, Y = esig_to_elbow(xE, zc)
    occ.translate([(1, c)], X, Y, 0)
    return occ.addCurveLoop([c])


loops = []
for xE, zb in inner_bot:
    top = -ES_W
    r = 0.5 * (top - zb)
    loops.append(ring_x(xE, r, top - r))
loft = occ.addThruSections(loops, makeSolid=True, makeRuled=False)

# spigot + inlet extension (E-Sig x from ES_X0 - INLET_EXT to 0)
X0, Y0 = esig_to_elbow(ES_X0 - INLET_EXT, ES_SPIGOT_ZC)
L = 0.0 - (ES_X0 - INLET_EXT)
spig = occ.addCylinder(X0, Y0, 0, 0, -L, 0, ES_SPIGOT_R)
# small-bore straight section (E-Sig x 89.8 -> 197.4) ends at the E-Sig end face
Xs, Ys = esig_to_elbow(ES_XS - 0.5, ES_SMALL_ZC)
small = occ.addCylinder(Xs, Ys, 0, 0, -(ES_XE - ES_XS + 0.5), 0, ES_SMALL_R)

# ------------------------------------------------ elbow
cav = occ.addCylinder(0, ES_END_Y, 0, 0, ELB_SHOULDER_Y - ES_END_Y, 0, THREAD_MINOR_R)
tor = occ.addTorus(-ELB_BEND_R, ELB_SHOULDER_Y, 0, ELB_BEND_R, ELB_BORE_R, angle=math.pi / 2)
occ.rotate([(3, tor)], -ELB_BEND_R, ELB_SHOULDER_Y, 0, 0, 0, 1, -math.pi / 2)
leg_in = occ.addCylinder(0, ELB_SHOULDER_Y + 0.5, 0, 0, -0.5, 0, ELB_BORE_R)
leg_out = occ.addCylinder(-ELB_BEND_R + 0.5, 0, 0, VENT_IN_X - (-ELB_BEND_R + 0.5), 0, 0, ELB_BORE_R)

# ------------------------------------------------ current Venturi (revolved bore)
st = mwm.load_stations(os.path.join(HERE, "..", "..", "geometry", "bore_stations.csv"), "current")
wall = mwm.build_wall(st, 0.3)
pts = [occ.addPoint(VENT_IN_X - z, r, 0) for z, r, _ in wall]
z_out = VENT_IN_X - wall[-1][0]
pa = [occ.addPoint(z_out, 0, 0), occ.addPoint(VENT_IN_X, 0, 0)]
allp = pts + pa
lines = [occ.addLine(allp[i], allp[(i + 1) % len(allp)]) for i in range(len(allp))]
prof = occ.addPlaneSurface([occ.addCurveLoop(lines)])
vent = occ.revolve([(2, prof)], 0, 0, 0, 1, 0, 0, 2 * math.pi)
vent = [dt for dt in vent if dt[0] == 3]
# small overlap into the elbow so the fuse is clean
vent_lip = occ.addCylinder(VENT_IN_X + 0.3, 0, 0, -0.3, 0, 0, ELB_BORE_R)

parts = [(3, spig), (3, small), (3, cav), (3, tor), (3, leg_in), (3, leg_out), (3, vent_lip)] + loft + vent
if os.environ.get("DBG"):
    for d,t in parts: print("part",t,round(occ.getMass(d,t)),[round(v,1) for v in occ.getBoundingBox(d,t)])
halves = []
for dt in parts:
    bx = occ.addBox(-400, -400, -400, 800, 1200, 400)
    h, _ = occ.cut([dt], [(3, bx)])
    if os.environ.get("DBG"): print("half", dt, [round(occ.getMass(d, t)) for d, t in h])
    halves += h
fluid, _ = occ.fuse(halves[:1], halves[1:])
if os.environ.get("DBG"): print("fused", fluid, [round(occ.getMass(d,t)) for d,t in fluid])
occ.synchronize()
vols = gmsh.model.getEntities(3)
print("volumes", vols, "volume mm3", [round(occ.getMass(3, t)) for _, t in vols])

# ------------------------------------------------ classify boundary faces
y_inlet = Y0
patches = {"inlet": [], "outlet": [], "symmetry": [], "wall_esig": [], "wall_elbow": [], "wall_venturi": []}
for _, s in gmsh.model.getBoundary(vols, oriented=False):
    s = abs(s)
    b = gmsh.model.getBoundingBox(2, s)
    cx, cy, cz = occ.getCenterOfMass(2, s)
    if abs(b[2]) < 1e-3 and abs(b[5]) < 1e-3:
        patches["symmetry"].append(s)
    elif abs(b[1] - y_inlet) < 1e-4 and abs(b[4] - y_inlet) < 1e-4:
        patches["inlet"].append(s)
    elif abs(b[0] - z_out) < 1e-4 and abs(b[3] - z_out) < 1e-4:
        patches["outlet"].append(s)
    elif cx < VENT_IN_X - 0.05:
        patches["wall_venturi"].append(s)
    elif cy > ES_END_Y + 0.05:
        patches["wall_esig"].append(s)
    else:
        patches["wall_elbow"].append(s)
print({k: len(v) for k, v in patches.items()})
for k, v in patches.items():
    gmsh.model.addPhysicalGroup(2, v, name=k)

gmsh.option.setNumber("Mesh.Algorithm", 1)   # MeshAdapt: robust on the lofted S-curve
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)   # one solid per physical group
gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)       # STL in metres
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.write(os.path.join(out_dir, "fluid.stl"))
gmsh.option.setNumber("Mesh.ScalingFactor", 1.0)
gmsh.write(os.path.join(out_dir, "fluid_mm.stl"))

port_z_below = mwm.key_stations(st)["port"]
info = dict(
    frame="elbow frame, mm; Venturi axis = x, flow -x (down, gravity -x); threaded leg = y, flow -y; z>=0 half",
    engage_mm=ENGAGE_34, esig_flat_side=ESIG_FLAT_SIDE, esig_end_face_y=ES_END_Y, venturi_inlet_x=VENT_IN_X, outlet_x=z_out,
    inlet_y=y_inlet, inlet_centre_x=X0, inlet_radius=ES_SPIGOT_R, inlet_extension_mm=INLET_EXT,
    spigot_entry_y=esig_to_elbow(ES_X0, 0)[1],
    port_x=VENT_IN_X - port_z_below, port_side="unknown; wall pressure sampled all round",
    cavity_radius=THREAD_MINOR_R, elbow_bend_radius=ELB_BEND_R, elbow_bore_radius=ELB_BORE_R,
    venturi_wall=[(z, r) for z, r, _ in wall], fluid_volume_mm3_half=sum(occ.getMass(3, t) for _, t in vols),
    assumptions=[
        "E-Sig 3/4in NPT male engaged 14.0 mm into the elbow female leg (L1 8.6 mm + 3 wrench turns)",
        "Current Venturi spigot fully seated in the elbow socket; inlet face on the socket shoulder",
        "Female-thread cavity between E-Sig end and elbow shoulder modelled as a 23.8 mm cylinder",
        "E-Sig flat side down (spigot end high, S-curve drops to the elbow); inlet extension of 130 mm at the 43.9 mm spigot bore",
        "Half domain: symmetric about the bend plane (the chlorine port is not modelled)"])
json.dump(info, open(os.path.join(out_dir, "assembly_info.json"), "w"), indent=1)
gmsh.finalize()
