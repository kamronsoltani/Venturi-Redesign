"""Rebuild the E-Sig spout reducer with a 1in NPT male end (gmsh/OpenCASCADE).
Keeps the E-Sig spigot, collar and fillets up to x=14 mm, reuses its eccentric S-curve,
and replaces the 3/4in end with a 1in Sch40 bore and an NPT-tapered blank (thread to be added in CAD).
Usage: python3 make_esig_1in.py <E-Sig 3/4in STEP>
"""
import gmsh, math, sys, json
gmsh.initialize(); gmsh.option.setNumber("General.Terminal",0)
gmsh.option.setString("Geometry.OCCTargetUnit","MM")
occ=gmsh.model.occ
occ.importShapes(sys.argv[1] if len(sys.argv)>1 else "esig34.step")  # the 3/4in E-Sig STEP; occ.synchronize()
V0=(3,1)
def section_bottom(X):
    r=occ.addRectangle(-100,-100,0,200,200); occ.rotate([(2,r)],0,0,0,0,1,0,-math.pi/2); occ.translate([(2,r)],X,0,0)
    out,_=occ.intersect([V0],[(2,r)],removeObject=False,removeTool=True); occ.synchronize()
    zb=min(gmsh.model.getBoundingBox(2,t)[2] for d,t in out)
    occ.remove(out,recursive=True); occ.synchronize(); return zb
# ---------------- parameters (mm) ----------------
W   = 2.54                 # wall, same as E-Sig (0.100 in)
XJ  = 14.0                 # join: keep original (spigot, collar, fillets) for x <= XJ
XS  = 89.8                 # end of eccentric transition (same as E-Sig)
XE  = 197.4                # end face (same overall length as E-Sig)
IDs = 26.64                # 1" Sch40 bore (1.049 in)
ODs = IDs+2*W              # 31.72
# 1" NPT male, ASME B1.20.1
P   = 25.4/11.5            # pitch 2.2087
E0  = 1.21363*25.4         # pitch dia at small end, 30.826
LT  = 0.9845*25.4          # L4 total thread length, 25.0
HT  = 1.60                 # modelled thread depth (nominal 0.8P = 1.77; slightly shallow for printing)
XT  = XE-LT
# eccentric S-curve: reuse E-Sig's own bottom-line shape
xs=[XJ+i*(XS-XJ)/24 for i in range(25)]
zb0=[section_bottom(x) for x in xs]
f=[(z-zb0[0])/(zb0[-1]-zb0[0]) for z in zb0]
R0=-zb0[0]/2                       # outer radius at the join (flat top at z=0)
Rs=ODs/2
def circ(x,r,zc):
    c=occ.addCircle(0,0,0,r); occ.rotate([(1,c)],0,0,0,0,1,0,math.pi/2); occ.translate([(1,c)],x,0,zc); return c
def ring(x,r):
    return occ.addCurveLoop([circ(x,r,-r)])
outer=[ring(x,R0+(Rs-R0)*fi) for x,fi in zip(xs,f)]
def ring_c(x,r,zc):
    return occ.addCurveLoop([circ(x,r,zc)])
inner=[ring_c(x,R0+(Rs-R0)*fi-W,-(R0+(Rs-R0)*fi)) for x,fi in zip(xs,f)]
vo=occ.addThruSections(outer,makeSolid=True,makeRuled=False)
vi=occ.addThruSections(inner,makeSolid=True,makeRuled=False)
zc=-Rs
# straight section and thread core (root cone)
pipe=occ.addCylinder(XS,0,zc,XT-XS,0,0,Rs)
rp=lambda x: E0/2+(XE-x)/32.0            # pitch radius (taper 1:16 on dia)
crest=lambda x: rp(x)+0.8*P/2
core_r=lambda x: crest(x)-HT
core=occ.addCone(XT-1.0,0,zc,XE-XT+1.0,0,0,crest(XT-1.0),crest(XE))   # NPT blank at major diameter
def chamtool():
    a=occ.addCone(XE-6,0,zc,6,0,0,crest(XE)+5,crest(XE)-1.0); b=occ.addCylinder(XT-3,0,zc,XE-6-(XT-3),0,0,crest(XE)+8)
    t,_=occ.fuse([(3,a)],[(3,b)]); return t
core_c,_=occ.intersect([(3,core)],chamtool())
thr=core_c
new,_=occ.fuse(vo,[(3,pipe)])
new,_=occ.fuse(new,thr)
new,_=occ.cut(new,vi)
bore=occ.addCylinder(XS-0.5,0,zc,XE-XS+2,0,0,IDs/2)
new,_=occ.cut(new,[(3,bore)])
# keep original up to the join
kbox=occ.addBox(-100,-100,-100,100+XJ,200,200)
keep,_=occ.intersect([V0],[(3,kbox)])
part,_=occ.fuse(keep,new)
occ.synchronize()
vols=gmsh.model.getEntities(3); print("volumes",vols)
for d,t in vols: print(t, round(occ.getMass(3,t)), [round(v,1) for v in gmsh.model.getBoundingBox(3,t)])
mass=sum(occ.getMass(3,t) for d,t in vols); print("solid volume mm3",round(mass))
print("bbox",[round(v,2) for v in gmsh.model.getBoundingBox(-1,-1)])
gmsh.write("TuriTap_ESig_reducer_1in_NPT.step")
gmsh.option.setNumber("Mesh.MeshSizeMax",1.2); gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",16)
gmsh.model.mesh.generate(2); gmsh.write("esig1.stl")
json.dump(dict(W=W,XJ=XJ,XS=XS,XE=XE,ID_small=IDs,OD_small=ODs,NPT_pitch=P,E0=E0,L4=LT,thread_depth=HT,thread_start=XT,
  R_join=R0,f=f,xs=xs),open("esig1_params.json","w"),indent=1)
gmsh.finalize()
