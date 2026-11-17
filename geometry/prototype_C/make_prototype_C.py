"""TuriTap handpump Venturi with integrated 1" long-radius elbow (prototype C).
Frame: Venturi axis = Z, flow in -Z. Bend exit plane at z=0. Inlet is horizontal, entering from +X.
"""
import cadquery as cq, numpy as np, math, json, sys

tan=lambda d: math.tan(math.radians(d))
P=dict(Dp=26.6, Rc=38.0, Ls=2.0, thc=20.0, Dt=14.0, Lt=14.0, thd=4.0, Dex=26.0, Ltail=0.89, wall=3.0,
       Lh=6.0, Lthr=18.0, Dboss=42.0, port_angle_deg=180.0, vane_t=1.2, vane_L=12.0)
# ---- 1" NPT internal thread (ASME B1.20.1): 11.5 TPI, taper 1:16 on diameter, E1 = 31.461 at the face
NPT=dict(P=25.4/11.5, E1=31.461, h=0.8*25.4/11.5, half=math.degrees(math.atan(1/32)))
# ---- injector port, as on the original part
RW0=6.131; PITCH=25.4/28; RMIN=2.7369; RMAJ=3.175; DEPTH=9.011; TILT=5.0
def fil(A,B,C,R):
    A,B,C=map(np.array,(A,B,C)); u1=(A-B)/np.linalg.norm(A-B); u2=(C-B)/np.linalg.norm(C-B)
    phi=math.acos(np.clip(u1@u2,-1,1)); t=R/math.tan(phi/2)
    T1,T2=B+u1*t,B+u2*t; b=(u1+u2)/np.linalg.norm(u1+u2); Cn=B+b*R/math.sin(phi/2); M=Cn+(B-Cn)/np.linalg.norm(B-Cn)*R
    return tuple(T1),tuple(M),tuple(T2)
def port_cutter():
    hole=cq.Workplane("XY").circle(RMIN).extrude(DEPTH+1.0)
    t30=tan(30); ri=RMIN-0.06; wo=PITCH/8; wi=wo+2*(RMAJ-ri)*t30; z0=-0.6
    helix=cq.Wire.makeHelix(PITCH,DEPTH+1.6,2.95,center=cq.Vector(0,0,z0))
    thr=cq.Workplane("XZ").polyline([(ri,z0-wi/2),(RMAJ,z0-wo/2),(RMAJ,z0+wo/2),(ri,z0+wi/2)]).close().sweep(cq.Workplane(obj=helix),isFrenet=True)
    thr=thr.intersect(cq.Workplane("XY").circle(RMAJ+0.3).extrude(DEPTH+1.0))
    cham=cq.Workplane("XY").workplane(offset=DEPTH-(3.47-RMIN)).circle(RMIN).workplane(offset=(3.47-RMIN)).circle(3.47).loft().union(
         cq.Workplane("XY").workplane(offset=DEPTH).circle(3.47).extrude(1.0))
    orf=cq.Workplane("XY").circle(0.8325).workplane(offset=-4.0).circle(0.8325-4.0*tan(1)).loft()
    return hole.union(thr).union(cham).union(orf).val()
def npt_cutter(L):
    """Local frame: +Z points out of the socket; z=0 at thread bottom, z=L at the face."""
    p,h,E1=NPT['P'],NPT['h'],NPT['E1']; k=1/32.0
    rmin=lambda z: (E1-h)/2-(L-z)*k; rmaj=lambda z: (E1+h)/2-(L-z)*k
    hole=cq.Workplane("XY").circle(rmin(0)).workplane(offset=L+1.0).circle(rmin(L+1.0)).loft()
    wo=0.033*p*2; wi0=wo+2*(h+0.08)*tan(30); z0=wi0/2+0.05; Hh=L+p-z0     # groove starts just above the thread bottom, runs out past the face
    helix=cq.Wire.makeHelix(p,Hh,(rmin(z0)+rmaj(z0))/2,center=cq.Vector(0,0,z0),angle=NPT['half'])
    ri=rmin(z0)-0.08; ro=rmaj(z0); wi=wo+2*(ro-ri)*tan(30)
    thr=cq.Workplane("XZ").polyline([(ri,z0-wi/2),(ro,z0-wo/2),(ro,z0+wo/2),(ri,z0+wi/2)]).close().sweep(cq.Workplane(obj=helix))
    cham=cq.Workplane("XY").workplane(offset=L-1.2).circle(rmin(L)-0.2).workplane(offset=1.2).circle(rmaj(L)+0.4).loft().union(
         cq.Workplane("XY").workplane(offset=L).circle(rmaj(L)+0.4).extrude(1.0))
    return [hole.val(),thr.val(),cham.val()], rmin(0)
def stations(p):
    Lc=(p['Dp']-p['Dt'])/2/tan(p['thc']); Ld=(p['Dex']-p['Dt'])/2/tan(p['thd'])
    z=np.cumsum([0,p['Ls'],Lc,p['Lt'],Ld,p['Ltail']]); r=np.array([p['Dp'],p['Dp'],p['Dt'],p['Dt'],p['Dex'],p['Dex']])/2
    return z,r
def build(p,name,vane=True):
    z,r=stations(p); L=z[-1]; Rc=p['Rc']; rp=p['Dp']/2; w=p['wall']; zp=z[2]+p['Lt']/2
    V=[(float(a),-float(b)) for a,b in zip(r,z)]
    # ---- fluid: Venturi (revolved) + bend (swept) + horizontal inlet bore
    wp=cq.Workplane("XZ").moveTo(0,0).lineTo(*V[0])
    for i,R in ((1,10.0),(2,14.0),(3,14.0)):
        T1,M,T2=fil(V[i-1],V[i],V[i+1],R); wp=wp.lineTo(*T1).threePointArc(M,T2)
    for v in V[4:]: wp=wp.lineTo(*v)
    vent=wp.lineTo(0,-L).close().revolve(360,(0,0,0),(0,1,0))
    s2=math.sqrt(0.5)
    path=cq.Workplane("XZ").moveTo(0,0).threePointArc((Rc*(1-s2),Rc*s2),(Rc,Rc))
    bend=lambda rad: cq.Workplane("XY").circle(rad).sweep(path)
    xin=Rc+p['Lh']+p['Lthr']                               # inlet face x
    horiz=lambda rad,x0,x1: cq.Workplane("YZ").workplane(offset=x0).center(0,Rc).circle(rad).extrude(x1-x0)
    cut_npt,rbot=npt_cutter(p['Lthr'])
    # 45 deg shoulder from thread bottom down to the flow bore
    sh=(rbot-rp)
    shoulder=cq.Workplane("YZ").workplane(offset=Rc+p['Lh']-sh).center(0,Rc).circle(rp).workplane(offset=sh).circle(rbot).loft()
    fluid=vent.union(bend(rp)).union(horiz(rp,Rc-0.01,Rc+p['Lh'])).union(shoulder)
    # ---- solid
    ro=lambda zz: float(np.interp(zz,z,r))+w
    wp=cq.Workplane("XZ").moveTo(0,0)
    for a,b in zip(r,z): wp=wp.lineTo(float(a)+w,-float(b))
    outer=wp.lineTo(0,-L).close().revolve(360,(0,0,0),(0,1,0))
    outer=outer.union(bend(rp+w)).union(horiz(rp+w,Rc-0.01,Rc+p['Lh']+2)).union(horiz(p['Dboss']/2,Rc+p['Lh']-sh-3.0,xin))
    # port boss + port, rotated to port_angle about Z (0 deg = +X = toward the inlet side)
    shift=p['Dt']/2-RW0; ct,st=math.cos(math.radians(TILT)),math.sin(math.radians(TILT))
    dn=cq.Vector(ct,0,st)                                   # inward direction for a port on the -X side
    wallpt=cq.Vector(-p['Dt']/2,0,-zp)                      # axis meets the bore wall at mid-throat
    seat=wallpt-dn*((7.9625-RW0)/ct); face=seat-dn*DEPTH; inner=wallpt+dn*1.2
    boss=cq.Workplane(cq.Plane(origin=face,xDir=cq.Vector(0,1,0),normal=dn)).circle(5.0).extrude((inner-face).Length).faces("<X").edges().fillet(0.5)
    loc=cq.Location(cq.Plane(origin=seat,xDir=cq.Vector(0,1,0),normal=dn*-1))
    pc=cq.Workplane(obj=port_cutter().moved(loc))
    rot=p['port_angle_deg']-180.0
    boss=boss.rotate((0,0,0),(0,0,1),rot); pc=pc.rotate((0,0,0),(0,0,1),rot)
    nloc=cq.Location(cq.Plane(origin=cq.Vector(xin-p['Lthr'],0,Rc),xDir=cq.Vector(0,1,0),normal=cq.Vector(1,0,0)))
    body=outer.union(boss).cut(fluid).cut(pc)
    for c in cut_npt: body=body.cut(cq.Workplane(obj=c.moved(nloc)))
    vanes=None
    if vane:
        for ang in (0,90):
            b=cq.Workplane("XY").workplane(offset=-L).rect(p['Dex']+1.0,p['vane_t']).extrude(p['vane_L']).rotate((0,0,0),(0,0,1),ang)
            vanes=b if vanes is None else vanes.union(b)
        body=body.union(vanes)
    s=body.val(); cq.exporters.export(body,name+'.step')
    # fluid volume for CFD: add the pipe bore inside the socket up to the inlet face, remove vanes
    fl=fluid.union(horiz(rp,Rc+p['Lh']-0.01,xin))
    if vane: fl=fl.cut(vanes)
    cq.exporters.export(fl,name+'_fluid.step'); cq.exporters.export(fl,name+'_fluid.stl',tolerance=0.02,angularTolerance=0.1)
    bb=s.BoundingBox()
    info=dict(name=name,valid=bool(s.isValid()),solids=len(body.solids().vals()),volume_cm3=round(s.Volume()/1000,2),
              fluid_valid=bool(fl.val().isValid()),fluid_volume_cm3=round(fl.val().Volume()/1000,2),
              inlet_face_x=round(xin,2),inlet_centreline_z=Rc,outlet_face_z=round(-L,2),height_centreline_to_outlet=round(Rc+L,2),
              port_axis_z=round(-zp,2),bbox=[round(v,1) for v in (bb.xmin,bb.xmax,bb.ymin,bb.ymax,bb.zmin,bb.zmax)],
              stations_z=[round(-float(a),2) for a in z],stations_D=[round(2*float(a),2) for a in r],params=p)
    print(json.dumps({k:v for k,v in info.items() if k!='params'}))
    return body,fl,info
if __name__=='__main__':
    body,fl,info=build(P,'TuriTap_HP_C_integrated_elbow_1in_crossvane',True)
    json.dump(info,open('prototype_C_parameters.json','w'),indent=1)
    half=body.cut(cq.Workplane("XY").box(300,100,400).translate((0,-50,0)))
    cq.exporters.export(half,'C_section.svg',opt={"projectionDir":(0,-1,0),"showHidden":False,"width":1100,"height":1500})
