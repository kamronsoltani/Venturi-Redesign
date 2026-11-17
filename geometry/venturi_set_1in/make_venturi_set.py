"""Four 1" handpump Venturis for the user's reducing elbow (socket dia 33.8 x 18 deep, bore dia 31.2 at the shoulder).
Frame: axis Z, inlet face z=0, flow in -Z, chlorine port on the -X side."""
import cadquery as cq, numpy as np, math, json
tan=lambda d: math.tan(math.radians(d))
# ---- chlorine port: 1/4-28 UNF, thread form moved out radially by 0.10 mm (0.2 on diameter)
CLR=0.10; PITCH=25.4/28; RMIN=2.7369+2*CLR; RMAJ=3.175+CLR; DEPTH=9.011; TILT=5.0; LAND=(7.9625-6.131)   # radial land, as original
RCH=3.47+CLR; RORF=0.835
def port_parts():
    """Same construction as the earlier 3/4" files: one cutter = hole + thread + chamfer + orifice.
    Local frame: z=0 flat seat, +z out of the port. Minor dia opened 0.4 so the thread crest keeps its original 0.29 mm width."""
    t30=tan(30); ri=RMIN-0.06; wo=PITCH/8+2*CLR*t30; wi=wo+2*(RMAJ-ri)*t30; z0=-0.6
    hole=cq.Workplane("XY").circle(RMIN).extrude(DEPTH+1.0)
    helix=cq.Wire.makeHelix(PITCH,DEPTH+1.6,2.95,center=cq.Vector(0,0,z0))
    thr=cq.Workplane("XZ").polyline([(ri,z0-wi/2),(RMAJ,z0-wo/2),(RMAJ,z0+wo/2),(ri,z0+wi/2)]).close().sweep(cq.Workplane(obj=helix),isFrenet=True)
    thr=thr.intersect(cq.Workplane("XY").circle(RMAJ+0.3).extrude(DEPTH+1.0))
    cham=cq.Workplane("XY").workplane(offset=DEPTH-(RCH-RMIN)).circle(RMIN).workplane(offset=(RCH-RMIN)).circle(RCH).loft().union(
         cq.Workplane("XY").workplane(offset=DEPTH).circle(RCH).extrude(1.0))
    orf=cq.Workplane("XY").workplane(offset=-5.0).circle(RORF).extrude(5.01)     # plain dia 1.67 through-hole
    whole=hole.union(thr).union(cham).union(orf).val()
    parts=[hole.val(),thr.val(),cham.val(),orf.val()]
    print('port cutter: union vol %.1f, parts %s, crest width %.3f'%(whole.Volume(),[round(c.Volume(),1) for c in parts],PITCH-(wo+2*(RMAJ-RMIN)*t30)))
    assert whole.Volume()>parts[0].Volume()+0.5*parts[1].Volume(), 'union dropped a part'
    return [('port',whole)]
PORT=port_parts()
def fil(A,B,C,R):
    A,B,C=map(np.array,(A,B,C)); u1=(A-B)/np.linalg.norm(A-B); u2=(C-B)/np.linalg.norm(C-B)
    phi=math.acos(np.clip(u1@u2,-1,1)); t=R/math.tan(phi/2)
    T1,T2=B+u1*t,B+u2*t; b=(u1+u2)/np.linalg.norm(u1+u2); Cn=B+b*R/math.sin(phi/2); M=Cn+(B-Cn)/np.linalg.norm(B-Cn)*R
    return tuple(T1),tuple(M),tuple(T2)
def make(name,Dt,Dex,L,blades,Din=31.20,thc=20.0,Lt=14.0,thd=4.0,w=3.0,Rs=16.70,Ls=18.0):
    Lc=(Din-Dt)/2/tan(thc); Ld=(Dex-Dt)/2/tan(thd)
    z=np.cumsum([0,Lc,Lt,Ld]); r=np.array([Din,Dt,Dt,Dex])/2
    assert z[-1]<=L+1e-6, (name,z[-1])
    if L-z[-1]>1e-6: z=np.append(z,L); r=np.append(r,Dex/2)
    V=[(float(a),-float(b)) for a,b in zip(r,z)]; zp=z[1]+Lt/2
    wp=cq.Workplane("XZ").moveTo(0,0).lineTo(*V[0])
    for i in (1,2):
        T1,M,T2=fil(V[i-1],V[i],V[i+1],14.0); wp=wp.lineTo(*T1).threePointArc(M,T2)
    for v in V[3:]: wp=wp.lineTo(*v)
    bore=wp.lineTo(0,-L).close().revolve(360,(0,0,0),(0,1,0))
    ro=lambda zz: float(np.interp(zz,z,r))+w
    O=[(Rs,0),(Rs,-Ls),(ro(Ls),-Ls)]+[(float(a)+w,-float(b)) for a,b in zip(r,z) if b>Ls+1e-6]
    wp=cq.Workplane("XZ").moveTo(0,0)
    for o in O: wp=wp.lineTo(*o)
    body=wp.lineTo(0,-L).close().revolve(360,(0,0,0),(0,1,0))
    ct,st=math.cos(math.radians(TILT)),math.sin(math.radians(TILT)); dn=cq.Vector(ct,0,st)      # inward
    wallpt=cq.Vector(-Dt/2,0,-zp); seat=wallpt-dn*(LAND/ct); face=seat-dn*DEPTH; inner=wallpt+dn*1.5
    boss=cq.Workplane(cq.Plane(origin=face,xDir=cq.Vector(0,1,0),normal=dn)).circle(5.0).extrude((inner-face).Length).faces("<X").edges().fillet(0.5)
    body=body.union(boss).cut(bore)
    loc=cq.Location(cq.Plane(origin=seat,xDir=cq.Vector(0,1,0),normal=dn*-1)); rem={}
    for nm,c in PORT:
        v0=body.val().Volume(); body=body.cut(cq.Workplane(obj=c.moved(loc))); rem[nm]=round(v0-body.val().Volume(),1)
    if blades:
        for k in range(blades//2):
            body=body.union(cq.Workplane("XY").workplane(offset=-L).rect(Dex+1.0,1.2).extrude(12).rotate((0,0,0),(0,0,1),k*360.0/blades))
    s=body.val(); out=dn*-1; side=cq.Vector(0,1,0)
    grid={}
    for d in (0.3,0.9,1.5):                      # inside the land: only the orifice may be open
        grid[d]=''.join('.' if not s.isInside(wallpt+out*(d/ct)+side*o) else '#' for o in (0,0.6,1.1,2.0,2.9))
    for d in (2.3,5.0,9.0):                      # above the seat: hole open, wall solid
        grid[d]=''.join('.' if not s.isInside(wallpt+out*(d/ct)+side*o) else '#' for o in (0,0.6,1.1,2.0,2.7,3.6))
    ok=all(grid[d]=='..###' for d in (0.3,0.9,1.5)) and all(grid[d]=='.....#' for d in (2.3,5.0,9.0)) and all(v>0 for v in rem.values())
    es=[e.Length() for f in s.Faces() if f.Center().x<-(Dt/2+0.5) and abs(f.Center().z+zp)<7 for e in f.Edges()]
    cq.exporters.export(body,name+'.step')
    info=dict(name=name,Dt=Dt,Dex=Dex,L=L,blades=blades,valid=bool(s.isValid()),solids=len(body.solids().vals()),port_ok=ok,removed=rem,
              min_port_edge_mm=round(min(es),3),port_z=round(zp,2),stations=[(round(float(a),2),round(2*float(b),2)) for a,b in zip(z,r)],
              column_below_port_mm=round(L-zp,1),xmin=round(s.BoundingBox().xmin,2))
    print(json.dumps(info)); return body,info
if __name__=='__main__':
    out=[]
    for a in (('TuriTap_HP_V1_throat14_exit26',14,26,125.0,4),('TuriTap_HP_V2_throat16_exit26',16,26,125.0,4),
              ('TuriTap_HP_V3_throat18_exit26',18,26,125.0,4),('TuriTap_HP_V4_throat14_exit30_long',14,30,153.0,6)):
        b,i=make(*a); out.append(i)
        if a[1]==14 and a[2]==26: B1=b
    json.dump(out,open('venturi_set_parameters.json','w'),indent=1)
    zp=out[0]['port_z']
    box=cq.Workplane("XY").box(16,60,14).translate((-12.5,0,-zp-0.3))
    half=B1.intersect(box).cut(cq.Workplane("XY").box(100,100,300).translate((0,-50,-60)))
    cq.exporters.export(half,'port_sec.svg',opt={"projectionDir":(0,-1,0),"showHidden":False,"width":1200,"height":1100})
