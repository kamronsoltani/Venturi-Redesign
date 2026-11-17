# showfluid.py <fluid_mm.stl> <out.png>  -- quick picture of a water volume from build_fluid.py:
# symmetry-plane section plus a 3-D view coloured by patch (check orientation before meshing).
import numpy as np,sys,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sol={}; cur=None; v=[]
for l in open(sys.argv[1]):
    l=l.strip()
    if l.startswith('solid'): cur=l.split()[1]; sol[cur]=[]
    elif l.startswith('vertex'): sol[cur].append(list(map(float,l.split()[1:])))
fig,ax=plt.subplots(1,2,figsize=(16,9),gridspec_kw=dict(width_ratios=[1,1.3]))
S=np.array(sol['symmetry']).reshape(-1,3,3)
ax[0].add_collection(PolyCollection(S[:,:,[1,0]],facecolor='#9cc3e6',edgecolor='none'))
ax[0].autoscale(); ax[0].set_aspect('equal'); ax[0].set_xlabel('y (mm, flow from right)'); ax[0].set_ylabel('x (mm, up)'); ax[0].set_title('symmetry plane (water region)'); ax[0].grid(alpha=.3)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ax3=fig.add_subplot(1,2,2,projection='3d'); ax[1].axis('off')
cols=dict(inlet='#2e7d32',outlet='#c62828',wall_esig='#4a6fa5',wall_elbow='#e0a030',wall_venturi='#7b5aa6')
for k,c in cols.items():
    T=np.array(sol[k]).reshape(-1,3,3)
    n=np.cross(T[:,1]-T[:,0],T[:,2]-T[:,0]); n/=np.linalg.norm(n,axis=1)[:,None]+1e-12
    sh=0.45+0.55*np.abs(n@np.array([0.4,0.5,0.75]))
    rgb=np.array(matplotlib.colors.to_rgb(c)); ax3.add_collection3d(Poly3DCollection(T[:,:,[1,2,0]],facecolors=np.clip(sh[:,None]*rgb,0,1),edgecolor='none'))
P=np.concatenate([np.array(v) for v in sol.values()]); c=(P.max(0)+P.min(0))/2; r=(P.max(0)-P.min(0)).max()/2
ax3.set_xlim(c[1]-r,c[1]+r); ax3.set_ylim(c[2]-r,c[2]+r); ax3.set_zlim(c[0]-r,c[0]+r); ax3.view_init(18,-120)
ax3.set_title('green inlet, blue E-Sig, orange elbow, purple Venturi, red outlet')
plt.tight_layout(); plt.savefig(sys.argv[2],dpi=90)
