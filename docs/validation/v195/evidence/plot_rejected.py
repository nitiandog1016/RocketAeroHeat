"""原箭体单步失败场诊断图；明确拒绝状态、壁邻身份和单位。"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
base=Path(__file__).resolve().parent;out=base/'figure02';out.mkdir(exist_ok=False)
s=json.loads((base/'rocket01/surface.json').read_text());v=np.array(s['vertices']);f=np.array(s['facesOneBased'])-1
d=np.genfromtxt(base/'audit01/rejected_original_wall_diagnostic.csv',delimiter=',',names=True)
assert np.array_equal(d['face_id'],np.arange(1,19753))
fig=plt.figure(figsize=(13,3.5),layout='constrained')
for i,(key,label,scale) in enumerate([('p_raw_Pa','Adjacent-cell pressure [kPa]',1e-3),('T_raw_fluid_K','Adjacent-cell fluid temperature [K]',1)]):
    ax=fig.add_subplot(1,2,i+1,projection='3d');values=d[key]*scale;norm=Normalize(values.min(),values.max())
    collection=Poly3DCollection(v[f],facecolors=plt.cm.viridis(norm(values)),edgecolors='none',rasterized=True)
    ax.add_collection3d(collection);lo=v.min(0);hi=v.max(0);ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_zlim(lo[2],hi[2])
    ax.set_box_aspect(np.ptp(v,axis=0));ax.view_init(elev=-25,azim=-72)
    ax.set_proj_type('ortho');ax.set_axis_off();ax.set_title(label,fontsize=10)
    ax.text2D(.1,.1,'Original geometry: length 2.878 m | nose +X (right)',transform=ax.transAxes,fontsize=9)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),ax=ax,shrink=.65,pad=.02)
fig.suptitle('M5 / alpha 20 deg | Rejected first-step diagnostic\nInvalid pseudo-time step; not converged; q and solid wall temperature not solved',fontsize=12)
fig.savefig(out/'rejected_rocket_first_step.png',dpi=160)
print(str(out/'rejected_rocket_first_step.png'))
