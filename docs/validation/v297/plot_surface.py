"""用同一份原面表绘制热敏珠的表面放热与珠体、引线的局部通量范围。"""
import io
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np

HERE=Path(__file__).resolve().parent
with zipfile.ZipFile(HERE.parents[2]/"cases/thermistor_b/leads3.zip") as archive:
    body=json.loads(archive.read("body.json"))
    wall=np.genfromtxt(io.BytesIO(archive.read("wall.csv")),names=True,delimiter=",")
v=np.asarray(body["V"])*1000
f=np.asarray(body["F"])-1
outward=-wall["heatFluxWm2"]/1000
assert np.array_equal(wall["faceID"],np.arange(len(f)))
plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
fig=plt.figure(figsize=(10.6,5.25),facecolor="white")
ax=fig.add_axes([.045,.255,.54,.54])
side=fig.add_axes([.695,.185,.26,.64])
fig.text(.05,.91,"带引线 B 型热敏珠的表面放热",fontsize=18,weight="bold")
fig.text(.05,.85,"3 m/s 空气，293.15 K    ·    等温壁面 305.15 K    ·    轴向对称稳态",color="#485464")
view=np.array([.24,.91,.34]);view/=np.linalg.norm(view)
horizontal=np.array([1.,0.,0.])-view[0]*view;horizontal/=np.linalg.norm(horizontal)
up=np.cross(horizontal,view)
depth=(v[f]@view).mean(axis=1);order=np.argsort(depth)
norm=Normalize(0,20);cmap=plt.get_cmap("viridis")
project=np.column_stack((v@horizontal,v@up))
ax.add_collection(PolyCollection(project[f[order]],facecolors=cmap(norm(outward[order])),
                                edgecolors="face",linewidths=.22,antialiased=True))
ax.autoscale();ax.set_aspect("equal");ax.set_axis_off()
ax.annotate("来流",xy=(.39,.39),xytext=(.69,.39),arrowprops=dict(arrowstyle="->",color="#384858"),
            va="center",ha="center",color="#384858")
bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.12,.17,.39,.032]),orientation="horizontal")
bar.set_ticks([0,5,10,15,20]);bar.set_label("向空气放热通量  (kW/m²)")
for i,name in enumerate(body["regionNames"],1):
    pick=wall["regionID"]==i
    values=outward[pick]
    mean=np.average(values,weights=wall["areaM2"][pick])
    side.plot([values.min(),values.max()],[i,i],color="#91a6ae",linewidth=3,solid_capstyle="round")
    side.scatter(mean,i,s=30,color="#184d5a",zorder=3)
side.set_yticks(range(1,8));side.set_yticklabels(body["regionNames"])
side.set_ylim(.4,7.6);side.set_xlim(0,20);side.set_xticks([0,5,10,15,20])
side.grid(axis="x",color="#e4e9eb");side.set_axisbelow(True)
side.set_xlabel("向空气放热通量  (kW/m²)")
side.set_title("各区域范围与面积加权平均值",fontsize=11,loc="left",pad=8)
for key in ("top","right","left"):side.spines[key].set_visible(False)
side.tick_params(axis="y",length=0);side.spines["bottom"].set_color("#c4cbce")
fig.text(.05,.04,"细网首层 0.937 µm    ·    10,272 个原面    ·    图中向空气放热取正，工作台入壁热流取负",fontsize=9,color="#485464")
fig.savefig(HERE/"evidence/heat_flux.png",dpi=155)
