"""用已归档的真实测站位置、两档剪切和试验 Nu 生成共用横坐标的对照图。"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
d=json.loads((HERE.parents[2]/"cases/short_sphere/supportedDirection65.json").read_bytes())
plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"font.size":10,"axes.unicode_minus":False})
fig,(heat,shear)=plt.subplots(2,1,sharex=True,figsize=(8.1,5.5),constrained_layout=True)
f=d["fine"]["stations"];x=[r["angleDeg"] for r in f]
heat.plot(x,[r["measuredNu"] for r in f],"D",ms=5,color="#994321",label="Short 表 III 实测")
for key,label,style,color in (("coarse","五层网","--o","#598080"),("fine","十层网","-+","#172f58")):
    stations=d[key]["stations"]
    heat.plot(x,[r["computedNu"] for r in stations],style,color=color,ms=6,lw=1.3,label=label)
    bracket=d[key]["transitions"][0]
    points=sorted(stations+[bracket["left"],bracket["right"]],key=lambda r:r["angleDeg"])
    shear.plot([r["angleDeg"] for r in points],[np.nan if r["shearPa"] is None else r["shearPa"] for r in points],
               style,color=color,ms=6,lw=1.3,label=label)
heat.set_title("Short Test 65   |   双支管银球 · Re = 882",loc="left",fontweight="bold")
heat.set_ylabel("当地努塞尔数 Nu");heat.legend(frameon=False,ncol=3)
heat.annotate("实测 5.5",(119.8,5.5),xytext=(-9,2),textcoords="offset points",ha="right")
heat.annotate("实测 8.8",(166.7,8.8),xytext=(0,7),textcoords="offset points",ha="center")
shear.axhline(0,color="#757575",lw=.9)
b=d["fine"]["transitions"][0]
lo,hi=b["left"]["angleDeg"],b["right"]["angleDeg"]
shear.axvspan(lo,hi,color="#d8d0c4",alpha=.8)
shear.annotate("原生节点变号区间\n143.44°–149.06°",((lo+hi)/2,.064),xytext=(103,.084),
               arrowprops={"arrowstyle":"-","color":"#746b61"},fontsize=9)
shear.text(2,.012,"0° 方向未定义",color="#555555",fontsize=9)
shear.set_ylabel("沿流壁剪切 (Pa)");shear.set_xlabel("从迎风点沿 +Z 测线的极角 (°)")
shear.set_xticks([0,30,60,90,119.8,150,166.7]);shear.set_xlim(-3,181)
for ax in (heat,shear):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="y",color="#e4e6e7",lw=.65);ax.set_axisbelow(True)
fig.savefig(HERE/"evidence/stations.png",dpi=180)
