"""绘制同一物理时刻的七站主解、独立梯度，以及相对主解的差。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
model=json.loads((HERE.parents[2]/"cases/short_sphere/stationGradients65.json").read_bytes())
rows=next(s["rows"] for s in model["sets"] if s["label"]=="half")
def values(field):return np.array([r[field] for r in rows])

plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
fig,(main,delta)=plt.subplots(2,1,figsize=(8.4,6.3),sharex=True,gridspec_kw={"height_ratios":[2.8,1.1]})
x=values("angleDeg");q=values("reactionNu");g=values("gradientNu");ref=values("experimentNu")
main.plot(x,q,color="#294d69",lw=2,label="壁反力")
main.plot(x,g,"o",fillstyle="none",color="#26795e",ms=7,mew=1.7,label="法向温度梯度")
main.scatter(x,ref,color="#995634",s=43,marker="D",label="Short 表 III",zorder=4)
main.set(ylabel="局部 Nusselt 数",ylim=(0,46))
main.legend(loc="upper right",frameon=False)
delta.axhline(0,color="#999999",lw=.8)
delta.plot(x,g-q,"o-",color="#26795e",lw=1.5,ms=4)
delta.set(xlabel="从迎风点沿 +Z 测线的极角 (°)",ylabel="梯度 Nu − 壁反力 Nu",ylim=(-.005,.056))
for ax in (main,delta):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="y",color="#dedede",lw=.6);ax.set_axisbelow(True)
delta.set_xticks(x);delta.set_xticklabels([f"{v:g}" for v in x]);delta.set_xlim(-5,174)
fig.suptitle("双支管银球七站换热  ·  Re = 882  ·  扰动后 2.565681 ms",fontweight="bold",fontsize=14)
fig.subplots_adjust(left=.14,right=.96,top=.91,bottom=.11,hspace=.18)
fig.savefig(HERE/"evidence/stations.png",dpi=170,facecolor="white")
