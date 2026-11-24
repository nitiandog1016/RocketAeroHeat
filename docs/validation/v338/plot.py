"""绘制背风两站的整圈排序、赤道缺口及七站取样范围。"""
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
with ZipFile(HERE.parents[2]/"cases/short_sphere/stationAzimuth65.zip") as z:
    summary=json.loads(z.read("summary.json"));rows=list(csv.DictReader(io.StringIO(z.read("points.csv").decode())))
plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False})
fig=plt.figure(figsize=(10.7,6.1),constrained_layout=True)
grid=fig.add_gridspec(2,2,width_ratios=[1.35,1]);rear=fig.add_subplot(grid[0,0]);eq=fig.add_subplot(grid[1,0]);span=fig.add_subplot(grid[:,1])
blue="#2b6081";brown="#aa5d31";green="#33765d"
for a,ax,color in ((119.8,rear,blue),(166.7,rear,brown),(90,eq,green)):
    band=[r for r in rows if float(r["angleDeg"])==a];phi=[float(r["phiDeg"]) for r in band]+[360.]
    nu=[float(r["nu"]) if r["nu"] else np.nan for r in band];nu.append(nu[0])
    ax.plot(phi,nu,label=f"{a:g}°",color=color,lw=1.8)
    station=next(s for s in summary["stations"] if s["angleDeg"]==a)
    ax.axhline(station["experimentNu"],color=color,ls="--",lw=.95)
proof=json.loads((HERE/"evidence/comparison.json").read_bytes())
for i,(a,b) in enumerate(proof["equatorGapDegrees"]):
    eq.axvspan(a,b,color="#d5dadd",alpha=.85,label="遮盖或归属不一致" if i==0 else None)
for ax in (rear,eq):
    ax.set(xlim=(0,360),xticks=np.arange(0,361,90),ylabel="局部 Nu");ax.grid(axis="y",color="#e0e4e7")
rear.set_ylim(2,22);rear.set_title("背风两站 · 虚线为实验 Nu",loc="left");rear.legend(frameon=False,ncol=2,loc="upper center")
eq.set_ylim(3,24);eq.set_title("90° 赤道截线",loc="left");eq.legend(frameon=False,loc="lower center",fontsize=9,ncol=2);eq.set_xlabel("方位角 (°)，从 +Z 转向 +Y")
for i,s in enumerate(summary["stations"]):
    lo=s["nu"]["minimum"]["value"];hi=s["nu"]["maximum"]["value"]
    span.plot([lo,hi],[i,i],lw=3,color=blue,label="已采到的范围" if i==0 else None,solid_capstyle="round")
    span.plot(s["plusZNu"],i,"o",ms=5,color=blue,label="+Z 测线" if i==0 else None)
    span.plot(s["experimentNu"],i,"D",color=brown,ms=5,label="Short 表 III" if i==0 else None)
span.set_yticks(range(7));span.set_yticklabels([f'{s["angleDeg"]:g}°' for s in summary["stations"]]);span.invert_yaxis()
span.set(xlabel="局部 Nu",ylabel="测站极角",xlim=(0,46));span.set_title("外露球面的周向换热范围",loc="left")
span.grid(axis="x",color="#e0e4e7");span.legend(frameon=False,loc="lower right")
for ax in (rear,eq,span):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False);ax.set_axisbelow(True)
fig.suptitle("Short Test 65 双支管银球    ·    Re = 882    ·    扰动后 2.565681 ms",fontsize=15,fontweight="bold")
fig.savefig(HERE/"evidence/azimuth.png",dpi=165,facecolor="white")
