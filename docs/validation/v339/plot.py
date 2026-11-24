"""比较双支管银球晚时刻的绝对入壁热流、同面时间差和背风站的周向变化。"""
import json
from pathlib import Path
import sys
from zipfile import ZipFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent;CASE=HERE.parents[2]/"cases/short_sphere";LATE=CASE/"later1"
sys.path.insert(0,str(HERE.parent/"v329"))
from plot_surface import view


def curve(file):
    with ZipFile(file) as z:return np.genfromtxt(z.open("points.csv"),delimiter=",",names=True,dtype=None,encoding="utf8")


if __name__=="__main__":
    pair=json.loads((LATE/"supportedTime65.json").read_bytes());source=json.loads((CASE/"supportedTime65.json").read_bytes())
    with ZipFile(CASE/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    with ZipFile(LATE/"supportedTime65Half.zip") as z:w=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")
    with ZipFile(LATE/"timeDifference.zip") as z:delta=np.genfromtxt(z.open("half.csv"),names=True,delimiter=",")
    pts=np.asarray(g["V"])*1000;tri=np.asarray(g["F"])-1
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False,"axes.unicode_minus":False})
    fig=plt.figure(figsize=(10.6,7),constrained_layout=True);grid=fig.add_gridspec(3,2,height_ratios=[1.1,.09,1.05])
    q=w["fluxInWm2"]/1000
    a=view(fig.add_subplot(grid[0,0]),pts,tri,q,[-.9,.22,.3],"10.263 ms · 细步全表面入壁热流",plt.Normalize(np.floor(q.min()),0))
    a.set_cmap("viridis");fig.colorbar(a,cax=fig.add_subplot(grid[1,0]),orientation="horizontal").set_label("入壁热流 (kW/m²)")
    dq=delta["deltaFluxInWm2"];limit=max(abs(dq))
    b=view(fig.add_subplot(grid[0,1]),pts,tri,dq,[-.9,.22,.3],"10.263 ms 减 2.566 ms · 同一细步分支",plt.Normalize(-limit,limit))
    b.set_cmap("RdBu_r");fig.colorbar(b,cax=fig.add_subplot(grid[1,1]),orientation="horizontal").set_label("入壁热流变化 (W/m²)")
    x=fig.add_subplot(grid[2,0]);colors=["#2c6282","#a9572d"]
    for (branch,dt),color in zip((("full",.25),("half",.125)),colors):
        chain=[r for r in pair["history"][branch] if r["elapsedDU"]>=.25]
        x.plot([r["elapsedS"]*1000 for r in chain],[(r["wallW"]-chain[0]["wallW"])*1e9 for r in chain],"o-",color=color,label=f"Δt U/D = {dt:g}",ms=4)
    x.set(xlabel="扰动后时间 (ms)",ylabel="向外放热变化 (nW)",title="相对各自 2.566 ms 终点的变化");x.legend(frameon=False)
    x=fig.add_subplot(grid[2,1]);late=curve(LATE/"stationAzimuth.zip");early=curve(CASE/"stationAzimuth65.zip")
    for angle,col in zip((119.8,166.7),colors):
        left=early[early["angleDeg"]==angle];right=late[late["angleDeg"]==angle]
        np.testing.assert_allclose(left["phiDeg"],right["phiDeg"],rtol=0,atol=1e-10)
        y=right["nu"]-left["nu"];x.plot(np.r_[right["phiDeg"],360],np.r_[y,y[0]],color=col,label=f"{angle:g}°")
    x.set(xlabel="方位角 (°)，从 +Z 转向 +Y",ylabel="局部 Nu 变化",title="两个背风极角 · 同方位细步时间差",xlim=(0,360),xticks=[0,90,180,270,360])
    x.ticklabel_format(axis="y",style="sci",scilimits=(-2,2));x.legend(frameon=False)
    fig.suptitle("Short Test 65    ·    Re = 882    ·    14D 远场",fontsize=15,fontweight="bold")
    fig.savefig(HERE/"evidence/later.png",dpi=165)
