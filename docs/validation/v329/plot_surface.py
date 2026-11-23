"""从发布包绘制两方向的同色标热流云图，与七个银球实测位置并列。"""
import io
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

HERE=Path(__file__).resolve().parent
DATA=HERE.parents[2]/"cases/short_sphere"


def view(ax,vertices,faces,values,direction,title,norm):
    towards=np.array(direction,dtype=float);towards/=np.linalg.norm(towards)
    horizontal=np.array([0.,1.,0.]);horizontal-=towards[1]*towards;horizontal/=np.linalg.norm(horizontal)
    vertical=np.cross(towards,horizontal)
    if vertical[2]<0:vertical=-vertical
    tri=vertices[faces]
    visible=(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])@towards)>0
    ids=np.flatnonzero(visible);ids=ids[np.argsort(tri[ids].mean(axis=1)@towards)]
    xy=(vertices@np.stack([horizontal,vertical],axis=1))[faces[ids]]
    cloud=PolyCollection(xy,cmap="viridis",norm=norm,edgecolors="face",linewidths=.14)
    cloud.set_array(values[ids]);ax.add_collection(cloud)
    ax.set(xlim=(-20,20),ylim=(-8.5,8.5),aspect="equal");ax.axis("off");ax.set_title(title,fontsize=12)
    ax.plot([-17,-12],[-7,-7],color="#434b50",lw=1.4);ax.text(-14.5,-6.6,"5 mm",ha="center",fontsize=8)
    return cloud


if __name__=="__main__":
    with ZipFile(DATA/"supported65.zip") as z:geom=json.loads(z.read("geometry.json"))
    with ZipFile(DATA/"supportedSurface65.zip") as z:
        w=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True)
        stations=np.genfromtxt(io.BytesIO(z.read("stations.csv")),delimiter=",",names=True)
        s=json.loads(z.read("summary.json"))
    v=np.array(geom["V"])*1000;f=np.array(geom["F"])-1
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(10.5,6.2),layout="constrained")
    grid=fig.add_gridspec(3,2,height_ratios=[1,.11,1])
    norm=plt.Normalize(0,7.5);q=-w["fluxInWm2"]/1000
    view(fig.add_subplot(grid[0,0]),v,f,q,[.9,.22,.3],"迎风侧 · 从 +X 侧观察",norm)
    cloud=view(fig.add_subplot(grid[0,1]),v,f,q,[-.9,.22,.3],"背风侧 · 从 −X 侧观察",norm)
    cb=fig.colorbar(cloud,cax=fig.add_subplot(grid[1,:]),orientation="horizontal")
    cb.set_label("向外放热通量 (kW/m²)");cb.set_ticks([0,1.5,3,4.5,6,7.5])
    ax=fig.add_subplot(grid[2,0]);ax.plot(stations["angleDeg"],stations["computedNu"],"-o",label="计算",color="#305b76",ms=4)
    ax.scatter(stations["angleDeg"],stations["reportedNu"],label="Short 表 III",marker="D",s=26,color="#ab5029")
    ax.set(xlabel="+Z 测线极角 (°)",ylabel="局部 Nu",xticks=[0,45,90,135,180],xlim=(-4,180))
    ax.legend(frameon=False);ax.grid(axis="y",alpha=.18)
    ax=fig.add_subplot(grid[2,1]);vals=[r["outW"] for r in s["regions"]];y=np.arange(5)
    ax.barh(y,vals,color="#517481");ax.set_yticks(y,["外露银球","+Y 管壁","−Y 管壁","+Y 端盖","−Y 端盖"])
    ax.invert_yaxis();ax.set(xlabel="分区放热 (W)",xlim=(0,1.37))
    for i,value in enumerate(vals):ax.text(value+.02,i,f"{value:.3f}",va="center",fontsize=9)
    fig.suptitle("Short Test 65   |   名义双支管银球   |   Re = 882",fontsize=14)
    out=HERE/"evidence";out.mkdir(exist_ok=True)
    fig.savefig(out/"surface.png",dpi=155);plt.close(fig)
