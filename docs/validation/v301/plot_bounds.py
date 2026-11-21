"""用保存帧绘制温界偏差与末帧低温节点的轴向、径向位置。"""
import argparse
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from transport_case import read
from thermistorLeadMesh import meridian
from thermistorBead import D


def plot(folder):
    report=read(folder/"frames.json")
    with gzip.open(folder/"out_of_bounds.csv.gz","rt",encoding="ascii") as f:
        nodes=np.atleast_1d(np.genfromtxt(f,names=True,delimiter=","))
    frames=report["frames"];low=nodes[nodes["candidate_theta"]<-1e-12]
    assert len(low)==frames[-1]["below"] and len(low)>0
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,(history,position)=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    fig.suptitle("3 m/s 带线热敏珠    温差标量下界检查",fontsize=14)
    for label,pts,color in (
        ("MUSCL + Barth–Jespersen",[(f["step"],-f["min"]["theta"]) for f in frames],"#a94738"),
        ("原一阶，同迭代",[(f["step"],-f["sameIteration"]["baseThetaRange"][0]) for f in frames if "sameIteration" in f],"#326b80")):
        pts=[(s,x) for s,x in pts if x>0]
        if pts:history.semilogy(*zip(*pts),marker="o",color=color,label=label)
    history.axhline(1e-12,ls="--",lw=1,color="#555",label="允许下界偏差")
    history.set(xlabel="原生稳态迭代编号",ylabel="低于零的最大幅度  |min θ|",title="完整保存帧")
    history.legend(frameon=False,fontsize=9);history.grid(axis="y",alpha=.18)
    # 同一子午坐标在各方位可投影重合，颜色给出该节点实际下界偏差的数量级。
    dots=position.scatter(low["x_m"]*1000,np.hypot(low["y_m"],low["z_m"])*1000,
                          c=np.log10(-low["candidate_theta"]),s=11,cmap="viridis")
    wall,_=meridian();position.plot(wall[:,0]*D*1000,wall[:,1]*D*1000,color="#252525",lw=1.4,label="试件母线")
    position.set(xlabel="X  [mm]（+X 为迎流端）",ylabel="到轴线的径向距离  [mm]",title=f"终帧低于容差的 {len(low)} 个节点")
    position.legend(frameon=False,fontsize=9)
    fig.colorbar(dots,ax=position,label=r"$\log_{10}|\theta|$",shrink=.8)
    fig.savefig(folder/"bounds.png",dpi=155);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("folder",type=Path)
    plot(p.parse_args().folder)
