"""绘制原论文七个角度的分布、对应误差及两种实验功率口径。"""
import argparse
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("surface",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();rows=np.genfromtxt(a.surface/"stations.csv",delimiter=",",names=True)
    s=json.loads((a.surface/"summary.json").read_text("utf8"))
    plt.rcParams.update({"font.sans-serif":["Microsoft YaHei","SimHei","DejaVu Sans"],"axes.unicode_minus":False,"font.size":10})
    fig,(ax,err)=plt.subplots(2,1,figsize=(7.3,5.4),sharex=True,gridspec_kw={"height_ratios":[2,1]},constrained_layout=True)
    x=rows["angleDeg"]
    ax.plot(x,rows["computedNu"],"-o",color="#335872",label="光滑球 · SU2 稳态基流")
    ax.scatter(x,rows["reportedNu"],marker="D",s=36,color="#ac5f32",label="风洞测站 · Short 表 III",zorder=3)
    ax.set(ylabel="局部 Nu（表温物性）",ylim=(0,49))
    ax.legend(frameon=False,loc="upper right");ax.set_title("12.7 mm 银球局部换热    |    Test 65，Re = 882",loc="left")
    err.axhline(0,color=".45",lw=.8);err.vlines(x,0,rows["deltaNu"],color="#335872",lw=2)
    err.scatter(x,rows["deltaNu"],color="#335872",s=20)
    err.set(xlabel="距迎风驻点的极角 / °",ylabel="计算 − 实验 / Nu",xticks=[0,30,60,90,120,150,180],xlim=(-4,180),ylim=(-12,10))
    for axis in (ax,err):
        axis.grid(axis="y",alpha=.22);axis.spines["top"].set_visible(False);axis.spines["right"].set_visible(False)
    fig.savefig(a.out,dpi=165);plt.close(fig)
