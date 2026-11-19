"""在同一物理截面上绘出原生壁压、热流及共同方位的离散差值。"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(result, path):
    cone=result["nativeSections"][1]
    a=np.asarray(cone["coarseAzimuthDeg"])
    b=np.asarray(cone["fineAzimuthDeg"])
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
    fig,axes=plt.subplots(2,2,figsize=(9.6,5.9),sharex="col",gridspec_kw={"height_ratios":[2,1]})
    for column,(key,scale,label) in enumerate((("pressureGauge",1000,"相对环境壁压 / mPa"),
                                               ("heatFlux",1,"入壁热流 / (W/m²)"))):
        coarse=np.asarray(cone[key]["coarse"])*scale
        fine=np.asarray(cone[key]["fine"])*scale
        top,bottom=axes[:,column]
        top.plot(a,coarse,"o--",color="#235b87",ms=4.2,mfc="white",label="半周 12 单元")
        top.plot(b,fine,".-",color="#a05125",ms=5,lw=1.4,label="半周 24 单元")
        top.set_ylabel(label)
        top.legend(loc="best",frameon=False)
        bottom.axhline(0,color="#8c9194",lw=.8)
        bottom.plot(a,fine[::2]-coarse,"o-",ms=3.6,color="#385a66")
        bottom.set_ylabel("24 减 12 / "+("mPa" if column==0 else "W/m²"))
        bottom.set_xlabel("方位角 / °")
        for ax in (top,bottom):
            ax.set_xlim(-94,94)
            ax.set_xticks(np.arange(-90,91,30))
            ax.grid(axis="y",color="#e2e5e6",lw=.7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    fig.suptitle("CM-011 锥面周向分布   ·   X = %.6f mm"%(cone["xM"]*1000),fontsize=14)
    fig.text(.5,.925,"来流 0.1 m/s   ·   攻角 +10°   ·   −90° 为迎流侧，+90° 为背流侧",ha="center",color="#455158")
    fig.tight_layout(rect=(.015,.01,.99,.90),h_pad=1.6,w_pad=2.6)
    fig.savefig(path,dpi=165,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    plot(json.loads(args.comparison.read_bytes()),args.output)
