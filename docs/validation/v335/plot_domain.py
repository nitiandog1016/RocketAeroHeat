"""并列绘制七站绝对换热和外域差值，以及五个区域的有符号放热变化。"""
from pathlib import Path
import argparse
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(source,dest):
    data=json.loads(source.read_text(encoding="utf8"));s=data["half"]
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.spines.top":False,
                         "axes.spines.right":False,"axes.unicode_minus":False,"font.size":10})
    fig=plt.figure(figsize=(10,5.8),layout="constrained")
    grid=fig.add_gridspec(2,2,width_ratios=[1.3,1],height_ratios=[2,1])
    ax=fig.add_subplot(grid[0,0]);dx=fig.add_subplot(grid[1,0],sharex=ax);bar=fig.add_subplot(grid[:,1])
    x=data["angleDeg"];q=np.asarray(s["outerNu"]);q0=np.asarray(s["baseNu"])
    ax.plot(x,q0,"-",color="#465465",lw=2,label="14D")
    ax.plot(x,q,"o",mfc="white",mec="#216657",ms=5,label="20D")
    ax.scatter(x,data["measuredNu"],color="#934524",marker="x",label="Short 实验",zorder=5)
    ax.set(ylabel="局部 Nu",title="七站换热");ax.tick_params(labelbottom=False);ax.legend(frameon=False,ncol=3)
    dx.axhline(0,lw=.8,color="0.6");dx.plot(x,q-q0,"o-",color="#216657",ms=4)
    dx.set(xlabel="测站极角 (°)",ylabel="Nu 差 (20D − 14D)",xticks=[0,30,60,90,120,150,180])
    labels=["外露银球","+Y 支管","−Y 支管","+Y 端盖","−Y 端盖"]
    watts=[r["deltaW"]*1000 for r in s["regions"]]
    bar.axvline(0,lw=.8,color="0.5");bar.barh(labels,watts,color="#416c82",height=.56);bar.invert_yaxis()
    bar.set(xlabel="向外放热差 (mW)，20D − 14D",title="五区放热变化")
    bar.margins(x=.24);bar.grid(axis="x",alpha=.2)
    fig.suptitle(f"Short Test 65   ·   Re = 882   ·   扰动后 {1000*data['elapsedS']:.6f} ms   ·   两半步",fontsize=13)
    fig.savefig(dest,dpi=160);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("comparison",type=Path);p.add_argument("image",type=Path)
    a=p.parse_args();plot(a.comparison,a.image)
