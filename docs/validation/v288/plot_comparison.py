"""全周首步的验收窗和首次装配的焓残差位置，按两种不同取样阶段绘图。"""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(data,out):
    here=Path(__file__).resolve().parent
    old=json.loads((here.parent/"v287/evidence/first_step.json").read_bytes())
    new=json.loads((data/"audit01/result.json").read_bytes())
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","Noto Sans CJK SC","DejaVu Sans"]
    fig,(a,b)=plt.subplots(1,2,figsize=(8.2,3.25),constrained_layout=True)
    keys=list(old["last20MaxLogRMS"])
    y=np.arange(len(keys))
    for row,rec in enumerate((old,new)):
        a.scatter([rec["last20MaxLogRMS"][k] for k in keys],y+.12*(2*row-1),
                  color=["#8B9298","#205980"][row],marker=["o","s"][row],s=29,
                  label=["3 m/s 初猜","6 m/s 初猜"][row])
    a.axvline(-8,color=".35",ls=":",lw=1)
    a.set(yticks=y,yticklabels=["压力","轴向速度","横向速度 V","横向速度 W","焓","温差标量"],
          xlabel="末 20 次中最大的 log$_{10}$ RMS",title="全周首步验收窗")
    a.invert_yaxis()
    a.legend(frameon=False,fontsize=8)
    for name,label,color in [("probe_ray01","单射线","#8B9298"),("probe_mirror01","五射线镜像","#205980")]:
        with (data/name/"restart_00001.csv").open(encoding="ascii") as f:
            header=next(csv.reader(f))
            x,h=np.loadtxt(f,delimiter=",",usecols=[header.index("x"),header.index("Residual_Enthalpy")],unpack=True)
        weights=h*h
        edges=np.arange(-21,22)
        fractions,_=np.histogram(x,bins=edges,weights=weights/np.sum(weights)*100)
        assert np.isclose(fractions.sum(),100)
        b.step(edges[:-1]+.5,fractions,where="mid",color=color,label=label,lw=1.6)
    b.set(xlabel="轴向位置 x/D（下游为负）",ylabel="各轴向区间占残差平方和 / %",title="首次内迭代的焓残差分布")
    b.legend(frameon=False,fontsize=8)
    for axis in (a,b):
        axis.spines["right"].set_visible(False)
        axis.spines["top"].set_visible(False)
    fig.savefig(out,dpi=165)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    plot(args.data,args.output)
