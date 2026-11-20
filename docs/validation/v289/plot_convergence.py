"""分别显示定常松弛和首个物理步的残差，直接使用随版本发布的原生历史。"""
import csv
import gzip
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def history(path):
    with gzip.open(path,"rt",encoding="ascii") as stream:
        keys=[k.strip() for k in next(csv.reader(stream,skipinitialspace=True))]
        values=np.loadtxt(stream,delimiter=",",ndmin=2)
    return dict(zip(keys,values.T))


def plot():
    folder=Path(__file__).resolve().parent/"evidence"
    prev=folder.parents[1]/"v288/evidence/first_history.csv.gz"
    steady=history(folder/"steady_history.csv.gz")
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","Noto Sans CJK SC","DejaVu Sans"]
    has_step=(folder/"step_history.csv.gz").is_file()
    fig,axes=plt.subplots(1,2 if has_step else 1,figsize=(8.6 if has_step else 5.3,3.3),
                          squeeze=False,constrained_layout=True)
    a=axes[0,0]
    for key,label,color in [("rms[h]","焓","#A34725"),("rms[P]","压力","#2F6480"),
                            ("rms[rho*Y_0]","强温差标量","#5A7161")]:
        a.plot(steady["Inner_Iter"],steady[key],lw=1.4,color=color,label=label)
    a.set(title="无扰动全周基流",xlabel="定常松弛迭代",ylabel="原生残差 log$_{10}$ RMS")
    if has_step:
        extra=folder/"window_history.csv.gz"
        now=history(extra if extra.is_file() else folder/"step_history.csv.gz")
        old=history(prev)
        b=axes[0,1]
        b.plot(old["Inner_Iter"],old["rms[h]"],"--",color="#7D8487",label="v288 焓",lw=1.3)
        b.plot(now["Inner_Iter"],now["rms[h]"],color="#A34725",label="全周松弛后 · 焓",lw=1.6)
        b.plot(now["Inner_Iter"],now["rms[P]"],color="#2F6480",label="全周松弛后 · 压力",lw=1.3)
        b.axvspan(now["Inner_Iter"][-20],now["Inner_Iter"][-1],color=".5",alpha=.12,label="末 20 次")
        if extra.is_file():
            b.axvline(179,color=".5",ls="--",lw=.8,label="原 180 次终点")
        b.set(title="加入同一扰动后的首步",xlabel="首个物理步内迭代",ylim=(-11.4,-.7))
    for ax in axes[0]:
        ax.axhline(-8,color=".4",linestyle=":",lw=1,label="验收门槛")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.legend(fontsize=7.5,frameon=False,ncol=2 if has_step and ax is b else 1,loc="upper right")
    fig.savefig(folder/"convergence.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    plot()
