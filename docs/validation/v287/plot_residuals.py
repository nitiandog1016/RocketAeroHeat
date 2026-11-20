"""从实际保存的内迭代记录绘制首步残差。合并同量级速度曲线时取较大的分量。"""
import csv
import gzip
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(folder):
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","Noto Sans CJK SC","DejaVu Sans"]
    with gzip.open(folder/"inner_history.csv.gz","rt",encoding="ascii") as f:
        keys=[s.strip() for s in next(csv.reader(f,skipinitialspace=True))]
        data=np.loadtxt(f,delimiter=",")
    fig,ax=plt.subplots(figsize=(7.4,3.6),constrained_layout=True)
    n=data[:,keys.index("Inner_Iter")]
    for label,fields,color,style in [
        ("焓",["h"],"#A94F22","-"),("压力",["P"],"#263B55","-"),
        ("轴向速度",["U"],"#467F9E","--"),("横向速度（两分量较大值）",["V","W"],"#5D7856","-."),
        ("温差标量",["rho*Y_0"],"#6E547C","-")]:
        y=data[:,[keys.index(f"rms[{v}]") for v in fields]].max(axis=1)
        ax.plot(n,y,style,label=label,color=color,linewidth=1.6)
    ax.axhline(-8,color=".35",linewidth=.9,linestyle=":",label="验收门槛")
    ax.axvspan(n[-20],n[-1],color=".5",alpha=.12,label="末 20 次验收窗")
    ax.set(xlabel="首个物理步的内迭代编号",ylabel="原生残差 log$_{10}$ RMS",xlim=(0,n[-1]))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.legend(fontsize=8,loc="upper right",ncol=2,frameon=False)
    fig.savefig(folder/"residuals.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    plot(Path(__file__).resolve().parent/"evidence")
