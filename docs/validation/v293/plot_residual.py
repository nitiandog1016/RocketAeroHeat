"""绘制真实内迭代历史和保存终点的逐层残差份额，用于数值诊断报告。"""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


def draw(history,summary,target):
    raw=gzip.decompress(history.read_bytes()).decode("ascii")
    stream=io.StringIO(raw);keys=[s.strip() for s in next(csv.reader(stream,skipinitialspace=True))]
    values=np.loadtxt(stream,delimiter=",");s=json.loads(summary.read_bytes())
    assert len(values)==s["innerRows"] and s["limitingEquation"]=="rms[h]"
    names=["rms[P]","rms[U]","rms[V]","rms[W]","rms[h]","rms[rho*Y_0]"]
    plt.rcParams["font.family"]=font_manager.FontProperties(fname=r"C:\Windows\Fonts\msyh.ttc").get_name()
    fig,(ax,bx)=plt.subplots(1,2,figsize=(10.5,4),gridspec_kw={"width_ratios":[1.27,1]},constrained_layout=True)
    colors=["#286480","#80a5b5","#649273","#9997aa","#ab4c29","#737a57"]
    for i,(name,label) in enumerate(zip(names,["连续性","速度 U","速度 V","速度 W","焓","温差标量"])):
        ax.plot(values[:,keys.index("Inner_Iter")]+1,values[:,keys.index(name)],color=colors[i],label=label,lw=2 if i==4 else 1.2)
    ax.axhline(-8,color="#252a30",ls="--",lw=1,label="验收门槛")
    ax.axvspan(len(values)-19,len(values),color="#d3d6da",alpha=.45,zorder=0)
    ax.set(xlabel="内迭代次数",ylabel=r"原生残差 $\log_{10}(\mathrm{RMS})$",title="六项残差与末 20 次验收窗口",xlim=(1,len(values)))
    ax.legend(frameon=False,ncol=3,fontsize=8)
    portion=100*np.asarray(s["limitingResidualByLayer"]);idx=np.arange(len(portion))
    bx.bar(idx,portion,color=np.where(idx<=16,"#286480","#a5adb3"),width=.86)
    bx.set(xlabel="壁到远场的节点层号（壁面为 0）",ylabel="焓残差平方和占比 / %",title="终点剩余误差的径向分布",xlim=(-1,57))
    bx.text(.98,.96,f"加密带（层 0–16）合计 {portion[:17].sum():.1f}%",transform=bx.transAxes,ha="right",va="top",fontsize=9)
    for a in (ax,bx):
        a.spines["right"].set_visible(False);a.spines["top"].set_visible(False);a.grid(axis="y",alpha=.18);a.set_axisbelow(True)
    fig.savefig(target,dpi=160);plt.close(fig)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("target",type=Path)
    out=cli.parse_args().target;assert not out.exists()
    evidence=Path(__file__).resolve().parent/"evidence"
    draw(evidence/"first_history.csv.gz",evidence/"first_state.json",out)
