"""绘制同一母平面的两套原生网格和最差标量节点绕支管轴投影的位置。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

HERE=Path(__file__).resolve().parent
if __name__=="__main__":
    data=json.loads((HERE.parents[2]/"cases/short_sphere/supportedGrid65.json").read_bytes())
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,1,figsize=(9.3,5.8),constrained_layout=True,sharex=True,sharey=True)
    for ax,name,title in zip(axes,("hybrid","structured"),("混合外区","结构外区")):
        d=data[name];s=d["section"];xy=np.array(s["pointsD"]);edges=np.array(s["edges"])-1
        ax.add_collection(LineCollection(xy[edges],colors="#667d87",linewidths=.45))
        frame=d["stages"][-1]["frames"][-1];x,y,z=frame["worstXYZD"]
        ax.plot(y,np.hypot(x,z),"o",color="#a3451c",mec="white",mew=1,ms=8,label="终帧最低值节点的径向投影")
        ax.set(xlim=(-1.8,2.1),ylim=(-.025,1.15),aspect="equal",ylabel="距支管轴 r/D")
        ax.set_title(f"{title}    ·    {d['geometry']['cells']:,} 个体单元    ·    末帧 {frame['below']} 个节点低于下限",fontsize=11)
        ax.legend(loc="upper right",framealpha=.95,fontsize=8)
    axes[-1].set_xlabel("沿支管轴 Y/D")
    fig.savefig(HERE/"evidence/mesh.png",dpi=160)
    plt.close(fig)
