"""显示同一物理位置的网边、最宽区间以及 16 条锁定和 32 条独立测线。"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from flip_bound import DEFAULT_BASE,read,meridian,load_mesh,np
from thermistorBead import D


def plot(base,candidate,audit,destination):
    rows=read(audit/"rays.json")["rows"];hold=[r for r in rows if r["fraction"]!=.5]
    worst=max(hold,key=lambda r:r["widthRatio"])
    strip=worst["strip"];wall=meridian()[0];t=wall[strip+1]-wall[strip];origin=wall[strip]+worst["fraction"]*t;t/=np.linalg.norm(t)
    n=np.array([-t[1],t[0]]);unit=D*1e6
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":9,"axes.unicode_minus":False})
    fig,(ax1,ax2,ax3)=plt.subplots(1,3,figsize=(11,3.9),gridspec_kw={"width_ratios":[.9,.9,1.25]})
    fig.subplots_adjust(left=.07,right=.98,bottom=.25,top=.76,wspace=.36)
    for ax,folder,key,title in ((ax1,base,"baseCutsD","第三档原网"),(ax2,candidate,"newCutsD","169 处换边")):
        _,xyz,marks=load_mesh(folder/"mesh.su2")
        edges=np.array(sorted({tuple(sorted((i,j))) for f in marks["sym0"] for i,j in zip(f,f[1:]+f[:1])}))
        pts=xyz[edges,:2]-origin;local=unit*np.stack((pts@t,pts@n),axis=2)
        keep=(local[:,:,0].min(axis=1)<20)&(local[:,:,0].max(axis=1)>-20)&(local[:,:,1].min(axis=1)<50)&(local[:,:,1].max(axis=1)>15)
        ax.add_collection(LineCollection(local[keep],colors="#678093",linewidths=.9));ax.axvline(0,color="#434c50",ls="--",lw=.8)
        v=np.array(worst[key])*unit;j=np.argmax(np.diff(v))
        ax.plot([0,0],v[j:j+2],"o-",lw=3,ms=4,color="#a94a27")
        ax.set(xlim=(-20,20),ylim=(15,50),aspect="equal",title=f"{title}　{v[j+1]-v[j]:.2f} µm",xlabel="局部切向 [µm]")
    ax1.set_ylabel("距壁法向坐标 [µm]")
    ax3.axhspan(-25,25,color="#e6eeee");ax3.axhline(25,color="#56706f",ls="--",lw=.8)
    for frac,marker,color,label in ((.25,"o","#185879","1/4"),(.5,".","#707070","1/2（锁定）"),(.75,"^","#9b522e","3/4")):
        r=[v for v in rows if v["fraction"]==frac]
        ax3.scatter([v["strip"]+frac for v in r],[100*(v["widthRatio"]-1) for v in r],s=16,marker=marker,color=color,label=label)
    ax3.set(title="逐线最宽间距变化",xlabel="条带编号 + 段内位置",ylabel="相对第三档 [%]",ylim=(-30,110))
    ax3.legend(loc="upper center",bbox_to_anchor=(.5,-.26),ncol=3,fontsize=8,frameon=False)
    ax3.set_xticks([54,58,62,66,70]);fig.suptitle("热敏珠前中段　·　65 号条带四分点与全带间距核对",fontsize=13)
    fig.savefig(destination,dpi=160,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("candidate","audit","destination"):p.add_argument(k,type=Path)
    p.add_argument("--base",type=Path,default=DEFAULT_BASE);a=p.parse_args();plot(a.base,a.candidate,a.audit,a.destination)
