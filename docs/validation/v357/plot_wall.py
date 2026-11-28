"""按原 STL 的轴向面带对同刻壁压及热流求面积平均，在同一坐标轴上绘出两套动量网。"""
import argparse
from pathlib import Path
from zipfile import ZipFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from refine_flow import ROOT,geometry


def plot(bundle,dest):
    case=ROOT/"cases/simple_launcher/wake"
    with ZipFile(case/"wake_heat_grid.zip") as z:
        old=np.loadtxt(z.open("fine_wall.csv"),delimiter=",",skiprows=1)
    with ZipFile(bundle) as z:
        new=np.loadtxt(z.open("wall.csv"),delimiter=",",skiprows=1)
    _,tri,_,_=geometry();x=tri[:,:,0]
    bands,ids=np.unique(np.c_[x.min(axis=1),x.max(axis=1)],axis=0,return_inverse=True)
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    np.testing.assert_allclose(area,new[:,1],rtol=0,atol=1e-15)
    np.testing.assert_array_equal(new[:,:2],old[:,:2])
    mass=np.bincount(ids,weights=area);loc=np.bincount(ids,weights=area*x.mean(axis=1))/mass
    base=(bands==0).all(axis=1)
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False})
    fig,axes=plt.subplots(2,1,figsize=(8.4,5.4),sharex=True,layout="constrained")
    for table,label,color,ls in ((old,"原动量网","#9a521b","--"),(new,"加密动量网","#1c5b83","-")):
        for ax,field,unit in zip(axes,((table[:,2]-101325)*1000,-table[:,3]*1000),
                ("相对环境壁压 / mPa","向外热流 / (mW/m²)")):
            mean=np.bincount(ids,weights=area*field)/mass
            ax.plot(loc[~base],mean[~base],c=color,ls=ls,lw=1.6,label=label)
            ax.scatter(loc[base],mean[base],s=23,facecolors="none",edgecolors=color,zorder=3)
            ax.set_ylabel(unit);ax.spines[["top","right"]].set_visible(False)
    axes[0].legend(frameon=False,ncol=2);axes[0].axhline(0,color="#777",lw=.5)
    axes[1].set_xlabel("原面带轴向位置 X / m（平底在 0，鼻端在右）")
    fig.suptitle("0.05 m/s 风洞 · 0.600 s · 同一细热网的原面带面积平均")
    fig.savefig(dest,dpi=155);plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle",type=Path);parser.add_argument("dest",type=Path)
    a=parser.parse_args();plot(a.bundle,a.dest)
