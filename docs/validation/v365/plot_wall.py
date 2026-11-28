"""按圆柱原面轴向环带画热流与壁剪切，标出末帧两档同时出现近壁反向的位置。"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DOC=Path(__file__).resolve().parent


def draw():
    raw=json.loads((DOC/"research/surface.json").read_bytes())
    geometry=DOC.parents[2]/"cases/simple_launcher/launcher-m.stl"
    dtype=np.dtype([("n","<f4",3),("v","<f4",(3,3)),("att","<u2")])
    sel=slice(2016,6912)
    tri=np.fromfile(geometry,dtype=dtype,offset=84)["v"][sel].astype(float)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    ends,band=np.unique(np.sort(tri[:,:,0],axis=1)[:,[0,-1]],axis=0,return_inverse=True)
    x=ends.mean(axis=1)
    groups=np.array(raw["faceGroup"])[sel]-1
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False,"font.size":10})
    fig,axes=plt.subplots(2,1,figsize=(8,5),sharex=True,layout="constrained")
    for key,color,style,label in (("normal4","#354f6c","--","首层 4 µm"),("normal2","#a44b18","-","首层 2 µm")):
        values=np.array(raw[key])[groups]
        for ax,column,scale in zip(axes,(1,2),(1e-3,1)):
            average=np.bincount(band,area*values[:,column])/np.bincount(band,area)
            ax.plot(x,average*scale,style,color=color,lw=1.6,label=label)
    for ax in axes:
        if raw["commonReverseCylinderFaces"]:
            ax.axvspan(*raw["commonReverseXRangeM"],color="#666666",alpha=.10)
        ax.axhline(0,lw=.8,color="#838383")
        ax.grid(axis="y",color="#dddddd",lw=.6)
        ax.spines[["top","right"]].set_visible(False)
    axes[0].set_title("M0.30 圆柱段 · 全壁层流试算末帧，稳态未收敛",loc="left")
    axes[0].set_ylabel("入壁热流  (kW/m²)")
    axes[0].legend(frameon=False,ncol=2)
    axes[1].set_ylabel("沿来流壁剪切  (Pa)")
    axes[1].set_xlabel("模型 X  (m)      来流从头罩侧流向收缩肩 →")
    axes[1].set_xlim(ends.max(),ends.min())
    fig.savefig(DOC/"cylinder.png",dpi=155)
    plt.close(fig)


if __name__=="__main__":
    draw()
