"""由交付网格绘制装置侧视图和正 Y 接头处的原球面遮盖比例。"""
import io
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from scipy.spatial import ConvexHull

HERE=Path(__file__).resolve().parent
CASES=HERE.parents[2]/"cases/short_sphere"


def plot():
    with ZipFile(CASES/"supported65.zip") as z:
        mesh=json.loads(z.read("geometry.json"))
        parents=np.genfromtxt(io.BytesIO(z.read("parents.csv")),delimiter=",",names=True)
    with ZipFile(CASES/"surface65.zip") as z:
        ref=json.loads(z.read("body.json"))
    v=1000*np.array(mesh["V"]);f=np.array(mesh["F"])-1;part=np.array(mesh["regionID"])
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig,(a,b)=plt.subplots(1,2,figsize=(10,4.1),gridspec_kw={"width_ratios":[1.55,1]},constrained_layout=True)
    for k,color in [(2,"#91a9b2"),(3,"#91a9b2"),(1,"#dee5e9")]:
        uv=v[np.unique(f[part==k])][:,[1,0]]
        hull=uv[ConvexHull(uv).vertices]
        a.fill(*hull.T,fc=color,ec="#445968",lw=.85)
    for key in ("joinPositive","joinNegative"):
        loop=v[np.array(mesh[key])[:,0]-1]
        a.plot(loop[:,1],loop[:,0],color="#a65327",lw=1.5)
    a.annotate("来流",xy=(0,7.0),xytext=(0,10.6),ha="center",arrowprops={"arrowstyle":"->","color":"#47545c"})
    a.set(xlim=(-21,21),ylim=(-10,12),aspect="equal",title="名义外形与两处接缝",xlabel="Y / mm",ylabel="X / mm")
    old=(1000*np.array(ref["V"]))[np.array(ref["F"])-1]
    near=(old[:,:,1].min(axis=1)>0)&(np.abs(old[:,:,[0,2]]).min(axis=1)<1.9).all(axis=1)
    fraction=parents["coveredM2"]/parents["originalM2"]
    patch=PolyCollection(old[near][:,:,[0,2]],edgecolors="#939fa7",linewidths=.45,cmap="Blues")
    patch.set_array(fraction[near]*100);patch.set_clim(0,100);b.add_collection(patch)
    loop=v[np.array(mesh["joinPositive"])[:,0]-1];circle=np.vstack([loop,loop[0]])
    b.plot(circle[:,0],circle[:,2],color="#a65327",lw=1.7)
    b.set(xlim=(-1.9,1.9),ylim=(-1.9,1.9),aspect="equal",title="+Y 接头处的原球面",xlabel="X / mm",ylabel="Z / mm")
    fig.colorbar(patch,ax=b,shrink=.74,label="原面被遮盖面积 / %")
    fig.savefig(HERE/"evidence/geometry.png",dpi=160)
    plt.close(fig)


if __name__=="__main__":
    plot()
