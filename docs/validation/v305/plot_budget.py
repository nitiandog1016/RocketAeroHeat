"""按原三角面绘制分带放热差，并用同一组源数据展示四个珠带的分项功率。"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from common_layer import GRIDS, ROOT, read


def plot(path):
    b=read(ROOT/"cases/thermistor_b/layerBands.json");bands=b["bands"];d=[k["delta"] for k in bands]
    ids=np.asarray(b["faceBand"])-1
    body=read(GRIDS[1]/"body.json")
    triangles=1000*np.asarray(body["V"])[np.asarray(body["F"])-1]
    visible=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])[:,2]>0
    density=np.array([v["faceOutW"]/r["areaM2"] for r,v in zip(bands,d)])
    bound=np.max(abs(density))
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False})
    fig,(left,right)=plt.subplots(1,2,figsize=(10.6,4.5),gridspec_kw={"width_ratios":[1,1.15]})
    fig.subplots_adjust(left=.07,right=.98,bottom=.28,top=.80,wspace=.34)
    cloud=PolyCollection(triangles[visible,:,:2],array=density[ids[visible]],edgecolors="face",linewidths=.15,cmap="RdBu_r")
    cloud.set_clim(-bound,bound);left.add_collection(cloud);left.autoscale_view();left.set_aspect("equal")
    left.set(title="分带壁面放热变化",xlabel="X [mm]（+X 迎流）",ylabel="Y [mm]")
    fig.colorbar(cloud,cax=fig.add_axes([.11,.15,.3,.028]),orientation="horizontal",label="放热差 / 带面积 [W/m²]")
    y=np.arange(4);rows=d[3:7]
    for offset,key,label,color in ((-.22,"outerConvW","外侧对流","#98623e"),(0,"outerDiffW","外侧扩散","#376684"),(.22,"faceOutW","壁面放热","#242d39")):
        right.barh(y+offset,[r[key]*1e6 for r in rows],height=.19,label=label,color=color)
    right.set_yticks(y,["后极侧","后中段","前中段","前极侧"]);right.invert_yaxis()
    right.axvline(0,color=".55",lw=.8);right.set_axisbelow(True);right.grid(axis="x",alpha=.16)
    right.set(title="珠面各带的功率差",xlabel="第三档减第二档 [µW]")
    right.legend(ncol=3,loc="upper center",bbox_to_anchor=(.5,-.21),frameon=False,fontsize=9)
    fig.suptitle("带线 B 型热敏珠　·　3 m/s　·　给定壁温 305.15 K",fontsize=14)
    fig.savefig(path,dpi=160,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path);plot(p.parse_args().output)
