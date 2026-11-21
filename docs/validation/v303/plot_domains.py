"""用同一物理色标显示两种热边界下的外域热流差。"""
import csv
import io
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from check_state import HERE,ROOT


def plot():
    e=HERE/"evidence"
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/outer30.zip") as z:body=json.loads(z.read("body.json"))
    xyz=np.asarray(body["V"])[np.asarray(body["F"])-1]*1000
    front=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])[:,2]>0
    diff=[]
    with zipfile.ZipFile(e/"face_comparison.zip") as z:
        for mode in ("givenWall","givenPower"):
            rows=csv.DictReader(io.StringIO(z.read(mode+".csv").decode("utf-8-sig")))
            diff.append(np.array([float(r["q30_into_W_m2"])-float(r["q20_into_W_m2"]) for r in rows]))
    limit=max(float(np.max(abs(d))) for d in diff)
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(1,2,figsize=(10,3.85))
    fig.subplots_adjust(left=.08,right=.98,top=.81,bottom=.26,wspace=.23)
    for ax,delta,title in zip(axes,diff,("给定壁温 305.15 K","给定整件净功率 4 mW")):
        surface=PolyCollection(xyz[front,:,:2],array=delta[front],cmap="RdBu_r",edgecolors="face",linewidths=.12)
        surface.set_clim(-limit,limit);ax.add_collection(surface);ax.autoscale_view();ax.set_aspect("equal")
        ax.set(title=title,xlabel="X  [mm]（+X 为迎流端）",ylabel="Y  [mm]")
    fig.suptitle("带线热敏珠    30D 减 20D 的入壁热流差",fontsize=15)
    fig.colorbar(surface,cax=fig.add_axes([.27,.085,.5,.04]),orientation="horizontal",label="入壁热流差  [W/m²]")
    fig.savefig(e/"domain_heat.png",dpi=155,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":plot()
