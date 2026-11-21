"""绘制真实珠体的两次网距细化热流差，以及同工况珠部平均换热与实验读数。"""
import io
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from audit_third_mesh import HERE, ROOT, read


def plot():
    proof=HERE/"evidence"
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/leadsThird3.zip") as z:
        body=json.loads(z.read("body.json"))
    tri=1000*np.asarray(body["V"])[np.asarray(body["F"])-1]
    front=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])[:,2]>0
    with zipfile.ZipFile(proof/"face_comparison.zip") as z:
        w=np.genfromtxt(io.BytesIO(z.read("givenWall.csv")),delimiter=",",names=True)
    q=np.column_stack([w[f"q{i}_into_W_m2"] for i in (1,2,3)])
    delta=np.diff(q,axis=1);bound=abs(delta).max()
    summary=read(proof/"comparison.json");h=[g["meanHWm2K"] for g in summary["grids"]]
    exp=read(HERE.parent/"v282/evidence/comparison.json")["experiment"]
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":9})
    fig,axes=plt.subplots(1,3,figsize=(11,3.7),gridspec_kw={"width_ratios":[1.1,1.1,.95]})
    fig.subplots_adjust(left=.065,right=.98,top=.77,bottom=.29,wspace=.48)
    for j,ax in enumerate(axes[:2]):
        coll=PolyCollection(tri[front,:,:2],array=delta[front,j],cmap="RdBu_r",edgecolors="face",linewidths=.12)
        coll.set_clim(-bound,bound);ax.add_collection(coll);ax.autoscale_view();ax.set_aspect("equal")
        ax.set(title=("第二档减第一档","第三档减第二档")[j],xlabel="X [mm]（+X 迎流）",ylabel="Y [mm]")
    fig.colorbar(coll,cax=fig.add_axes([.13,.16,.43,.035]),orientation="horizontal",label="入壁热流差 [W/m²]")
    ax=axes[-1];ax.set_axisbelow(True);ax.grid(axis="y",alpha=.25)
    ax.axhspan(*exp["graphReadRangeWm2K"],color="#e4d7cd",alpha=.65)
    ax.axhline(exp["hWm2K"],color="#8b593c",ls="--",label="实验图读值及读图范围")
    ax.plot([1,2,3],h,"o-",color="#205b70")
    for j,value in enumerate(h,1):ax.annotate(f"{value:.1f}",(j,value),xytext=(0,5),textcoords="offset points",ha="center")
    ax.set_xticks([1,2,3],["第一档","第二档","第三档"])
    ax.set(xlim=(.75,3.25),ylim=(min(*h,exp["hWm2K"])-12,max(h)+17),title="珠部含接头平均换热",ylabel="h [W/(m²·K)]")
    ax.legend(loc="upper center",bbox_to_anchor=(.5,-.18),frameon=False)
    fig.suptitle("带线 B 型热敏珠　·　3 m/s　·　给定壁温 305.15 K",fontsize=13)
    fig.savefig(proof/"grid_heat.png",dpi=160,bbox_inches="tight");plt.close(fig)


if __name__=="__main__":
    plot()
