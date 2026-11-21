"""由带线原三角面绘制网距引起的热流变化，并按论文面积比较四个整珠平均量。"""
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from mesh_audit import HERE,ROOT,read
from compare_meshes import table


def plot():
    data=read(HERE/"evidence/comparison.json")
    cases=ROOT/"cases/thermistor_b"
    with zipfile.ZipFile(cases/"leadsMidbody3.zip") as z:
        body=json.loads(z.read("body.json"));fine=table(z.read("wall.csv"))
    with zipfile.ZipFile(cases/"leads3.zip") as z:base=table(z.read("wall.csv"))
    tri=np.asarray(body["V"])[np.asarray(body["F"])-1]*1000
    normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    show=normal[:,2]>0
    delta=(base["heatFluxWm2"]-fine["heatFluxWm2"])/1000
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False,
                         "axes.spines.right":False,"axes.spines.top":False})
    fig,(model,avg)=plt.subplots(1,2,figsize=(10,4.3),gridspec_kw={"width_ratios":[1.35,1]},constrained_layout=True)
    surface=PolyCollection(tri[show,:,:2],array=delta[show],cmap="RdBu_r",edgecolors="face",linewidths=.12)
    bound=max(abs(delta));surface.set_clim(-bound,bound)
    model.add_collection(surface);model.autoscale_view();model.set_aspect("equal")
    model.set(xlabel="X  [mm]   ·   +X 为迎流端",ylabel="Y  [mm]",title="同一带线试件的网距影响")
    fig.colorbar(surface,ax=model,orientation="horizontal",shrink=.84,pad=.17,label="向外热流变化，加密后减原网  [kW/m²]")
    y=np.arange(4);vals=[data["cases"][k]["meanHWm2K"] for k in ("baseSmooth","fineSmooth","baseWired","fineWired")]
    avg.barh(y,vals,height=.5,color=["#64788a","#278571","#64788a","#278571"])
    exp=read(HERE.parent/"v282/evidence/comparison.json")["experiment"]
    avg.axvspan(*exp["graphReadRangeWm2K"],color="#c39441",alpha=.20)
    avg.axvline(exp["hWm2K"],color="#845e24",ls="--",lw=1.1)
    for i,v in enumerate(vals):avg.text(v+4,i,f"{v:.1f}",va="center")
    avg.text(exp["hWm2K"],3.62,"实验图读 218.4",ha="center",fontsize=9)
    avg.set(yticks=y,yticklabels=["光滑原网","光滑加密","带线原网","带线加密"],xlim=(0,375),ylim=(4,-.6),
            xlabel="平均换热系数  [W/(m²·K)]",title="珠部平均换热")
    fig.suptitle("B 型热敏珠    3 m/s 空气  ·  给定壁温 305.15 K",fontsize=14)
    fig.savefig(HERE/"evidence/wired_heat.png",dpi=155);plt.close(fig)


if __name__=="__main__":
    plot()
