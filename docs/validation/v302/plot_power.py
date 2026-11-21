"""从工作台交付的原面表绘制入壁热流，并标出独立温差计算对功率反演的核对点。"""
import csv
import io
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from half_case import HERE,ROOT,read


def plot():
    e=HERE/"evidence";d=read(e/"workbench.json");half=read(e/"half.json");t=d["ui"]["thermal"]
    with zipfile.ZipFile(e/"power4mW.zip") as z:
        rows=list(csv.DictReader(io.StringIO(z.read("faces.csv").decode("utf-8-sig"))))
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/leadsMidbody3.zip") as z:b=json.loads(z.read("body.json"))
    tri=np.asarray(b["V"])[np.asarray(b["F"])-1]*1000
    front=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])[:,2]>0
    q=np.array([float(r["into_wall_heat_flux_W_m2"]) for r in rows])/1000
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,(model,law)=plt.subplots(1,2,figsize=(10,4.2),constrained_layout=True)
    fig.suptitle("带线 B 型热敏珠    3 m/s 空气  ·  整件净对流功率",fontsize=14)
    surface=PolyCollection(tri[front,:,:2],array=q[front],cmap="viridis",edgecolors="face",linewidths=.13)
    surface.set_clim(q.min(),0);model.add_collection(surface);model.autoscale_view();model.set_aspect("equal")
    model.set(xlabel="X  [mm]（+X 为迎流端）",ylabel="Y  [mm]",title="4 mW 工况的全表面解")
    fig.colorbar(surface,ax=model,orientation="horizontal",pad=.16,shrink=.86,label="入壁热流  [kW/m²]（负值向空气放热）")
    law.plot([0,1000*t["referenceReleaseW"]],[t["ambientTemperatureK"],t["referenceWallTemperatureK"]],
             color="#4f6778",label="按整件空气侧热导反演")
    law.scatter(4,t["wallTemperatureK"],s=42,color="#14796a",zorder=3,label="工作台 4 mW")
    law.scatter(1000*half["wallReleaseW"],299.15,s=70,facecolors="none",edgecolors="#ae5826",lw=2,zorder=4,label="独立冷态，给定 299.15 K")
    law.annotate("297.563 K",(4,t["wallTemperatureK"]),xytext=(8,-18),textcoords="offset points")
    law.set(xlabel="整件净对流功率  [mW]",ylabel="整件等温壁温  [K]",title="功率与稳态壁温")
    law.grid(axis="y",alpha=.2);law.legend(frameon=False,fontsize=9)
    fig.savefig(e/"power_check.png",dpi=155);plt.close(fig)


if __name__=="__main__":plot()
