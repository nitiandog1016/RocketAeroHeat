"""以原面热流、功率括区和七站实验比较说明银球的热平衡解。"""
import json
from pathlib import Path
import sys
from zipfile import ZipFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent;CASE=HERE.parents[2]/"cases/short_sphere"
sys.path.insert(0,str(HERE.parent/"v329"))
from plot_surface import view


def draw():
    proof=json.loads((HERE/"evidence/comparison.json").read_bytes());p=json.loads((CASE/"calorimetry65.json").read_bytes())
    with ZipFile(CASE/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    with ZipFile(CASE/"calorimetry65.zip") as z:
        w=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",");s=np.genfromtxt(z.open("stations.csv"),names=True,delimiter=",")
    navy="#245977";orange="#b45728"
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False,"axes.unicode_minus":False})
    fig=plt.figure(figsize=(11,7),constrained_layout=True);grid=fig.add_gridspec(2,2)
    ax=fig.add_subplot(grid[0,0]);out=-w["fluxInWm2"]/1000
    cloud=view(ax,np.array(g["V"])*1000,np.array(g["F"])-1,out,[.8,.25,.25],"向空气放热的原面热流密度",plt.Normalize(0,np.ceil(out.max())))
    fig.colorbar(cloud,ax=ax,orientation="horizontal",fraction=.09,pad=.02,label="向外热流 (kW/m²)")
    ax=fig.add_subplot(grid[0,1])
    for key,mark,label,col in (("low","o","定温原生解",navy),("high","o",None,navy),("balance","D","独立复算",orange)):
        r=p[key];ax.plot(r["temperatureK"],r["earlyRegionW"][0],mark,color=col,label=label)
    ax.plot([p[t]["temperatureK"] for t in ("low","high")],[p[t]["earlyRegionW"][0] for t in ("low","high")],color=navy,alpha=.55,zorder=0)
    ax.axhline(p["targetSphereW"],color="#4a4e52",ls="--",lw=1,label="量热功率")
    ax.set(xlabel="均匀球壁温度 (K)",ylabel="外露球面放热 (W)",title="球体功率与壁温");ax.legend(frameon=False)
    ax=fig.add_subplot(grid[1,0]);names=["外露球面","+Y 管","−Y 管","+Y 端盖","−Y 端盖"]
    vals=np.array(p["balance"]["earlyRegionW"]);bars=ax.barh(names[::-1],vals[::-1],color=[navy]*4+[orange])
    for bar,val in zip(bars,vals[::-1]):ax.text(val+.008,bar.get_y()+bar.get_height()/2,f"{val:.4f}",va="center",fontsize=9)
    ax.set(xlim=(0,vals.max()*1.18),xlabel="向外放热 (W)",title="球面量热与独立加热的支管")
    ax=fig.add_subplot(grid[1,1]);ax.plot(s["angleDeg"],s["computedNu"],"o-",color=navy,label="计算 · +Z 测线")
    ax.plot(s["angleDeg"],s["reportedNu"],"D--",color=orange,label="实验 · 表 III")
    ax.set(xlabel="测站极角 (°)",ylabel="局部 Nu",title="七个球面测站",xticks=[0,45,90,135,180]);ax.legend(frameon=False)
    fig.suptitle(f"Short Test 65    ·    Re = 882    ·    球面量热 {p['targetSphereW']:.6f} W    ·    反求壁温 {p['temperatureK']:.3f} K",
                 fontsize=15,fontweight="bold")
    assert proof["faces"]==len(w);fig.savefig(HERE/"evidence/calorimetry.png",dpi=160)


if __name__=="__main__":draw()
