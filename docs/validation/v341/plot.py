"""将膜温自洽解与固定膜温解按相同参考面、五个区域和七个测站比较。"""
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
    p=json.loads((HERE/"evidence/comparison.json").read_bytes())
    with ZipFile(CASE/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    with ZipFile(CASE/"calorimetryFilm65.zip") as z:w=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")
    with ZipFile(CASE/"calorimetry65.zip") as z:old=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")
    np.testing.assert_array_equal(w["faceID"],old["faceID"])
    points=np.array(g["V"])*1000;faces=np.array(g["F"])-1;direction=[-.85,.2,.28]
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False,"axes.unicode_minus":False})
    fig=plt.figure(figsize=(10.8,7),constrained_layout=True);grid=fig.add_gridspec(3,2,height_ratios=[1,.10,.92])
    heat=w["fluxInWm2"]/1000
    cloud=view(fig.add_subplot(grid[0,0]),points,faces,heat,direction,"随壁温更新物性 · 入壁热流",plt.Normalize(np.floor(min(heat)),0))
    fig.colorbar(cloud,cax=fig.add_subplot(grid[1,0]),orientation="horizontal").set_label("入壁热流 (kW/m²)")
    delta=w["fluxInWm2"]-old["fluxInWm2"];limit=max(abs(delta))
    cloud=view(fig.add_subplot(grid[0,1]),points,faces,delta,direction,"更新膜温减固定膜温 · 同面差",plt.Normalize(-limit,limit))
    cloud.set_cmap("RdBu_r");fig.colorbar(cloud,cax=fig.add_subplot(grid[1,1]),orientation="horizontal").set_label("入壁热流变化 (W/m²)")
    ax=fig.add_subplot(grid[2,0]);names=["外露银球","+Y 管壁","−Y 管壁","+Y 端盖","−Y 端盖"]
    powers=np.asarray(p["regionDeltaW"])*1000
    ax.barh(names[::-1],powers[::-1],color=["#a45630" if q<0 else "#386780" for q in powers[::-1]])
    ax.axvline(0,c="#474e51",lw=.8);ax.set(xlabel="向外放热变化 (mW)",title="相同量热输入下的分区功率差")
    ax=fig.add_subplot(grid[2,1]);angles=[0,41.3,60.2,75.8,90,119.8,166.7]
    ax.plot(angles,p["experimentalNu"],"D--",color="#aa582e",label="实验");ax.plot(angles,p["nu"],"o-",color="#2b607b",label="膜温自洽")
    ax.set(xlabel="极角 (°)",ylabel="局部 Nu",title="试验与 +Z 测线",xticks=[0,45,90,135,180]);ax.legend(frameon=False)
    f=p["transport"][p["final"]]
    fig.suptitle(f"Short Test 65   ·   Re = 882   ·   平衡壁温 {f['wallK']:.3f} K   ·   膜温 {f['filmK']:.3f} K",fontsize=15,fontweight="bold")
    fig.savefig(HERE/"evidence/film.png",dpi=165)


if __name__=="__main__":draw()
