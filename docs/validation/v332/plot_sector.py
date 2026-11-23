"""绘制同一视角、同一色标的两套热流云图，以及七站换热和反向面积。"""
import json
from pathlib import Path
import sys
from zipfile import ZipFile
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE=Path(__file__).resolve().parent;DATA=HERE.parents[2]/"cases/short_sphere"
sys.path.insert(0,str(HERE.parent/"v329"))
from plot_surface import view


if __name__=="__main__":
    p=json.loads((DATA/"supportedAngular65.json").read_bytes())
    with ZipFile(DATA/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    v=np.array(g["V"])*1000;f=np.array(g["F"])-1
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(9.8,6.4),constrained_layout=True)
    grid=fig.add_gridspec(3,2,height_ratios=[1,.09,1.25])
    limit=np.ceil(-min(p["sharedFluxRangeKWm2"])*2)/2;norm=plt.Normalize(0,limit)
    for j,name in enumerate(("supportedFine65.zip","supportedAngularSurface65.zip")):
        with ZipFile(DATA/name) as z:wall=np.genfromtxt(z.open("wall.csv"),delimiter=",",names=True)
        cloud=view(fig.add_subplot(grid[0,j]),v,f,-wall["fluxInWm2"]/1000,[-.88,.22,.36],str(64*2**j)+" 扇区 · 背风侧",norm)
    fig.colorbar(cloud,cax=fig.add_subplot(grid[1,:]),orientation="horizontal").set_label("向外热流 (kW/m²)")
    ax=fig.add_subplot(grid[2,0]);angles=p["anglesDeg"]
    ax.plot(angles,p["coarseNu"],"--o",color="#4c7775",ms=4,label="64 扇区")
    ax.plot(angles,p["fineNu"],"-+",color="#263e60",ms=7,label="128 扇区")
    ax.scatter(angles,p["measuredNu"],marker="D",s=22,color="#a24b29",label="Short 表 III")
    ax.set(xlabel="从迎风点沿 +Z 测线的极角 (°)",ylabel="当地 Nu",xticks=[0,45,90,135,180])
    ax.legend(frameon=False,fontsize=8);ax.grid(axis="y",alpha=.2)
    ax=fig.add_subplot(grid[2,1]);x=np.arange(3)
    for j,(k,label,color) in enumerate((("baseDirection","64 扇区","#698b8b"),("direction","128 扇区","#304b6a"))):
        vals=[1e6*p[k]["regions"][str(n)]["reverse"]["areaM2"] for n in (1,2,3)]
        ax.bar(x+(j-.5)*.34,vals,width=.32,color=color,label=label)
    ax.set_xticks(x);ax.set_xticklabels(["外露球面","+Y 管壁","−Y 管壁"])
    ax.set_ylabel("局部反向面积 (mm²)");ax.yaxis.set_major_locator(MaxNLocator(5));ax.grid(axis="y",alpha=.2)
    ax.set_axisbelow(True);ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Short Test 65   |   十层近壁网的周向间隔二分   |   Re = 882",fontsize=14)
    fig.savefig(HERE/"evidence/sector.png",dpi=160);plt.close(fig)
