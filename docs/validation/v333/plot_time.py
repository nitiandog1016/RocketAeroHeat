"""绘制同刻热流、步长差以及功率之差与空气储热率的对照。"""
import json
from pathlib import Path
import sys
from zipfile import ZipFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
HERE=Path(__file__).resolve().parent;DATA=HERE.parents[2]/"cases/short_sphere"
sys.path.insert(0,str(HERE.parent/"v329"))
from plot_surface import view

if __name__=="__main__":
    p=json.loads((DATA/"supportedTime65.json").read_bytes())
    with ZipFile(DATA/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    q=[]
    for branch in ("Full","Half"):
        with ZipFile(DATA/("supportedTime65"+branch+".zip")) as z:
            q.append(-np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")["fluxInWm2"])
    v=np.array(g["V"])*1000;f=np.array(g["F"])-1;delta=q[1]-q[0]
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(10.2,6.7),constrained_layout=True);grid=fig.add_gridspec(3,2,height_ratios=[1,.10,1.12])
    a=view(fig.add_subplot(grid[0,0]),v,f,q[1]/1000,[-.9,.22,.3],"两半步 · 背风侧",plt.Normalize(0,np.ceil(q[1].max()/1000)))
    scale=np.max(abs(delta));b=view(fig.add_subplot(grid[0,1]),v,f,delta,[-.9,.22,.3],"同一原面的向外热流差",plt.Normalize(-scale,scale))
    b.set_cmap("RdBu_r");fig.colorbar(a,cax=fig.add_subplot(grid[1,0]),orientation="horizontal").set_label("向外放热通量 (kW/m²)")
    fig.colorbar(b,cax=fig.add_subplot(grid[1,1]),orientation="horizontal").set_label("两半步减整步 (W/m²)")
    ax=fig.add_subplot(grid[2,0]);x=p["angles"]
    ax.plot(x,p["fullNu"],"--o",ms=5,label="整步",color="#52847f");ax.plot(x,p["halfNu"],"-+",ms=7,label="两半步",color="#253d58")
    ax.scatter(x,p["experimentNu"],s=22,marker="D",color="#a8502f",label="Short 表 III")
    ax.set(xlabel="从迎风点沿 +Z 测线的极角 (°)",ylabel="当地 Nu",xticks=[0,45,90,135,180]);ax.legend(frameon=False)
    ax=fig.add_subplot(grid[2,1]);es=[p[k] for k in ("fullEnergy","firstHalfEnergy","halfEnergy")];xs=np.arange(3)
    ax.bar(xs-.17,[1e6*(e["wallW"]-e["farW"]) for e in es],.33,label="壁放热减远场输出",color="#64867f")
    ax.bar(xs+.17,[1e6*e["storageW"] for e in es],.33,label="由前后温差场积分的储热率",color="#304b65")
    ax.set_xticks(xs);ax.set_xticklabels(["整步","第一半步","第二半步"]);ax.set_ylabel("功率 (µW)")
    ax.set_ylim(0,1.3e6*max(e["storageW"] for e in es));ax.legend(loc="upper center",frameon=False,fontsize=9)
    ax.grid(axis="y",alpha=.2);ax.set_axisbelow(True)
    fig.suptitle(f'Short Test 65   |   扰动后 {1000*p["elapsedS"]:.6f} ms   |   Re = 882')
    fig.savefig(HERE/"evidence/time.png",dpi=155);plt.close(fig)
