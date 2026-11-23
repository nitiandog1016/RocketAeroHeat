"""用相同观察方向、色标和原面绘制两套解及其热流差。"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import json
from zipfile import ZipFile

from compare import DATA, HERE, load_surface
sys.path.insert(0,str(HERE.parent/"v329"))
from plot_surface import view


if __name__=="__main__":
    with ZipFile(DATA/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    p=json.loads((DATA/"supportedNearWall65.json").read_bytes())
    _,c,_=load_surface(DATA/"supportedSurface65.zip");_,f,_=load_surface(DATA/"supportedFine65.zip")
    vertices=np.array(g["V"])*1000;faces=np.array(g["F"])-1;coarse=-c["fluxInWm2"];fine=-f["fluxInWm2"]
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(10.5,6.8),layout="constrained")
    grid=fig.add_gridspec(4,2,height_ratios=[1,.1,1,.1])
    direction=[.9,.22,.3];norm=plt.Normalize(0,7.5)
    view(fig.add_subplot(grid[0,0]),vertices,faces,coarse/1000,direction,"五层 · 首层 25.4 µm",norm)
    same=view(fig.add_subplot(grid[0,1]),vertices,faces,fine/1000,direction,"十层 · 首层 12.7 µm",norm)
    bar=fig.colorbar(same,cax=fig.add_subplot(grid[1,:]),orientation="horizontal");bar.set_label("向外热流 (kW/m²) · 两网共用色标")
    delta=view(fig.add_subplot(grid[2,0]),vertices,faces,fine-coarse,direction,"同一原面的向外热流差",plt.Normalize(-140,140))
    delta.set_cmap("RdBu_r");bar=fig.colorbar(delta,cax=fig.add_subplot(grid[3,0]),orientation="horizontal")
    bar.set_label("十层减五层 (W/m²)")
    ax=fig.add_subplot(grid[2:,1]);changes=[100*r["relativePowerChange"] for r in p["regions"]]
    ax.barh(np.arange(5),changes,color=["#a54e30" if x>0 else "#466a80" for x in changes])
    ax.set_yticks(range(5),["外露银球","+Y 管壁","−Y 管壁","+Y 端盖","−Y 端盖"]);ax.invert_yaxis()
    ax.set(xlim=(-1.36,.48),xlabel="分区放热变化 (%)");ax.axvline(0,color="#34383c",lw=.7)
    for k,x in enumerate(changes):ax.text(x+(.025 if x>0 else -.025),k,f"{x:+.3f}",ha="left" if x>0 else "right",va="center",fontsize=8)
    fig.suptitle("Short Test 65   |   近壁间隔对半   |   Re = 882",fontsize=14)
    fig.savefig(HERE/"evidence/nearwall.png",dpi=155);plt.close(fig)
