"""绘制三项同温差平均换热和原面放热通量的轴向分布。"""
import io
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from check_mesh import HERE,ROOT


def rings(ax,body,table,style,label,selected=None):
    tri=np.array(body["V"])[np.array(body["F"])-1]
    # 旧显示网的周向环数可变，以 0.05 mm 的固定轴向区间加权比较。
    x=tri[:,:,0].mean(axis=1)*1000
    selected=np.ones(len(x),bool) if selected is None else selected
    coord=[];heat=[]
    for lo in np.arange(-.575,.55,.05):
        mask=(x>=lo)&(x<lo+.05)&selected
        if mask.any():
            w=table["areaM2"][mask]
            coord.append(np.average(x[mask],weights=w))
            heat.append(np.average(-table["heatFluxWm2"][mask]/1000,weights=w))
    ax.plot(coord,heat,style,lw=2,label=label)


def main():
    cases=ROOT/"cases/thermistor_b"
    with zipfile.ZipFile(cases/"surface.zip") as z:body=json.loads(z.read("body.json"))
    with zipfile.ZipFile(cases/"surface_7p5.zip") as z:old=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    with zipfile.ZipFile(cases/"paired3.zip") as z:
        result=json.loads(z.read("comparison.json"));pair=json.loads(z.read("body.json"))
        smooth=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    with zipfile.ZipFile(cases/"leads3.zip") as z:
        wired=json.loads(z.read("body.json"));leads=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,(left,right)=plt.subplots(1,2,figsize=(11,4.75),gridspec_kw={"width_ratios":[.94,1.45]})
    fig.subplots_adjust(left=.14,right=.97,bottom=.19,top=.79,wspace=.42)
    fig.text(.055,.91,"B 型热敏珠：引线几何与数值离散的配对对照",fontsize=17,weight="bold")
    fig.text(.055,.84,"3 m/s 空气，293.15 K   ·   给定壁温 305.15 K   ·   轴向对称稳态",color="#4c5964")
    labels=["原网光滑","配对光滑","带线珠部"]
    colors=["#556d83","#178274","#b6592b"]
    values=[result["meanHWm2K"][k] for k in ("oldSmooth","pairedSmooth","wiredBeadAndJoints")]
    left.barh([2,1,0],values,height=.48,color=colors)
    for y,value in zip([2,1,0],values):left.text(value+5,y,f"{value:.1f}",va="center",fontsize=9)
    exp=result["experimentHWm2K"];left.axvline(exp,color="#222",ls=":",lw=1.4)
    left.text(exp-4,2.57,f"实验 {exp:.1f}",ha="right",fontsize=9)
    left.set(xlim=(0,382),ylim=(-.5,2.9),yticks=[2,1,0],yticklabels=labels,xlabel="珠部平均换热系数  [W/(m²·K)]")
    right.set_prop_cycle(color=colors)
    rings(right,body,old,"-",labels[0]);rings(right,pair,smooth,"--",labels[1])
    rings(right,wired,leads,"-.", "带线共同珠面",leads["regionID"]==4)
    right.set(xlabel="轴向位置 X  [mm]",ylabel="向空气放热通量  [kW/m²]")
    right.grid(alpha=.2);right.legend(frameon=False,loc="upper left",fontsize=9)
    fig.text(.055,.045,"平均系数按论文名义面积和 12 K 温差换算。曲线为各轴向区间的面积加权平均，+X 为迎流端。",fontsize=9,color="#4c5964")
    fig.savefig(HERE/"evidence/paired_heat.png",dpi=155)


if __name__=="__main__":main()
