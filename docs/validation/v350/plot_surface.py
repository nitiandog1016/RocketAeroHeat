"""从发布 STL 与逐面 CSV 绘制两张共用几何视角的定量云图。"""
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.font_manager import FontProperties, fontManager
import numpy as np


def main():
    root=Path(__file__).resolve().parents[3];case=root/"cases/simple_launcher"
    shape=json.loads((case/"geometry.json").read_bytes())
    dtype=np.dtype([("normal","<f4",3),("points","<f4",(3,3)),("attribute","<u2")])
    tr=np.fromfile(case/shape["file"],dtype=dtype,offset=84)["points"].astype(float)
    with ZipFile(case/"wind_tunnel.zip") as z:
        summary=json.loads(z.read("summary.json"))
        with z.open("wall.csv") as f:rows=np.loadtxt(f,delimiter=",",skiprows=1)
    assert len(tr)==len(rows) and np.array_equal(rows[:,0],np.arange(1,len(tr)+1))
    pressure=1000*(rows[:,2]-summary["condition"]["ambientPressurePa"]);heat=1000*rows[:,3]
    view=np.array([-.49,-.79,.38]);view/=np.linalg.norm(view)
    right=np.cross([0.,0.,1.],view);right/=np.linalg.norm(right);up=np.cross(view,right)
    polygons=np.stack((tr@right,tr@up),axis=2);order=np.argsort((tr@view).mean(axis=1))
    font_path=Path("C:/Windows/Fonts/msyh.ttc")
    if font_path.exists():
        fontManager.addfont(str(font_path));plt.rcParams["font.family"]=FontProperties(fname=str(font_path)).get_name()
    plt.rcParams["axes.unicode_minus"]=False
    fig=plt.figure(figsize=(12,5.2),facecolor="white")
    fig.text(.055,.89,"低速风洞 · 全表面计算",size=21,weight="bold",color="#203644")
    fig.text(.055,.81,"轴向来流 0.05 m/s    |    空气 300 K    |    壁面 300.001 K    |    9,216 个原面",size=11,color="#48606d")
    panels=[("相对环境壁压","mPa",pressure,"coolwarm",TwoSlopeNorm(0.,pressure.min(),pressure.max())),
            ("入壁热流密度","mW/m²",heat,"viridis",Normalize(heat.min(),heat.max()))]
    for i,(title,unit,values,cmap,norm) in enumerate(panels):
        left=.055+.49*i
        fig.text(left,.694,f"{title}  ({unit})",size=13,weight="bold",color="#213b47")
        ax=fig.add_axes([left,.31,.405,.33]);ax.set_aspect("equal");ax.axis("off")
        colors=plt.get_cmap(cmap)(norm(values[order]))
        ax.add_collection(PolyCollection(polygons[order],facecolors=colors,edgecolors="face",linewidths=.38))
        ax.set_xlim(polygons[:,:,0].min()-.028,polygons[:,:,0].max()+.028)
        ax.set_ylim(polygons[:,:,1].min()-.018,polygons[:,:,1].max()+.018)
        cax=fig.add_axes([left,.232,.403,.023])
        ticks=[values.min(),0.,values.max()] if i==0 else [values.min(),-5.,-2.5,values.max()]
        bar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal",ticks=ticks)
        bar.set_ticklabels([f"{x:.2f}" for x in ticks]);bar.outline.set_visible(False)
        bar.ax.tick_params(labelsize=9,length=2,color="#526773")
    fig.text(.055,.075,"环境参考压力 101325 Pa。入壁为正；本工况下壁面向空气放热。",size=10,color="#425864")
    fig.savefig(case/"wind-tunnel.png",dpi=150,facecolor="white");plt.close(fig)


if __name__=="__main__":main()
