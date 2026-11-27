"""以同一正投影绘制低速风洞两张原面云图，色标采用有符号入壁热流。"""
from zipfile import ZipFile

from refined import DOC, read, sha
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from simpleLauncherTunnel import CASE, geometry


def plot():
    _,tri,_,_=geometry()
    folder=CASE/"wake"
    receipt=read(folder/"receipt.json")
    fields=[]
    for bundle,key,name in (("wind_tunnel.zip","sha256","wall.csv"),
                            ("wake_heat_grid.zip","wakeHeatGridSHA256","fine_wall.csv")):
        assert sha(folder/bundle)==receipt[key]
        with ZipFile(folder/bundle) as z:
            wall=np.loadtxt(z.open(name),delimiter=",",skiprows=1)
        np.testing.assert_array_equal(wall[:,0],np.arange(1,len(tri)+1))
        fields.append(1000*wall[:,3])
    coarse,fine=fields
    # 稍从平底一侧观察，按三角形深度由远到近绘制。色彩只编码数值。
    eye=np.array([-.42,-1.,.48]);eye/=np.linalg.norm(eye)
    right=np.cross([0.,0.,1.],eye);right/=np.linalg.norm(right)
    up=np.cross(eye,right)
    uv=tri@np.column_stack((right,up))
    order=np.argsort(tri.mean(axis=1)@eye)
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False})
    fig=plt.figure(figsize=(9,4.7),facecolor="white")
    fig.text(.055,.93,"简化箭体全表面换热",fontsize=19,color="#213444")
    fig.text(.055,.865,"来流 0.05 m/s     距弱扰动初场 0.600 s     给定壁温 300.001 K     空气 300 K",fontsize=10,color="#364958")
    delta=fine-coarse
    for row,(values,limits,cmap,title) in enumerate((
        (fine,(fine.min(),fine.max()),"viridis","近壁细网 · 入壁热流"),
        (delta,(-abs(delta).max(),abs(delta).max()),"RdBu_r","同面差值 · 细网减粗网"))):
        bottom=.50-row*.395
        ax=fig.add_axes([.055,bottom,.66,.31])
        poly=PolyCollection(uv[order],cmap=cmap,edgecolors="face",linewidths=.33)
        poly.set_array(values[order]);poly.set_clim(*limits)
        ax.add_collection(poly);ax.autoscale();ax.set_aspect("equal");ax.set_axis_off()
        fig.text(.735,bottom+.22,title,fontsize=11,color="#213444")
        cb=fig.colorbar(poly,cax=fig.add_axes([.735,bottom+.143,.205,.03]),orientation="horizontal")
        cb.set_ticks(np.linspace(*limits,3));cb.ax.set_xticklabels([f"{x:.2f}" for x in np.linspace(*limits,3)])
        cb.ax.tick_params(labelsize=9,length=2);cb.outline.set_linewidth(.5)
        fig.text(.735,bottom+.047,"mW/m²",fontsize=9,color="#364958")
    fig.text(.055,.035,"负热流表示箭体向空气放热。每个三角面按原面值着色。",fontsize=9,color="#364958")
    fig.savefig(DOC/"surface.png",dpi=155);plt.close(fig)


if __name__=="__main__":
    plot()
