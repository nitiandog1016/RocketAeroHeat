"""从两份发布的原面表绘制整周加权的轴向热流与四带放热变化。"""
import json
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from check_grid import HERE, ROOT, read
from compare_meshes import axial_position, table


def main():
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":10,"axes.unicode_minus":False,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,(flux,bands)=plt.subplots(1,2,figsize=(10,4),constrained_layout=True,gridspec_kw={"width_ratios":[1.5,1]})
    fig.suptitle("光滑 B 型热敏珠    3 m/s 空气  ·  给定壁温 305.15 K",fontsize=14)
    for name,label,color in (("paired3","加密前","#596b7a"),("midbody3","中段加密","#19766b")):
        with zipfile.ZipFile(ROOT/f"cases/thermistor_b/{name}.zip") as z:
            x=axial_position(json.loads(z.read("body.json")))
            w=table(z.read("wall.csv"));strip=table(z.read("flow.csv"))["meridian_strip_0"]
        centers=[];heat=[]
        for j in np.unique(strip):
            mask=strip==j
            centers.append(np.average(x[mask],weights=w["areaM2"][mask]))
            heat.append(np.average(-w["heatFluxWm2"][mask]/1000,weights=w["areaM2"][mask]))
        flux.plot(centers,heat,label=label,color=color,lw=2)
    flux.set(xlabel="轴向位置 X/a（+1 为迎流端）",ylabel="周向平均向外热流  [kW/m²]",xlim=(-1,1))
    flux.legend(frameon=False);flux.grid(axis="y",alpha=.18)
    rows=read(HERE/"evidence/comparison.json")["nativeBands"]
    delta=[1000*(r["outwardW"]["refined"]-r["outwardW"]["baseMixed"]) for r in rows]
    bands.barh(range(4),delta,color="#19766b",height=.5)
    bands.set(xlabel="放热变化，加密后减加密前  [mW]",yticks=range(4),
              yticklabels=["−1 至 −0.5","−0.5 至 0","0 至 0.5","0.5 至 1"],ylabel="轴向位置 X/a",xlim=(-1.1,0))
    for i,value in enumerate(delta):
        bands.text(value-.02,i,f"{value:+.3f}",ha="right",va="center",fontsize=9)
    bands.invert_yaxis()
    fig.savefig(HERE/"evidence/axial_heat.png",dpi=160)
    plt.close(fig)


if __name__=="__main__":
    main()
