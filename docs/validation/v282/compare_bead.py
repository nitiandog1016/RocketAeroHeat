"""比较三套网格的逐面计算与图 7 的单个整珠换热读数。实验数据只参与后处理。"""
from itertools import product
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "outputs/v282_validation"
SOURCE = "c4c81c99f21ab0b42037d09aaadead1aa28b80e387ab0e6dc238c4e23f16daed"


def digitize():
    crop = RUN / "source01/figure7.png"
    assert hashlib.sha256(crop.read_bytes()).hexdigest() == SOURCE
    # 像素以裁图左上角为原点，刻度在测点相邻处定位，保留印刷线宽和手工定位范围。
    x_ref, y_ref, point = (92.5, 588.), (363.5, 264.5), (122.5, 270.5)
    factor = 4184 / 3600  # 热化学千卡；每小时换成每秒
    def axis(pixel, ref, span):
        return span * (pixel-ref[0]) / (ref[1]-ref[0])
    def limits(pixel, ref, span):
        vals = [axis(pixel+d, (ref[0]+e, ref[1]+f), span)
                for d,e,f in product((-2.,2.), (-.75,.75), (-.75,.75))]
        return [min(vals), max(vals)]
    h = axis(point[1], y_ref, 200)
    return dict(doi="10.1299/kikai1938.25.810", page=815, figure=7, specimen="B",
                speedMps=3., cropSHA256=SOURCE, cropPixels=[627,457], origin="top-left",
                xPixelsAt0and200K=x_ref, yPixelsAt0and200kcal=y_ref, pointXY=point,
                pointHalfWidthPx=2, tickHalfWidthPx=.75, deltaTK=axis(point[0],x_ref,200),
                deltaTReadRangeK=limits(point[0],x_ref,200), hKcalPerM2HourK=h,
                joulesPerKcal=4184, hWm2K=h*factor,
                graphReadRangeWm2K=[v*factor for v in limits(point[1],y_ref,200)],
                reportedPointUncertaintyWm2K=None)


def compare(out):
    out.mkdir(parents=True, exist_ok=True)
    obs = digitize()
    grids = {}
    rows = {}
    for grid in ("coarse","fine","outer"):
        folder = RUN/f"result_{grid}01"
        grids[grid] = json.loads((folder/"summary.json").read_bytes())
        assert hashlib.sha256((folder/"wall.csv").read_bytes()).hexdigest()==grids[grid]["wallTableSHA256"]
        rows[grid] = np.genfromtxt(folder/"wall.csv", delimiter=",", names=True)
    body = json.loads((RUN/"strong_outer01/body.json").read_bytes())
    faces = np.array(body["F"])-1
    rings = np.asarray(body["nodeRow"])[faces].min(axis=1)
    assert np.array_equal(np.unique(rings), np.arange(64))
    pairs = {}
    for a,b in (("coarse","fine"),("fine","outer")):
        old,new = rows[a],rows[b]
        assert np.array_equal(old["faceID"],new["faceID"]) and np.array_equal(old["areaM2"],new["areaM2"])
        pairs[a+"To"+b] = dict(hRelativeToFirst=grids[b]["nominalExperimentAreaCoefficientWm2K"] /
                                grids[a]["nominalExperimentAreaCoefficientWm2K"]-1,
                                hDifferenceRelativeToSecond=abs(grids[a]["powerW"]/grids[b]["powerW"]-1),
                                maxLocalQChangeRelativeToFirst=float(max(abs(new["heatFluxWm2"]/old["heatFluxWm2"]-1))),
                                maxLocalPressureChangePa=float(max(abs(new["pressurePa"]-old["pressurePa"]))))
    predicted = grids["outer"]["nominalExperimentAreaCoefficientWm2K"]
    lo,hi = obs["graphReadRangeWm2K"]
    result = dict(experiment=obs, predictedNominalAreaHWm2K=predicted, relativeDifference=predicted/obs["hWm2K"]-1,
                  relativeDifferenceAtReadLimits=[predicted/hi-1,predicted/lo-1], gridComparisons=pairs)
    (out/"comparison.json").write_text(json.dumps(result,indent=2)+"\n", encoding="utf8")
    # 左图展示沿全珠母线的环带面积平均；右图的误差棒只覆盖图读数范围。
    fig,(ax,avg) = plt.subplots(1,2,figsize=(9,3.65),gridspec_kw={"width_ratios":[1.6,1]})
    for grid,color,style in (("coarse","#949da3",":"),("fine","#c56429","--"),("outer","#175577","-")):
        q = rows[grid]["heatFluxWm2"]
        area = rows[grid]["areaM2"]
        ring_flux = np.bincount(rings, weights=area*q)/np.bincount(rings,weights=area)
        ax.plot((63.5-np.arange(64))*180/64,ring_flux/1000,style,color=color,label=grid)
    ax.set(xlabel="Angle from upstream pole (deg)",ylabel="Heat flux into bead (kW/m²)",xlim=(0,180))
    ax.legend(frameon=False, ncol=3, fontsize=8)
    avg.errorbar(0, obs["hWm2K"], yerr=[[obs["hWm2K"]-lo],[hi-obs["hWm2K"]]],fmt="o",capsize=5,color="#913e24")
    avg.plot(1,predicted,"o",color="#175577")
    avg.set(xticks=[0,1], xticklabels=["Fig. 7 reading","Simulation"], xlim=(-.6,1.6),ylim=(180,280),
            ylabel="Whole-bead h (W/m²/K)")
    for a in (ax,avg):
        a.grid(axis="y",alpha=.23)
        a.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out/"distribution.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    compare(Path(sys.argv[1]))
