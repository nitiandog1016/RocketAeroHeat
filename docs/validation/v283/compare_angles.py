"""在同一原面上比较两种周向离散，保存整珠热量、局部差值和可复核的统计。"""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts"))
from thermistorBead import D, read, sha


def compare_field(coarse, refined, area):
    difference = refined-coarse
    face = int(np.argmax(abs(difference)))
    return dict(areaWeightedRMS=float(np.sqrt(np.average(difference**2, weights=area))),
                minChange=float(difference.min()), maxChange=float(difference.max()),
                maxAbsChangeFaceID=face, coarseAtFace=float(coarse[face]), refinedAtFace=float(refined[face]))


def compare(out):
    base=ROOT/"outputs/v282_validation"
    run=ROOT/"outputs/v283_validation"
    old=read(base/"result_outer01/summary.json")
    new=read(run/"result01/summary.json")
    assert new["bodySHA256"]==old["bodySHA256"] and new["configSHA256"]==old["configSHA256"]
    assert new["angularCells"]==4 and new["fullSurfaceComputed"]
    paths=[base/"result_outer01/wall.csv",run/"result01/wall.csv"]
    for path,summary in zip(paths,[old,new]):
        assert sha(path)==summary["wallTableSHA256"]
    a,b=[np.genfromtxt(p,names=True,delimiter=",") for p in paths]
    area=a["areaM2"]
    np.testing.assert_array_equal(a["faceID"],np.arange(5216))
    np.testing.assert_array_equal(a["faceID"],b["faceID"])
    np.testing.assert_array_equal(area,b["areaM2"])
    np.testing.assert_allclose([area@t["heatFluxWm2"] for t in (a,b)],
                               [old["powerW"],new["powerW"]],rtol=2e-14,atol=0)
    assert sha(run/"strong01/body.json")==new["bodySHA256"]
    body=read(run/"strong01/body.json")
    faces=np.array(body["F"])-1
    centers=np.array(body["V"])[faces].mean(axis=1)
    rings=np.asarray(body["nodeRow"])[faces].min(axis=1)
    energy=read(run/"energy01/energy.json")
    assert energy["inputFrameSHA256"]==new["sourceFrameSHA256"] and energy["sampledThetaExact"]
    assert abs(energy["wallReleaseW"]+new["powerW"])<1e-11 and abs(energy["relativeDifference"])<.005
    assert abs(energy["farOutW"]+new["probeAllResidualW"])<1e-11
    observation=read(ROOT/"docs/validation/v282/evidence/comparison.json")["experiment"]
    h0,h1=[s["nominalExperimentAreaCoefficientWm2K"] for s in (old,new)]
    dq=b["heatFluxWm2"]-a["heatFluxWm2"]
    peak=[int(np.argmin(t["heatFluxWm2"])) for t in (a,b)]
    note=dict(baselineSHA256=old["sourceFrameSHA256"],refinedSHA256=new["sourceFrameSHA256"],
              faceCount=len(area),areaM2=area.sum(),newlyCoveredFaces=0,verifiedJointAreaM2=0,
              meanH15=h0,meanH7p5=h1,experimentH=observation["hWm2K"],hChangeFraction=h1/h0-1,
              experimentalDeviation15=h0/observation["hWm2K"]-1,
              experimentalDeviation7p5=h1/observation["hWm2K"]-1,
              peakFluxFaces=peak,peakFluxCentroidsM=centers[peak].tolist(),
              peakReleaseFluxWm2=[float(-t["heatFluxWm2"].min()) for t in (a,b)],
              heatFluxWm2=compare_field(a["heatFluxWm2"],b["heatFluxWm2"],area),
              maxLocalRelativeFluxChange=float(max(abs(dq/a["heatFluxWm2"]))),
              gaugePressurePa=compare_field(a["pressureGaugePa"],b["pressureGaugePa"],area),
              firstNodeTemperatureK=compare_field(a["nearWallCellTemperatureK"],b["nearWallCellTemperatureK"],area))
    for tag,folder in (("baseline",base/"strong_outer01"),("refined",run/"strong01")):
        geometry=read(folder/"geometry.json")
        with np.load(folder/"nodes.npz") as grid:
            ids=grid["ids"]; coords=grid["xyz"]
        m=65
        first=D*np.linalg.norm(coords[ids[:,m:2*m]]-coords[ids[:,:m]],axis=2)
        note[tag+"Grid"]=dict(firstLayerMicron=[first.min()*1e6,first.max()*1e6],
                               nodes=len(coords),cells=geometry["hexahedra"]+geometry["prisms"],
                               fullChordWallAreaM2=geometry["fullChordWallAreaM2"])
    # 输出文件只含本次比较。实验点及测量定义跟随 v282 原读数。
    out.mkdir(parents=True,exist_ok=False)
    (out/"comparison.json").write_text(json.dumps(note,indent=2)+"\n",encoding="utf8")
    np.savetxt(out/"angular_difference.csv",np.column_stack((a["faceID"],area,dq,b["pressureGaugePa"]-a["pressureGaugePa"])),
               delimiter=",",fmt=["%d"]+["%.17g"]*3,comments="",
               header="nativeFaceID,areaM2,heatFluxChangeWm2,gaugePressureChangePa")
    angle=(63.5-np.arange(64))*180/64
    def ring_mean(values):
        return np.bincount(rings,weights=area*values)/np.bincount(rings,weights=area)
    fig,(flow,delta)=plt.subplots(2,1,figsize=(7.6,5.8),sharex=True,gridspec_kw={"height_ratios":[1.35,1]})
    for values,label,color,style in [(a["heatFluxWm2"],"15°","#bd672c","--"),
                                      (b["heatFluxWm2"],"7.5°","#1e5975","-")]:
        flow.plot(angle,ring_mean(values)/1000,style,color=color,label=label)
    flow.legend(frameon=False);flow.set_ylabel("Heat flux into bead (kW/m²)")
    low=np.array([dq[rings==j].min() for j in range(64)])
    high=np.array([dq[rings==j].max() for j in range(64)])
    delta.fill_between(angle,low,high,color="#b6ced9",label="Range over faces in ring")
    delta.plot(angle,ring_mean(dq),color="#1e5975",label="Area-weighted ring mean")
    delta.axhline(0,color="0.45",linewidth=.7);delta.legend(frameon=False,fontsize=8)
    delta.set(xlabel="Angle from upstream pole (deg)",ylabel="7.5° − 15° (W/m²)",xlim=(0,180))
    for ax in (flow,delta):
        ax.grid(axis="y",alpha=.22);ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout();fig.savefig(out/"angular_difference.png",dpi=160);plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    compare(parser.parse_args().output)
