"""从原生表面独立重建原面值，与页面、周向 CSV 和功率逐项比对，绘制同一截面的两种读数。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens
from cm011Tunnel import D, RHO, U, CP, DT, P0, meridian, write_json
from runEmbeddedSU2 import sha


def verify(flow,work,section):
    stats=json.loads((work/"audit/summary.json").read_bytes())
    page=work/"workbench/case-data.js"
    text=page.read_text(encoding="utf8").strip()
    assert text.startswith("window.ROCKET_CASE=") and text.endswith(";")
    model=json.loads(text[len("window.ROCKET_CASE="):-1])
    fields={f["key"]:f for f in model["fields"]}
    body=json.loads((flow/"body.json").read_bytes())
    vertices,faces=np.array(model["vertices"]),np.array(model["faces"])
    assert np.array_equal(vertices,body["V"]) and np.array_equal(faces+1,body["F"])
    wall,_=meridian(2)
    physical=wall*D
    distance,row=cKDTree(physical).query(np.column_stack((vertices[:,0],np.hypot(vertices[:,1],vertices[:,2]))))
    assert max(distance)<1e-13
    s=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(physical,axis=0),axis=1))]
    phi=np.angle(abs(vertices[:,1])+1j*vertices[:,2])
    rays=np.linspace(-np.pi/2,np.pi/2,13)
    surf,values=vtk_tokens(flow/"surface.vtk")
    grid=wall[::2]
    lookup=cKDTree(surf)
    surface_ids=[]
    for theta in rays:
        expected=np.column_stack((grid[:,0],grid[:,1]*np.cos(theta),grid[:,1]*np.sin(theta)))
        dist,ids=lookup.query(expected)
        assert max(dist)<2e-6
        surface_ids.append(ids)
    assert len(np.unique(surface_ids))==len(surf)
    points=np.column_stack((s[row],phi))
    max_differences={}
    for key,native,scale,offset in (("heatFlux","Heat_Flux",RHO*U*CP*DT/1000,0),
                                    ("pressureGauge","Pressure",RHO*U**2,0),
                                    ("pressure","Pressure",RHO*U**2/1000,P0/1000),
                                    ("wallFluidTemperature","Temperature",DT,300-DT)):
        samples=values[native][:,0][np.array(surface_ids)].T
        get=RegularGridInterpolator((s[::2],rays),samples)
        expected=get(points)[faces].mean(axis=1)*scale+offset
        f=fields[key]
        max_differences[key]=float(np.max(abs(np.array(f["values"])-expected)))
        np.testing.assert_allclose(f["values"],expected,rtol=2e-13,atol=1e-12)
        assert np.all(np.array(f["status"])==1) and "研究" in f["computedStateLabel"]
    table=np.genfromtxt(work/"audit/wall.csv",names=True,delimiter=",")
    np.testing.assert_allclose(fields["nearWallCellTemperature"]["values"],table["nearWallCellTemperatureK"],atol=1e-12,rtol=0)
    for k in ("wallSolved","edgeTemperature","recoveryTemperature"):
        assert set(fields[k]["status"])=={0}
    assert set(fields["wallInput"]["status"])=={2} and set(fields["wallInput"]["values"])=={300.1}
    area=.5*np.linalg.norm(np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]),axis=1)
    np.testing.assert_allclose(area,model["faceArea"],atol=1e-17,rtol=1e-13)
    power=area@(1000*np.array(fields["heatFlux"]["values"]))
    np.testing.assert_allclose(power,model["wallCoverage"]["wholeSurfaceHeatPowerW"],atol=1e-16,rtol=0)
    energy=model["nativeEnergyBalance"]
    assert sorted(energy["sourceFrameSHA256"].values())==sorted(stats["frameSHA256"].values())
    assert energy["sourceRestartSHA256"]==sha(flow/"restart.csv") and energy["belowHalfPercent"]
    assert model["resultCondition"]["alphaDeg"]==10 and model["resultCondition"]["betaDeg"]==0
    np.testing.assert_allclose(model["flowDirection"],[-np.cos(np.pi/18),0,np.sin(np.pi/18)],rtol=0,atol=1e-14)
    curve=json.loads((section/"curve.json").read_bytes())
    assert abs(curve["xM"]-physical[28,0])<1e-16
    for key in ("pressureGauge","heatFlux"):
        f=fields[key]
        with (section/(key+".csv")).open(encoding="utf-8-sig",newline="") as stream:
            records=list(csv.DictReader(stream))
        assert len(records)==curve["fields"][key]["segments"]
        for r in records:
            face=int(r["面元编号"])-1
            y=float(r[f["label"]+" / "+f["unit"]])
            assert y==f["values"][face] and r["算法"]==f["methodNames"][f["methodIndex"][face]]
            assert abs(float(r["截面 X / m"])-physical[28,0])<1e-16 and "+10°" in r["工况"]
    return dict(originalFaces=len(faces),areaM2=float(area.sum()),independentPowerW=float(power),
                maxDifferencesInViewerUnits=max_differences,fieldsChecked=list(fields),payloadSHA256=sha(page),
                sameRingX_M=curve["xM"],sectionSegments=curve["fields"]["heatFlux"]["segments"],
                curveUnits={key:curve["fields"][key]["unit"] for key in curve["fields"]},browserInteractionTested=False)


def figure(stats,curve,target):
    ring=stats["coneRing"]
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
    fig,(ax_p,ax_q)=plt.subplots(2,1,figsize=(8,5.7),sharex=True)
    for ax,key,node_key,scale,sample_scale,label in [
            (ax_p,"pressureGauge","gaugePressurePa",1e3,1e3,"相对环境壁压 / mPa"),
            (ax_q,"heatFlux","heatFluxWm2",1,1e3,"入壁热流 / (W/m²)")]:
        x=np.array(ring["azimuthDeg"])
        values=np.array(ring[node_key])*scale
        pairs=sorted({(round(a%360,7),v) for a,v in zip(np.r_[x,180-x],np.tile(values,2))})
        a,b=map(list,zip(*pairs))
        ax.plot(a+[360],b+[b[0]],"-o",color="#12647b",linewidth=1,markersize=3.8,label="求解网格壁节点（含镜像）")
        for i,p in enumerate(curve["fields"][key]["pieces"]):
            ax.plot([p["a"],p["b"]],[p["value"]*sample_scale]*2,color="#ae582d",linewidth=1.6,
                    label="工作台截面（面元常值）" if i==0 else None)
        for angle in (90,270): ax.axvline(angle,color="#7e858a",linestyle=":",linewidth=.8)
        ax.set_ylabel(label); ax.set_xlim(0,360)
        ax.grid(axis="y",color="#dadfe1",linewidth=.6)
        ax.spines[["top","right"]].set_visible(False)
    ax_p.set_title("CM-011 低速试件 · 攻角 +10°",loc="left",size=15,pad=14)
    ax_p.legend(loc="best",frameon=True,fontsize=9)
    ax_q.set_xticks([0,45,90,135,180,225,270,315,360],
                    ["0","45","90\n背流 +Z","135","180","225","270\n迎流 −Z","315","360"])
    ax_q.set_xlabel("周向角 / °（从 +Y 向 +Z）")
    fig.text(.98,.97,f"锥面截面 x = {ring['xM']*1000:.3f} mm · 来流 0.1 m/s",ha="right",va="top",size=9)
    fig.subplots_adjust(left=.15,right=.97,bottom=.16,top=.88,hspace=.27)
    fig.savefig(target,dpi=160)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("flow","workbench","section","new_output"): parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    report=verify(args.flow,args.workbench,args.section)
    args.new_output.mkdir(parents=True,exist_ok=False)
    write_json(args.new_output/"export.json",report)
    figure(json.loads((args.workbench/"audit/summary.json").read_bytes()),
           json.loads((args.section/"curve.json").read_bytes()),args.new_output/"cone_ring.png")
