"""在同一批三角面上核对两域压力、壁热与首层温度，逐段验证周向曲线并绘图。"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZIP_DEFLATED

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v311")]
from verify_delivery import cut_values
from thermistorBead import read,sha
from thermistorReaction import terminal


def reference(name):
    path=ROOT/"cases/thermistor_b"/(name+".zip");key=read(path.with_suffix(".json"))
    assert sha(path)==key["archiveSHA256"]
    with ZipFile(path) as z:
        assert z.testzip() is None
        raw=z.read("wall.csv");s=json.loads(z.read("summary.json"))
        import hashlib
        assert hashlib.sha256(raw).hexdigest()==s["wallSHA256"]
        return np.genfromtxt(io.BytesIO(raw),delimiter=",",names=True),s,json.loads(z.read("body.json"))


def measure_delta(values,area):
    i=int(np.argmax(abs(values)))
    return dict(min=float(values.min()),max=float(values.max()),maxAbs=float(abs(values[i])),faceID=i+1,
                areaRMS=float(np.sqrt(area@(values**2)/area.sum())))


def plot(profiles,out):
    """共同色标容易掩盖小变化，所以同时给出原量及同面差值。"""
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    fig,axes=plt.subplots(2,2,figsize=(10,5.8),sharex=True,constrained_layout=True)
    fig.suptitle("带线热敏珠侧吹  ·  X = 0.172 mm  ·  来流 3 m/s",fontsize=14)
    for row,(base,expanded,label,diff_label,factor) in enumerate((
            ("pressureGauge","outerGauge30","相对环境壁压（Pa）","表压变化（mPa）",1000),
            ("heatFlux","outerHeat30","入壁热流（kW/m²）","入壁热流变化（W/m²）",1000))):
        a=np.asarray(profiles[base]);b=np.asarray(profiles[expanded]);left,right=axes[row]
        left.plot(a[:,0],a[:,2],color="#305a74",lw=2.1,label="20D")
        left.plot(b[:,0],b[:,2],color="#b65938",lw=1.1,ls="--",label="30D")
        right.plot(a[:,0],factor*(b[:,2]-a[:,2]),color="#604f79",lw=1.2,marker=".",ms=3)
        left.set_ylabel(label);right.set_ylabel(diff_label)
    axes[0,0].set_title("两域同面结果");axes[0,1].set_title("30D 减 20D");axes[0,0].legend(frameon=False)
    for ax in axes.flat:
        ax.grid(axis="y",alpha=.2);ax.axvline(180,color=".7",ls=":",lw=.8);ax.set_xlim(0,360)
        ax.set_xticks([0,60,120,180,240,300,360]);ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    for ax in axes[-1]:ax.set_xlabel("周向角（°）  0° 迎风，180° 背风")
    fig.savefig(out,dpi=165);plt.close(fig)


def verify(folder,ui,out):
    folder=folder.resolve();ui=ui.resolve();prep=folder.parent/"prepared01"
    a,s,body=reference("sidewind3");b,t,b2=reference("sidewindOuter30")
    assert body==b2 and s["sourceFrameSHA256"]==read(prep/"mesh_audit.json")["seedSHA256"]
    rec=terminal(folder);extra=terminal(folder/"plus20");react=terminal(folder/"reaction0")
    assert rec["outputSHA256"]["restart.csv"]==sha(folder/"restart.csv")==t["sourceFrameSHA256"]
    assert react["outputSHA256"]["restart.csv"]==t["sourceReactionSHA256"] and t==read(folder/"surface/summary.json")
    assert t["relativeChange20"]<.002 and abs(t["energy"]["relativeDifference"])<.005
    assert extra["inputSHA256"][str(folder/"plus20/initial_flow.csv")]==t["sourceFrameSHA256"]
    for k in ("areaM2","regionID","faceID"):np.testing.assert_array_equal(a[k],b[k])
    area=b["areaM2"];n=len(area);assert n==10272 and abs(float(area.sum())-3.19665061256752e-6)<1e-17
    src=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1)
    with np.load(folder/"nodes.npz") as z:wall=z["wallMap"];near=z["firstMap"]
    with (folder/"reaction0/restart.csv").open(encoding="ascii") as stream:
        header=next(csv.reader(stream));rhs=np.loadtxt(stream,delimiter=",",usecols=(0,header.index("Residual_Species_0")))
    np.testing.assert_array_equal(rhs[:,0],np.arange(len(src)))
    c=read(folder/"case.json");deltaT=c["wallK"]-c["airK"]
    assert c["outerRadiusOverD"]==30 and c["velocityMps"]==[0,-3,0] and deltaT==6
    weights=np.zeros(len(src))
    for v in wall.T:np.add.at(weights,v,area/3)
    q=(rhs[wall,1]/weights[wall]).mean(axis=1)*c["rho"]*c["speedMps"]*c["cp"]*deltaT*c["diameterM"]**2
    gauge=src[wall,4].mean(axis=1)*c["rho"]*c["speedMps"]**2
    first=c["airK"]+deltaT*src[near,9].mean(axis=1)
    def close(x,y,tol):np.testing.assert_allclose(x,y,rtol=0,atol=tol)
    close(q,b["heatFluxWm2"],1e-8);close(gauge,b["pressureGaugePa"],1e-12);close(first,b["nearWallCellTemperatureK"],1e-12)
    close(b["pressurePa"],c["airPa"]+gauge,1e-10)
    dq=b["heatFluxWm2"]-a["heatFluxWm2"];dp=b["pressurePa"]-a["pressurePa"];dt=first-a["nearWallCellTemperatureK"]
    raw=(ui/"workbench/case-data.js").read_text(encoding="utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    page=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in page["fields"]}
    receipt=read(ui/"debugger_process.json");u=read(ui/"ui.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not (receipt["stopReason"] or receipt["supervisorError"] or receipt["survivingProcessIDs"])
    assert u["webglError"]==0 and u["visibleSamples"]>0 and u["faces"]==n
    assert page["resultCondition"]["physics"]["velocityMps"]==[0,-3,0] and page["resultCondition"]["wallUniformValueK"]==299.15
    expected={"outerDeltaHeat":dq,"outerDeltaPressure":dp,"outerDeltaFirst":dt,"heatFlux":a["heatFluxWm2"]/1000,
              "outerHeat30":b["heatFluxWm2"]/1000,"pressureGauge":a["pressureGaugePa"],"outerGauge30":b["pressureGaugePa"],
              "nearWallCellTemperature":a["nearWallCellTemperatureK"],"outerFirst30":first}
    for key,v in expected.items():close(fields[key]["values"],v,1e-10)
    pair=page["outerDomainComparison"];side=read(ui/"deltaHeat.json")
    assert side["outerDomainComparison"]==pair and side["csvSHA256"]==sha(ui/"deltaHeat.csv")
    assert pair["baseFrame"]==s["sourceFrameSHA256"] and pair["expandedFrame"]==t["sourceFrameSHA256"]
    assert pair["expandedReaction"]==t["sourceReactionSHA256"] and pair["initialInteriorFrame"]==s["sourceFrameSHA256"]
    assert sum([max(s["residualWindow"].values())<-8,max(t["residualWindow"].values())<-8,max(t["plus20Window"].values())<-8])==3
    with (ui/"deltaHeat.csv").open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==n
    for i,row in enumerate(rows):
        assert int(row["face_id"])==i+1 and row["unit"]=="W/m²" and row["qualification_code"]=="0" and row["availability_code"]=="2"
        assert abs(float(row["value"])-dq[i])<1e-9
    regions=[]
    for j,r in enumerate(pair["regions"],1):
        use=b["regionID"]==j;old=float(-area[use]@a["heatFluxWm2"][use]);new=float(-area[use]@b["heatFluxWm2"][use])
        close([old,new],[r["outward20W"],r["outward30W"]],1e-13)
        regions.append(dict(name=r["name"],faces=int(use.sum()),areaM2=float(area[use].sum()),outward20W=old,outward30W=new,
                            deltaW=new-old,relativeChange=new/old-1))
    pw0=float(-area@a["heatFluxWm2"]);pw1=float(-area@b["heatFluxWm2"])
    close([pw0,pw1],[s["totalOutwardW"],t["totalOutwardW"]],1e-13)
    tri=np.asarray(body["V"])[np.asarray(body["F"])-1];section={}
    for key,unit in (("pressureGauge","Pa"),("heatFlux","kW/m²"),("outerGauge30","Pa"),("outerHeat30","kW/m²")):
        section[key]=cut_values(ui/(key+"_section.csv"),tri,.000172,expected[key],key,unit)
    indices=[x[:2] for x in section["heatFlux"]]
    assert all([x[:2] for x in v]==indices for v in section.values())
    sides=[]
    for angle in (0,180):
        chosen=min((v for v in indices if v[0]>=angle),key=lambda v:v[0])
        theta,id=chosen;i=id-1
        sides.append(dict(angleDeg=theta,faceID=id,gauge20Pa=float(a["pressureGaugePa"][i]),gauge30Pa=float(b["pressureGaugePa"][i]),
                          intoWall20Wm2=float(a["heatFluxWm2"][i]),intoWall30Wm2=float(b["heatFluxWm2"][i])))
    peak0=int(np.argmin(a["heatFluxWm2"]));peak1=int(np.argmin(b["heatFluxWm2"]))
    assert peak0+1==pair["peak20Face"]==s["peakFaceID"] and peak1+1==pair["peak30Face"]==t["peakFaceID"]
    assert u["probeFace"]==int(np.argmax(abs(dq)))+1
    stats=dict(faces=n,areaM2=float(area.sum()),outward20W=pw0,outward30W=pw1,deltaW=pw1-pw0,relativeChange=pw1/pw0-1,
               regions=regions,heatFluxDifferenceWm2=measure_delta(dq,area),pressureDifferencePa=measure_delta(dp,area),
               firstLayerDifferenceK=measure_delta(dt,area),peak20Face=peak0+1,peak30Face=peak1+1,peak20IntoWallWm2=float(a["heatFluxWm2"][peak0]),
               peak30IntoWallWm2=float(b["heatFluxWm2"][peak1]),heatPeakSeparationM=float(np.linalg.norm(tri[peak1].mean(axis=0)-tri[peak0].mean(axis=0))),
               sectionXM=.000172,segments=len(indices),sectionSides=sides,source20=pair["baseFrame"],source30=pair["expandedFrame"],
               validatedJointAreaM2=0,ledger=t["energy"],relativeHeatChange20=t["relativeChange20"],residualWindow=t["residualWindow"],
               uiSeconds=receipt["elapsedSeconds"],solverSeconds=rec["elapsedSeconds"],surfaceSHA256=t["wallSHA256"],
               pageSHA256=sha(ui/"workbench/case-data.js"),meshAuditSHA256=sha(prep/"mesh_audit.json"))
    out.mkdir(parents=True,exist_ok=False)
    (out/"comparison.json").write_bytes((json.dumps(stats,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    plot(section,out/"section_domains.png")
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def add(dir,names,tag):
            for name in names:z.write(dir/name,tag+"/"+name)
        add(prep,("mesh_audit.json","case.cfg","case.json"),"prepared")
        for stage in sorted(folder.parent.glob("steady*")):
            if (stage/"field_audit.json").is_file():
                add(stage,("engine_request.json","engine_process.json","field_audit.json","resources.json","history.csv","case.cfg"),stage.name)
                script=ROOT/"scripts/thermistorSidewindExtent.py";original=stage/"executed_extent.py"
                file=original if original.is_file() else script
                assert sha(file)==read(stage/"engine_process.json")["inputSHA256"][str(script)]
                z.write(file,stage.name+"/thermistorSidewindExtent.py")
        for name in ("plus20","reaction0","reaction20"):add(folder/name,("engine_process.json","resources.json","case.cfg"),name)
        add(folder/"ledger",("energy.json",),"ledger")
        for name in ("control","observed"):add(folder/"ledger"/name,("step_process.json",),"ledger/"+name)
        add(ui,("ui.json","debugger_process.json","deltaHeat.json",*[k+"_section.csv" for k in section]),"ui")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("run","ui","out"):p.add_argument(key,type=Path)
    a=p.parse_args();verify(a.run,a.ui,a.out)
