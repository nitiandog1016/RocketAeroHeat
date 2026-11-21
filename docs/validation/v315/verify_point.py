"""从两份原生边界反力重建测点场，并核对工作台原面、共同温度、风洞口径及图读数据。"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CASES=ROOT/"cases/thermistor_b"
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v314")]
from audit_angular import independent_overlap
sys.path.insert(0,str(HERE.parent/"v311"))
from verify_delivery import cut_values
from thermistorBead import read,sha
from thermistorReaction import terminal


def check(run6,run12,ui,dest):
    proof=read(ui/"ui.json");rec=read(ui/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["survivingProcessIDs"] or rec["stopReason"] or rec["supervisorError"])
    c=read(run12/"case.json");D=c["diameterM"];body=read(run12/"body.json")
    tri=np.asarray(body["V"])[np.asarray(body["F"])-1];groups=np.asarray(body["regionID"])
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    with np.load(run12/"nodes.npz") as z:
        xyz,wall,first,pi,ci=(z[k] for k in ("xyz","wallMap","originalFirstMap","overlapParent","overlapChild"))
    small=xyz[wall]*D;a=np.linalg.norm(np.cross(small[:,1]-small[:,0],small[:,2]-small[:,0]),axis=1)/2
    w=independent_overlap(tri,small,pi,ci)*a[ci,None]
    lump=np.bincount(wall.ravel(),np.repeat(a/3,3),minlength=len(xyz));wet=np.flatnonzero(lump)
    def project(v):return np.bincount(pi,(w*v[wall[ci]]).sum(axis=1),minlength=len(area))/area
    def reconstruct(folder,unit):
        ledger=terminal(folder);s=read(folder/"surface/summary.json");cc=read(folder/"case.json")
        assert sha(folder/"nodes.npz")==s["nodeMapSHA256"]==sha(run12/"nodes.npz")
        assert ledger["outputSHA256"]["restart.csv"]==s["sourceFrameSHA256"]==s["energy"]["inputFrameSHA256"]
        assert terminal(folder/"reaction0")["outputSHA256"]["restart.csv"]==s["sourceReactionSHA256"]
        assert s["energy"]["temperatureScaleK"]==unit and abs(s["energy"]["relativeDifference"])<.005 and s["relativeChange20"]<.002
        assert max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
        with (folder/"reaction0/restart.csv").open(encoding="ascii") as f:
            h=next(csv.reader(f));rhs=np.loadtxt(f,delimiter=",",usecols=h.index("Residual_Species_0"))
        st=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1,usecols=(4,9))
        density=np.zeros(len(xyz));density[wet]=rhs[wet]*c["rho"]*3*c["cp"]*D**2*unit/lump[wet]
        assert np.all(st[wet,1]*unit==cc["wallK"]-cc["airK"])
        return project(density),project(c["rho"]*9*st[:,0])+c["airPa"],c["airK"]+unit*st[first,1].mean(axis=1),s
    q6,p6,t6,s6=reconstruct(run6,6);q12,p12,t12,s12=reconstruct(run12,24)
    pred=2*q6;tPred=c["airK"]+2*(t6-c["airK"]);frame0=read(run12/"field_audit.json")["sourceFrame"]
    assert frame0==sha(run6/"initial_flow.csv") and read(run6/"field_audit.json")["steps"]==read(run12/"field_audit.json")["steps"]==40
    def close(x,y,tol):np.testing.assert_allclose(x,y,rtol=0,atol=tol)
    with ZipFile(CASES/"sidewindAngular96.zip") as z:previous=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True)
    with ZipFile(CASES/"leadsThird3.zip") as z:
        axial=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True);axS=json.loads(z.read("summary.json"))
    e=read(HERE.parent/"v282/evidence/comparison.json")["experiment"];B=3.094e-6;bead=np.isin(groups,[3,4,5])
    axW=float(-area[bead]@axial["heatFluxWm2"][bead]);close(axW,axS["beadAndJointsW"],1e-13)
    key=read(CASES/"windPoint12.json");assert key["archiveSHA256"]==sha(CASES/"windPoint12.zip")
    with ZipFile(CASES/"windPoint12.zip") as z:
        assert z.testzip() is None
        for prefix,folder in (("six_",run6),("",run12)):
            for n in ("wall.csv","summary.json"):
                assert z.read(prefix+n)==(folder/"surface"/n).read_bytes()
    cases={}
    for action,q,p,T,s in (("solve",pred,p6,tPred,s6),("archive",q12,p12,t12,s12)):
        folder=ui/action;raw=(folder/"workbench/case-data.js").read_text("utf8").strip();head="window.ROCKET_CASE="
        assert raw.startswith(head) and raw.endswith(";");data=json.loads(raw[len(head):-1]);m=data["measuredPoint"]
        assert proof[action]["frame"]==s["sourceFrameSHA256"] and proof[action]["pixels"]>0 and proof[action]["webglError"]==0
        fs={f["key"]:f for f in data["fields"]};assert data["resultCondition"]["wallUniformValueK"]==305.15 and m["action"]==action
        for k,values,tol in (("pressure",p/1000,1e-11),("pressureGauge",p-c["airPa"],1e-9),("heatFlux",q/1000,1e-9),
                             ("nearWallCellTemperature",T,1e-10),("wallInput",np.full(len(area),305.15),0)):
            close(fs[k]["values"],values,tol)
        if "wallSolved" in fs:assert all(v is None for v in fs["wallSolved"]["values"])
        meta=read(folder/"heat.json");assert meta["measuredPoint"]==m and meta["csvSHA256"]==sha(folder/"heat.csv")
        with (folder/"heat.csv").open(encoding="utf-8-sig",newline="") as stream:
            rows=list(csv.DictReader(stream));assert len(rows)==len(area)
            for j,row in enumerate(rows):
                assert int(row["face_id"])==j+1 and row["qualification_code"]=="0" and row["availability_code"]=="2"
                assert row["unit"]=="kW/m²" and abs(float(row["value"])-q[j]/1000)<1e-9
        with (folder/"wall.csv").open(encoding="utf-8-sig",newline="") as stream:
            rows=list(csv.DictReader(stream));assert len(rows)==len(area)
            for j,row in enumerate(rows):
                assert int(row["face_id"])==j+1 and row["reference_frame_sha256"]==s["sourceFrameSHA256"] and row["experiment_doi"]==e["doi"]
                assert abs(float(row["pressure"])*1000-p[j])<1e-8 and abs(float(row["heat_flux"])*1000-q[j])<1e-7
                assert row["pressure_unit"]=="kPa" and row["heat_flux_unit"]=="kW/m²" and float(row["temperature_difference_ratio"])==(2 if action=="solve" else 1)
        profiles=cut_values(folder/"section.csv",tri,proof[action]["sectionXM"],q/1000,"heatFlux","kW/m²")
        regions=[float(-area[groups==i]@q[groups==i]) for i in range(1,8)]
        close(regions,[item["outwardW"] for item in m["regions"]],1e-13)
        power=sum(regions);qb=sum(regions[2:5]);h=qb/(B*12);he=e["hWm2K"]
        close([power,qb,h,h/he-1,axW/(B*12),he*B*12],[m["wholeW"],m["beadW"],m["beadHWm2K"],
              m["relativeDifference"],m["axial"]["hWm2K"],m["equivalentExperimentPowerW"]],1e-10)
        assert m["experiment"]["reportedPointUncertaintyWm2K"] is None and m["experiment"]["measuredAxisAngleDeg"] is None
        assert m["experiment"]["graphReadRangeWm2K"]==e["graphReadRangeWm2K"] and m["jointVerifiedM2"]==0
        corr=m["correlation"];fluid=corr["air"];nu=fluid["dynamicViscosityPas"]/fluid["densityKgm3"];el=.53
        shape=16*el/(3*(np.sqrt(1-el**2)+np.arcsin(el)/el)**2*np.log((1+el)/(1-el)))
        formula=(shape+.29*np.sqrt(3*.343e-3/nu))*fluid["conductivityWmK"]/.343e-3
        close(formula,corr["meanHWm2K"],1e-10)
        i=int(np.argmin(q));assert i+1==m["peakFaceID"]==proof[action]["peakFace"]
        cases[action]=dict(frame=s["sourceFrameSHA256"],qPeakFace=i+1,qPeakWm2=q[i],wallPowerW=power,beadW=qb,leadsW=power-qb,
                           hWm2K=h,relative=h/he-1,regionsW=regions,segments=len(profiles),pressurePa=[min(p),max(p)],firstK=[min(T),max(T)],
                           sectionSHA256=sha(folder/"section.csv"),pageSHA256=sha(folder/"workbench/case-data.js"))
    matched=dict(maxPressurePa=float(max(abs(p12-p6))),maxFirstK=float(max(abs(t12-tPred))),maxHeatWm2=float(max(abs(q12-pred))),
                 maxRelativeHeat=float(max(abs((q12-pred)/pred))),relativePower=abs(s12["totalOutwardW"]/(2*s6["totalOutwardW"])-1))
    for k,val in matched.items():close(val,m["sourceCheck"][k],1e-8)
    assert matched["maxPressurePa"]<.001 and matched["maxFirstK"]<.001 and matched["maxRelativeHeat"]<.0005 and matched["relativePower"]<.0005
    unmatched=max(abs(p12-previous["pressurePa"]));drift=max(abs(p6-previous["pressurePa"]))
    close([unmatched,drift],[m["unmatched"]["maxPressurePa"],m["sixKPressureDriftPa"]],1e-9)
    result=dict(cases=cases,areaM2=float(area.sum()),faces=len(area),jointValidatedM2=0,experiment=m["experiment"],
                axial=m["axial"],correlation=m["correlation"],firstUnmatched=m["unmatched"],sixKPressureDriftPa=drift,matched=matched,
                graphRangeRelative=m["relativeAtGraphLimits"],equivalentExperimentPowerW=e["hWm2K"]*B*12,
                originalSixKFrame=frame0,energy6=s6["energy"],energy12=s12["energy"],window6=s6["residualWindow"],window12=s12["residualWindow"],
                stability6=s6["relativeChange20"],stability12=s12["relativeChange20"],uiSeconds=rec["elapsedSeconds"])
    dest.mkdir(parents=True,exist_ok=False)
    (dest/"comparison.json").write_bytes((json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    fig,ax=plt.subplots(figsize=(8,3.35),constrained_layout=True)
    ax.axvspan(*e["graphReadRangeWm2K"],color="#dedfd9",label="图 7 印刷读图范围");ax.axvline(e["hWm2K"],color="#555b5a",lw=1)
    values=[cases["solve"]["hWm2K"],m["axial"]["hWm2K"],formula,e["hWm2K"]]
    for y,(v,col) in enumerate(zip(values,["#316582","#597c96","#96663e","#454748"])):
        ax.scatter(v,y,s=48,color=col);ax.annotate(f"{v:.2f}",(v,y),xytext=(6,4),textcoords="offset points")
    ax.set(yticks=range(4),yticklabels=["名义侧吹 · 96 档","名义轴向 · 第三档","原文拟合式（28）","图 7 实验点"],ylim=(3.5,-.6),
           xlim=(205,290),xlabel="珠部平均换热系数  h（W/(m²·K)）",title="B 型热敏珠  ·  3 m/s  ·  温差约 12 K")
    ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False);ax.grid(axis="x",alpha=.15);ax.legend(loc="lower right",frameon=False)
    fig.savefig(dest/"mean_heat.png",dpi=170);plt.close(fig)
    with ZipFile(dest/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def save(folder,names,tag):
            for n in names:z.write(folder/n,tag+"/"+n)
        for path in (run6,run12):
            tag=path.name;save(path,("engine_process.json","engine_request.json","resources.json","history.csv","field_audit.json","case.cfg"),tag)
            for n in ("reaction0","reaction20","plus20"):save(path/n,("engine_process.json","resources.json","case.cfg"),tag+"/"+n)
            save(path/"ledger",("energy.json",),tag+"/ledger")
            for n in ("control","observed"):save(path/"ledger"/n,("step_process.json",),tag+"/ledger/"+n)
        save(ui,("ui.json","debugger_process.json"),"ui")
        for action in cases:save(ui/action,("section.csv","heat.json"),"ui/"+action)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("six","twelve","ui","dest"):p.add_argument(name,type=Path)
    a=p.parse_args();check(a.six.resolve(),a.twelve.resolve(),a.ui.resolve(),a.dest.resolve())
