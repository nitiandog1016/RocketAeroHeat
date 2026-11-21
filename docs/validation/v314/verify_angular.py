"""从两套实际终帧和独立热反力复算原面、七区、迎背风截面及工作台导出。"""
import argparse
import csv
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_angular import independent_overlap
sys.path.insert(0,str(HERE.parent/"v311"))
from verify_delivery import cut_values
from readCM011TunnelSector import vtk_tokens
from thermistorBead import read,sha
from thermistorReaction import terminal


def native(folder):
    """按真实原生三角面积恢复壁节点热流密度，并读取同帧压力、标量及壁摩擦。"""
    rec=terminal(folder);s=read(folder/"surface/summary.json");c=read(folder/"case.json")
    assert rec["outputSHA256"]["restart.csv"]==s["sourceFrameSHA256"]==s["energy"]["inputFrameSHA256"]
    assert sha(folder/"surface/wall.csv")==s["wallSHA256"] and sha(folder/"surface.vtk")==rec["outputSHA256"]["surface.vtk"]
    assert terminal(folder/"reaction0")["outputSHA256"]["restart.csv"]==s["sourceReactionSHA256"]
    assert s["relativeChange20"]<.002 and max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    assert read(folder/"ledger/energy.json")==s["energy"] and abs(s["energy"]["relativeDifference"])<.005
    with (folder/"reaction0/restart.csv").open(encoding="ascii") as f:
        h=next(csv.reader(f));rhs=np.loadtxt(f,delimiter=",",usecols=(h.index("Residual_Species_0"),))
    state=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1,usecols=(4,9))
    with np.load(folder/"nodes.npz") as z:
        xyz,wall,first=(z[k] for k in ("xyz","wallMap","firstMap"))
    tri=xyz[wall]*c["diameterM"];n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=np.linalg.norm(n,axis=1)/2;n/=2*area[:,None]
    lump=np.bincount(wall.ravel(),np.repeat(area/3,3),minlength=len(xyz))
    mass=np.zeros(len(xyz));wet=np.flatnonzero(lump);dt=c["wallK"]-c["airK"]
    mass[wet]=rhs[wet]*c["rho"]*c["speedMps"]*c["cp"]*dt*c["diameterM"]**2/lump[wet]
    np.testing.assert_allclose(-lump@mass,s["totalOutwardW"],rtol=0,atol=1e-13)
    nodes,flow=vtk_tokens(folder/"surface.vtk");dist,ix=cKDTree(nodes).query(xyz[wet])
    assert len(wet)==len(nodes)==len(np.unique(ix)) and max(dist)<2e-6
    friction=flow["Skin_Friction_Coefficient"][ix[np.searchsorted(wet,wall)]]
    height=np.einsum("fvi,fi->fv",xyz[first]-xyz[wall],n).mean(axis=1)*c["diameterM"]
    yplus=c["rho"]*c["speedMps"]/c["mu"]*height*np.sqrt(np.linalg.norm(friction,axis=2).mean(axis=1)/2)
    assert min(height)>0 and np.isfinite(yplus).all()
    return dict(wall=wall,first=first,xyz=xyz,tri=tri,area=area,qnode=mass,pnode=state[:,0]*c["rho"]*c["speedMps"]**2,
                temp=c["airK"]+dt*state[:,1],summary=s,yplus=[float(min(yplus)),float(max(yplus))])


def graph(profiles,out):
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    fig,axes=plt.subplots(2,2,figsize=(10,6),sharex=True,constrained_layout=True)
    fig.suptitle("带线热敏珠  ·  3 m/s 侧吹  ·  X = 0.172 mm")
    for row,(left,right,label,delta_label,factor) in enumerate((
            ("tipFineGauge","angular96Gauge","相对环境壁压（Pa）","表压变化（Pa）",1),
            ("tipFineHeat","angular96Heat","入壁热流（kW/m²）","热流变化（W/m²）",1000))):
        a,b=np.asarray(profiles[left]),np.asarray(profiles[right])
        orig,diff=axes[row]
        orig.plot(a[:,0],a[:,2],lw=2,label="48 档",color="#345f78")
        orig.plot(b[:,0],b[:,2],lw=1.2,ls="--",label="96 档",color="#ac5034")
        diff.plot(a[:,0],factor*(b[:,2]-a[:,2]),color="#695478",lw=1.2,marker=".",ms=2)
        orig.set_ylabel(label);diff.set_ylabel(delta_label)
    axes[0,0].set_title("两档在同一批原面上的结果");axes[0,1].set_title("96 档减 48 档")
    axes[0,0].legend(frameon=False)
    for ax in axes.flat:
        ax.grid(axis="y",alpha=.22);ax.axvline(180,ls=":",lw=.8,color=".6")
        ax.set_xlim(0,360);ax.set_xticks([0,60,120,180,240,300,360]);ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    for ax in axes[-1]:ax.set_xlabel("周向角（°）  0° 迎风，180° 背风")
    fig.savefig(out,dpi=160);plt.close(fig)


def verify(run,ui,out):
    prep=run.parent/"prepared01";inputs=read(prep/"inputs.json");coarse=Path(inputs["seed"])
    a=native(coarse);b=native(run);base=a["summary"];s=b["summary"]
    tri=a["tri"];area=a["area"];wall=a["wall"];n=len(area)
    with np.load(run/"nodes.npz") as z:
        pi,ci,edges=(z[k] for k in ("overlapParent","overlapChild","angularEdges"))
        np.testing.assert_array_equal(wall,z["originalWallMap"]);np.testing.assert_array_equal(a["first"],z["originalFirstMap"])
    np.testing.assert_array_equal(a["xyz"][a["first"]],b["xyz"][a["first"]])
    fractions=independent_overlap(tri,b["tri"],pi,ci);m=fractions*b["area"][ci,None]
    def original(v):
        loads=(v[b["wall"][ci]]*m).sum(axis=1)
        return np.bincount(pi,loads,minlength=n)/area
    def prolonged(v):
        return np.r_[v,v[edges].mean(axis=1)]
    q0=a["qnode"][wall].mean(axis=1);p0=a["pnode"][wall].mean(axis=1)
    q1=original(b["qnode"]);p1=original(b["pnode"])
    t0=a["temp"][a["first"]].mean(axis=1);t1=b["temp"][a["first"]].mean(axis=1)
    # 单独将原节点密度沿周向线性延拓；只量化壁面读数重构，不作为新网求解或误差上界。
    transferQ=original(prolonged(a["qnode"]))-q0
    transferP=original(prolonged(a["pnode"]))-p0
    def close(x,y,atol):
        np.testing.assert_allclose(x,y,rtol=0,atol=atol)
    for folder,p,q,t in ((coarse,p0,q0,t0),(run,p1,q1,t1)):
        rows=np.genfromtxt(folder/"surface/wall.csv",delimiter=",",names=True)
        np.testing.assert_array_equal(rows["faceID"],np.arange(n))
        close(rows["areaM2"],area,1e-19);close(rows["heatFluxWm2"],q,5e-8)
        close(rows["pressureGaugePa"],p,5e-11);close(rows["pressurePa"],p+101325,1e-9)
        close(rows["nearWallCellTemperatureK"],t,1e-12)
    receipt=read(ui/"debugger_process.json");u=read(ui/"ui.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not (receipt["survivingProcessIDs"] or receipt["stopReason"] or receipt["supervisorError"])
    assert u["faces"]==10272==n and u["webglError"]==0 and u["visibleSamples"]>0
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip()
    prefix="window.ROCKET_CASE=";assert raw.startswith(prefix) and raw.endswith(";")
    page=json.loads(raw[len(prefix):-1]);g=page["angularSpacingComparison"];fs={f["key"]:f for f in page["fields"]}
    assert g["baseFrame"]==u["baseFrame"]==base["sourceFrameSHA256"] and g["fineFrame"]==u["fineFrame"]==s["sourceFrameSHA256"]
    assert g["baseReaction"]==base["sourceReactionSHA256"] and g["fineReaction"]==s["sourceReactionSHA256"]
    assert g["mappingSHA256"]==sha(run/"nodes.npz") and g["mappingAuditSHA256"]==sha(prep/"mesh_audit.json")
    vals=dict(tipFineHeat=q0/1000,tipFineGauge=p0,tipFineFirst=t0,angular96Heat=q1/1000,angular96Gauge=p1,angular96First=t1,
              angularDeltaHeat=q1-q0,angularDeltaPressure=p1-p0,angularDeltaFirst=t1-t0)
    for k,v in vals.items():close(v,fs[k]["values"],5e-8)
    for k in ("angularDeltaHeat","angularDeltaPressure","angularDeltaFirst"):
        path=ui/(k+".csv");meta=read(path.with_suffix(".json"))
        assert meta["angularSpacingComparison"]==g and meta["csvSHA256"]==sha(path)
        with path.open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        assert len(rows)==n
        for j,r in enumerate(rows):
            assert int(r["face_id"])==j+1 and r["field_key"]==k and r["unit"]==fs[k]["unit"] and r["method"]
            assert r["availability_code"]=="2" and r["qualification_code"]=="0" and abs(float(r["value"])-vals[k][j])<5e-8
    profiles={}
    for k in ("tipFineHeat","angular96Heat","tipFineGauge","angular96Gauge"):
        profiles[k]=cut_values(ui/(k+"_section.csv"),tri,u["sectionXM"],vals[k],k,fs[k]["unit"])
    segments=[v[:2] for v in profiles["tipFineHeat"]]
    assert all([v[:2] for v in p]==segments for p in profiles.values())
    sides=[]
    for angle in (0,180):
        theta,id=min((v for v in segments if v[0]>=angle),key=lambda v:v[0]);j=id-1
        sides.append(dict(angleDeg=theta,faceID=id,gauge48Pa=p0[j],gauge96Pa=p1[j],intoWall48=q0[j],intoWall96=q1[j],
                          common48K=t0[j],common96K=t1[j],transferOnlyDeltaHeatWm2=transferQ[j],transferOnlyDeltaPressurePa=transferP[j]))
    body=read(run/"body.json");groups=np.asarray(body["regionID"]);regions=[]
    for j,r in enumerate(g["regions"],1):
        use=groups==j;old=float(-area[use]@q0[use]);new=float(-area[use]@q1[use])
        close([old,new],[r["outward48W"],r["outward96W"]],1e-13)
        regions.append(dict(name=r["name"],faces=int(use.sum()),areaM2=area[use].sum(),outward48W=old,outward96W=new,
                            deltaW=new-old,relativeChange=new/old-1))
    power0=-area@q0;power1=-area@q1
    close([power0,power1],[base["totalOutwardW"],s["totalOutwardW"]],1e-13)
    def measure(v):
        i=int(np.argmax(abs(v)));return dict(faceID=i+1,signedValue=v[i],areaRMS=np.sqrt(area@(v**2)/area.sum()))
    def ranking(v):
        ids=np.argsort(v)[:2]
        return [dict(faceID=int(i+1),region=body["regionNames"][groups[i]-1],q=v[i],centerM=tri[i].mean(axis=0).tolist()) for i in ids]
    # 浮点末位会改变几乎相等的两面的次序，正式面号按归档原值核验，独立重建按距离和数值核验。
    peaks=[]
    for v,k,src in ((q0,"base",base),(q1,"fine",s)):
        pair=ranking(v);ids=[g[k+"PeakFace"],g[k+"SecondFace"]]
        assert set(ids)=={z["faceID"] for z in pair} and ids[0]==src["peakFaceID"]
        peaks.append([next(z for z in pair if z["faceID"]==i) for i in ids])
    assert u["probeFace"]==g["finePeakFace"]
    hot_ids=sorted({r["faceID"] for pair in peaks for r in pair})
    hot=[dict(faceID=i,region=body["regionNames"][groups[i-1]-1],intoWall48=q0[i-1],intoWall96=q1[i-1],
              pressure48Pa=101325+p0[i-1],pressure96Pa=101325+p1[i-1],common48K=t0[i-1],common96K=t1[i-1],
              transferOnlyDeltaHeatWm2=transferQ[i-1]) for i in hot_ids]
    stats=dict(reportFaces=n,areaM2=area.sum(),newFaces=0,newAreaM2=0,jointExperimentalM2=0,regions=regions,
               outward48W=power0,outward96W=power1,powerRelativeChange=power1/power0-1,
               peaks48=peaks[0],peaks96=peaks[1],hotFacesCompared=hot,deltaHeatWm2=measure(q1-q0),deltaPressurePa=measure(p1-p0),
               deltaCommonK=measure(t1-t0),transferOnlyDeltaHeatWm2=measure(transferQ),transferOnlyDeltaPressurePa=measure(transferP),
               transferOnlyPowerChangeW=-area@transferQ,baseYPlus=a["yplus"],fineYPlus=b["yplus"],
               native96MinIntoWallWm2=float(b["qnode"][b["wall"]].mean(axis=1).min()),
               sectionXM=u["sectionXM"],segments=len(segments),sides=sides,baseFrame=g["baseFrame"],fineFrame=g["fineFrame"],
               meshAuditSHA256=g["mappingAuditSHA256"],residualWindow=s["residualWindow"],heatStability=s["relativeChange20"],
               energy=s["energy"],uiSeconds=receipt["elapsedSeconds"],pageSHA256=sha(ui/"workbench/case-data.js"))
    out.mkdir(parents=True,exist_ok=False)
    (out/"comparison.json").write_bytes((json.dumps(stats,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    graph(profiles,out/"section.png")
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def add(folder,names,tag):
            for name in names:z.write(folder/name,tag+"/"+name)
        add(prep,("inputs.json","mesh_audit.json","case.cfg","case.json"),"prepared")
        add(prep.parent,("prepared01_prepare.json","prepared01_audit.json"),"resources")
        for stage in sorted(run.parent.glob("steady*")):
            if (stage/"engine_process.json").is_file():
                names=["engine_request.json","engine_process.json","resources.json","history.csv","case.cfg"]
                if (stage/"field_audit.json").is_file():names.append("field_audit.json")
                add(stage,names,stage.name)
        for k in ("reaction0","plus20","reaction20"):add(run/k,("engine_process.json","resources.json","case.cfg"),k)
        add(run/"ledger",("energy.json",),"ledger")
        for k in ("control","observed"):add(run/"ledger"/k,("step_process.json",),"ledger/"+k)
        add(ui,("ui.json","debugger_process.json",*[k+"_section.csv" for k in profiles],
                "angularDeltaHeat.json","angularDeltaPressure.json","angularDeltaFirst.json"),"ui")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("run","ui","dest"):p.add_argument(name,type=Path)
    a=p.parse_args();verify(a.run.resolve(),a.ui.resolve(),a.dest.resolve())
