"""从热敏珠两套终帧与壁反力独立恢复场量，核验端面近壁对照、面元导出及引线截面。"""
import argparse
import csv
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(HERE.parent/"v311"))
from verify_delivery import cut_values  # v311 的三角形／平面独立求交
from readCM011TunnelSector import vtk_tokens
from thermistorBead import read,sha
from thermistorReaction import terminal


def reconstruct(folder):
    """在每个原壁顶点恢复强标量壁反力，使用各自实网和终帧读取温度、压力、y+。"""
    record=terminal(folder)
    summary=read(folder/"surface/summary.json");c=read(folder/"case.json")
    body=read(folder/"body.json");faces=np.asarray(body["F"])-1
    tri=np.asarray(body["V"])[faces];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=np.linalg.norm(normal,axis=1)/2;normal/=2*area[:,None]
    assert record["outputSHA256"]["restart.csv"]==summary["sourceFrameSHA256"]
    last=read(folder/"field_audit.json")
    assert last["allResidualsBelowMinus8"] and last["lastFrame"]==summary["sourceFrameSHA256"]
    assert sha(folder/"surface.vtk")==record["outputSHA256"]["surface.vtk"]
    assert sha(folder/"surface/wall.csv")==summary["wallSHA256"]
    assert terminal(folder/"reaction0")["outputSHA256"]["restart.csv"]==summary["sourceReactionSHA256"]
    state=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1)
    with (folder/"reaction0/restart.csv").open(encoding="ascii") as stream:
        header=next(csv.reader(stream))
        rows=np.loadtxt(stream,delimiter=",",usecols=(0,header.index("Residual_Species_0")))
    np.testing.assert_array_equal(rows[:,0],np.arange(len(state)))
    with np.load(folder/"nodes.npz") as z:
        wall,first,xyz=(z[k] for k in ("wallMap","firstMap","xyz"))
        common=z["commonMap"] if "commonMap" in z else first
    np.testing.assert_allclose(xyz[wall]*c["diameterM"],tri,rtol=0,atol=1e-16)
    weights=np.bincount(wall.ravel(),weights=np.repeat(area/3,3),minlength=len(state))
    scale=c["rho"]*c["speedMps"]*c["cp"]*(c["wallK"]-c["airK"])*c["diameterM"]**2
    q=scale*(rows[wall,1]/weights[wall]).mean(axis=1)
    gauge=c["rho"]*c["speedMps"]**2*state[wall,4].mean(axis=1)
    firstK=c["airK"]+(c["wallK"]-c["airK"])*state[first,9].mean(axis=1)
    matched=c["airK"]+(c["wallK"]-c["airK"])*state[common,9].mean(axis=1)
    height=c["diameterM"]*np.einsum("fvi,fi->fv",xyz[first]-xyz[wall],normal).mean(axis=1)
    nodes,flow=vtk_tokens(folder/"surface.vtk");wet=np.unique(wall)
    dist,ix=cKDTree(nodes).query(xyz[wet])
    assert len(nodes)==len(np.unique(ix))==len(wet) and dist.max()<2e-6
    friction=flow["Skin_Friction_Coefficient"][ix[np.searchsorted(wet,wall)]]
    yplus=height*c["rho"]*c["speedMps"]/c["mu"]*np.sqrt(np.linalg.norm(friction,axis=2).mean(axis=1)/2)
    table=np.genfromtxt(folder/"surface/wall.csv",delimiter=",",names=True)
    np.testing.assert_array_equal(table["faceID"],np.arange(len(q)))
    for actual,col,tol in ((q,"heatFluxWm2",1e-8),(gauge,"pressureGaugePa",1e-12),
                           (firstK,"nearWallCellTemperatureK",1e-12),(area,"areaM2",1e-19)):
        np.testing.assert_allclose(actual,table[col],rtol=0,atol=tol)
    return dict(body=body,tri=tri,area=area,q=q,p=gauge,first=firstK,matched=matched,commonXYZ=xyz[common],
                height=height,yplus=yplus,summary=summary,pabs=table["pressurePa"])


def draw_cap(tri, base, fine, out):
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    fig,axes=plt.subplots(1,3,figsize=(10.8,4.4),constrained_layout=True)
    fig.suptitle("前引线端面放热  ·  3 m/s 侧吹  ·  来流沿 −Y",fontsize=14)
    lo=min(base.min(),fine.min())/1000;hi=max(base.max(),fine.max())/1000
    delta=fine-base;bound=max(abs(delta));points=tri[:,:,1:]*1e6
    for i,(ax,value,title) in enumerate(zip(axes,(base/1000,fine/1000,delta),("原网","局部细网","细网减原网"))):
        limits=(-bound,bound) if i==2 else (lo,hi)
        polys=PolyCollection(points,array=value,edgecolors="none",cmap="RdBu_r" if i==2 else "viridis",clim=limits)
        ax.add_collection(polys);ax.autoscale();ax.set(aspect="equal",title=title,xlabel="Y（µm）")
        marked=np.argsort(value)[-2:] if i<2 else [np.argmax(abs(value))]
        ax.add_collection(PolyCollection(points[marked],facecolors="none",edgecolors="#222222",linewidths=.75))
        if i==0:ax.set_ylabel("Z（µm）")
        if i==1:fig.colorbar(polys,ax=axes[:2],orientation="horizontal",shrink=.75,label="向空气放热（kW/m²）")
        if i==2:fig.colorbar(polys,ax=ax,orientation="horizontal",label="放热变化（W/m²）")
    fig.savefig(out,dpi=160);plt.close(fig)


def verify(run,ui,dest):
    prepared=run.parent/"prepared01"
    source=Path(read(prepared/"inputs.json")["seed"])
    a=reconstruct(source);b=reconstruct(run);g=b["summary"];base=a["summary"]
    assert a["body"]==b["body"] and np.array_equal(a["commonXYZ"],b["commonXYZ"])
    tri=b["tri"];area=b["area"];n=len(area);assert n==10272
    dq=b["q"]-a["q"];dp=b["pabs"]-a["pabs"];dt=b["matched"]-a["matched"]
    model=read(ui/"ui.json");rec=read(ui/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert model["faces"]==n and model["webglError"]==0 and model["visibleSamples"]>0
    raw=(ui/"workbench/case-data.js").read_text(encoding="utf8").strip()
    prefix="window.ROCKET_CASE=";assert raw.startswith(prefix) and raw.endswith(";")
    page=json.loads(raw[len(prefix):-1]);pair=page["tipSpacingComparison"]
    assert pair["baseFrame"]==model["baseFrame"]==base["sourceFrameSHA256"]
    assert pair["fineFrame"]==model["fineFrame"]==g["sourceFrameSHA256"]
    assert pair["fineReaction"]==g["sourceReactionSHA256"] and pair["baseReaction"]==base["sourceReactionSHA256"]
    fields={f["key"]:f for f in page["fields"]}
    expected={"tipFineHeat":b["q"]/1000,"tipDeltaHeat":dq,"tipFineGauge":b["p"],"tipDeltaPressure":dp,
              "tipFineFirst":b["first"],"tipMatched":b["matched"],"tipDeltaMatched":dt,"tipHeight":b["height"]*1e6}
    for k,val in expected.items():
        np.testing.assert_allclose(fields[k]["values"],val,rtol=0,atol=2e-10)
    for k in ("tipDeltaHeat","tipDeltaMatched"):
        path=ui/(k+".csv");side=read(path.with_suffix(".json"))
        assert side["csvSHA256"]==sha(path) and side["tipSpacingComparison"]==pair
        with path.open(encoding="utf-8-sig",newline="") as fi:
            rows=list(csv.DictReader(fi))
        assert len(rows)==n
        for j,r in enumerate(rows):
            assert int(r["face_id"])==j+1 and r["availability_code"]=="2" and r["qualification_code"]=="0"
            assert r["method"] and r["field_key"]==k and r["unit"]==fields[k]["unit"]
            assert abs(float(r["value"])-expected[k][j])<2e-10
    section=cut_values(ui/"tip_section.csv",tri,model["sectionXM"],expected["tipFineHeat"],"tipFineHeat","kW/m²")
    groups=np.asarray(b["body"]["regionID"]);cap=groups==7
    assert all(groups[id-1]==6 for _,id,_ in section)
    peak0=int(np.argmin(a["q"]));peak1=int(np.argmin(b["q"]))
    assert peak0+1==pair["basePeakFace"]==10267 and peak1+1==pair["finePeakFace"]==model["probeFace"]
    ranking=[]
    for field,tag in ((a,"base"),(b,"fine")):
        top=np.argsort(field["q"])[:2]
        assert top[1]+1==pair[tag+"SecondFace"]
        gap=float(np.diff(field["q"][top])[0])
        assert abs(gap-pair[tag+"PeakGapWm2"])<3e-11
        ranking.append(dict(faceIDs=(top+1).tolist(),gapWm2=gap,centersM=tri[top].mean(axis=1).tolist()))
    mirrored=tri[peak0]*[1,1,-1]
    mirrorError=float(cKDTree(tri[peak1]).query(mirrored)[0].max())
    assert mirrorError<1e-15 and peak0!=peak1 and min(a["q"][peak0],b["q"][peak1])<0
    faces=np.asarray(b["body"]["F"]);neighbors=np.isin(faces,faces[peak0]).sum(axis=1)>=2
    assert neighbors.sum()==4
    local=[dict(faceID=int(i+1),region=b["body"]["regionNames"][groups[i]-1],baseQ=a["q"][i],fineQ=b["q"][i],
                deltaQ=dq[i],deltaCommonK=dt[i]) for i in np.flatnonzero(neighbors)]
    regs=[]
    for j,r in enumerate(pair["regions"],1):
        mask=groups==j
        old=float(-area[mask]@a["q"][mask]);new=float(-area[mask]@b["q"][mask])
        np.testing.assert_allclose([old,new],[r["baseOutW"],r["fineOutW"]],rtol=0,atol=1e-13)
        faceIDs=np.flatnonzero(mask)
        i0=faceIDs[np.argmin(a["q"][mask])];i1=faceIDs[np.argmin(b["q"][mask])]
        regs.append(dict(name=r["name"],faces=int(sum(mask)),baseW=old,fineW=new,deltaW=new-old,relative=new/old-1,
                         basePeakFace=int(i0+1),finePeakFace=int(i1+1),basePeakQ=a["q"][i0],finePeakQ=b["q"][i1]))
    total0=float(-area@a["q"]);total1=float(-area@b["q"])
    np.testing.assert_allclose([total0,total1],[base["totalOutwardW"],g["totalOutwardW"]],rtol=0,atol=1e-13)
    def delta(values):
        ix=int(np.argmax(abs(values)))
        return dict(faceID=ix+1,value=values[ix],maxAbs=abs(values[ix]),areaRMS=np.sqrt(area@(values**2)/area.sum()))
    def limits(v):return [float(v.min()),float(v.max())]
    result=dict(faces=n,areaM2=area.sum(),newFaces=0,newAreaM2=0,jointExperimentVerifiedM2=0,
                baseFrame=pair["baseFrame"],fineFrame=pair["fineFrame"],baseOutW=total0,fineOutW=total1,
                totalRelative=total1/total0-1,regions=regs,basePeakFace=peak0+1,finePeakFace=peak1+1,
                basePeakQ=a["q"][peak0],finePeakQ=b["q"][peak1],peakSeparationM=np.linalg.norm(tri[peak0].mean(axis=0)-tri[peak1].mean(axis=0)),
                topTwo=ranking,mirrorVertexMaxDistanceM=mirrorError,neighborsOfBasePeak=local,peakPatchAreaM2=area[neighbors].sum(),
                basePatchW=-area[neighbors]@a["q"][neighbors],finePatchW=-area[neighbors]@b["q"][neighbors],
                deltaHeatWm2=delta(dq),deltaPressurePa=delta(dp),deltaCommonK=delta(dt),
                basePeakFirstK=a["first"][peak0],fineFirstAtBasePeakK=b["first"][peak0],fineCommonAtBasePeakK=b["matched"][peak0],
                baseCapYPlus=limits(a["yplus"][cap]),fineCapYPlus=limits(b["yplus"][cap]),fineFullYPlus=limits(b["yplus"]),
                capFineYPlusAtBasePeak=b["yplus"][peak0],sectionXM=model["sectionXM"],segments=len(section),
                ledger=g["energy"],relativeHeatChange20=g["relativeChange20"],residualWindow=g["residualWindow"],
                steps=g["cumulativeSteps"],uiSeconds=rec["elapsedSeconds"],pageSHA256=sha(ui/"workbench/case-data.js"))
    dest.mkdir(parents=True,exist_ok=False)
    (dest/"comparison.json").write_bytes((json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    draw_cap(tri[cap],-a["q"][cap],-b["q"][cap],dest/"cap_map.png")
    with ZipFile(dest/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def add(folder,names,tag):
            for name in names:z.write(folder/name,tag+"/"+name)
        add(prepared,("inputs.json","mesh_audit.json","case.json","case.cfg"),"mesh")
        add(prepared.parent,("prepared01_prepare.json","prepared01_audit.json"),"preparation")
        for stage in sorted(run.parent.glob("steady*")):
            if (stage/"field_audit.json").is_file():
                add(stage,("engine_request.json","engine_process.json","field_audit.json","history.csv","resources.json","case.cfg"),stage.name)
        for stage in ("reaction0","reaction20","plus20"):add(run/stage,("engine_process.json","resources.json","case.cfg"),stage)
        for stage in ("observed","control"):add(run/"ledger"/stage,("step_process.json",),"ledger/"+stage)
        add(run/"ledger",("energy.json",),"ledger")
        add(ui,("ui.json","debugger_process.json","tip_section.csv","tipDeltaHeat.json","tipDeltaMatched.json"),"ui")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("run","ui","out"):p.add_argument(k,type=Path)
    a=p.parse_args();verify(a.run.resolve(),a.ui.resolve(),a.out.resolve())
