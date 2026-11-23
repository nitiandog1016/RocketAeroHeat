"""逐原面复核两个测温入口、实验比较分母、导出及原生回执。"""
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
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v311")]
from thermistorBead import read,sha
from verify_delivery import cut_values


def load_wall(stem):
    with ZipFile(CASES/(stem+".zip")) as z:
        return np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=","),json.loads(z.read("summary.json"))


def verify(ui,run,out):
    old,ref=load_wall("longWire12");direct,ds=load_wall("resistance12")
    k=read(CASES/"resistance12.json");input_rec=read(run/"request.json");f=k["factor"];air=input_rec["ambientK"]
    assert sha(CASES/"resistance12.zip")==k["archiveSHA256"] and sha(run/"surface/wall.csv")==k["wallSHA256"]
    assert ds["sourceFrameSHA256"]==k["frameSHA256"] and ref["sourceFrameSHA256"]==k["referenceFrame"]
    u=read(ui/"ui.json");receipt=read(ui/"debugger_process.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not(receipt["stopReason"] or receipt["supervisorError"] or receipt["survivingProcessIDs"])
    assert u["originalFrame"]==ref["sourceFrameSHA256"]
    body=read(run/"prepared/body.json");tri=np.array(body["V"])[np.array(body["F"])-1]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    group=np.array(body["regionID"]);bead=np.isin(group,(3,4,5));outcomes={}
    def close(a,b,tol):np.testing.assert_allclose(a,b,rtol=0,atol=tol)
    for action,source,multiplier in (("solve",old,f),("archive",direct,1)):
        d=ui/action;proof=u[action]
        raw=(d/"workbench/case-data.js").read_text("utf8").strip();head="window.ROCKET_CASE="
        assert raw.startswith(head) and raw.endswith(";");page=json.loads(raw[len(head):-1])
        w=page["thermistorWire"];m=w["thermometry"];fields={field["key"]:field for field in page["fields"]}
        assert proof["sourceFrame"]==w["sourceFrameSHA256"]==(k["referenceFrame"] if action=="solve" else k["frameSHA256"])
        assert proof["faces"]==10272 and proof["webglError"]==0 and proof["pixels"]>0 and w["jointExperimentVerifiedM2"]==0
        assert page["resultCondition"]["wallInputType"]=="perFace" and page["resultCondition"]["wallUniformValueK"] is None
        assert page["resultCondition"]["physics"]["resistanceTemperatureK"]==m["internalK"]==305.15
        assert "wallSolved" not in fields or all(v is None for v in fields["wallSolved"]["values"])
        q=source["heatFluxWm2"]*multiplier
        first=air+multiplier*(source["nearWallCellTemperatureK"]-air);tw=air+multiplier*(source["wallInputK"]-air)
        diff=q-old["heatFluxWm2"]
        for field,val in (("pressure",source["pressurePa"]/1000),("pressureGauge",source["pressureGaugePa"]),("heatFlux",q/1000),
                          ("nearWallCellTemperature",first),("wallInput",tw),("thermometryHeatFluxChange",diff/1000),
                          ("wireHeatFluxChange",(q-f*source["uniformHeatFluxWm2"])/1000)):
            close(fields[field]["values"],val,1e-9)
        watts=-area*q;bead_w=float(watts[bead].sum());datum=w["experiment"]["hWm2K"];meas_w=datum*12*3.094e-6
        close(m["apparentHWm2K"],bead_w/(12*3.094e-6),1e-9)
        close(m["surfaceHWm2K"],bead_w/((m["surfaceK"]-air)*3.094e-6),1e-9)
        close(m["remainingGapW"],bead_w-meas_w,1e-12)
        close(m["fractionOfPowerGap"],(ref["beadAndJointsW"]-bead_w)/(ref["beadAndJointsW"]-meas_w),1e-10)
        bands=w["regions"]
        for i,band in enumerate(bands,1):close(band["outwardW"],watts[group==i].sum(),1e-13)
        close(w["totalOutwardW"],watts.sum(),1e-13);close(w["energy"]["wallReleaseW"],watts.sum(),1e-11)
        with (d/"wall.csv").open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        assert len(rows)==len(area)
        for i,r in enumerate(rows):
            assert int(r["face_id"])==i+1 and r["thermometry_reference_frame_sha256"]==k["referenceFrame"]
            assert r["source_frame_sha256"]==proof["sourceFrame"] and r["wall_boundary_sha256"]==k["boundarySHA256"]
            assert float(r["resistance_temperature_K"])==305.15 and abs(float(r["inferred_bead_surface_K"])-m["surfaceK"])<1e-10
            assert abs(float(r["source_field_temperature_ratio"])-multiplier)<1e-13 and abs(float(r["uniform_reference_temperature_ratio"])-f)<1e-13
            assert abs(float(r["wall_input_K"])-tw[i])<1e-10 and abs(float(r["pressure"])*1000-source["pressurePa"][i])<1e-8
            assert abs(float(r["heat_flux"])*1000-q[i])<1e-7 and abs(float(r["thermometry_change_kW_m2"])-diff[i]/1000)<1e-10
        change=read(d/"change.json");assert change["csvSHA256"]==sha(d/"change.csv") and change["thermistorWire"]==w
        with (d/"change.csv").open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        assert len(rows)==len(q)
        for i,r in enumerate(rows):
            assert int(r["face_id"])==i+1 and r["unit"]=="kW/m²" and r["availability_code"]=="2" and r["qualification_code"]=="0"
            assert abs(float(r["value"])-diff[i]/1000)<1e-9
        segments=cut_values(d/"section.csv",tri,proof["sectionXM"],diff/1000,"thermometryHeatFluxChange","kW/m²")
        assert proof["peakFace"]==w["peakFaceID"]==int(np.argmin(q))+1
        outcomes[action]=dict(wholeW=float(watts.sum()),beadW=bead_w,apparentH=m["apparentHWm2K"],surfaceH=m["surfaceHWm2K"],
                              regionW=[b["outwardW"] for b in bands],peakFace=proof["peakFace"],peakWm2=min(q),
                              firstK=[min(first),max(first)],wallK=[min(tw),max(tw)],segments=len(segments),source=proof["sourceFrame"],
                              pageSHA256=sha(d/"workbench/case-data.js"),sectionSHA256=sha(d/"section.csv"))
    m=u["solve"]["model"];reading=m["figure12"];assert reading["range"][0]<=m["ratio"]<=reading["range"][1]
    exp=w["experiment"]
    result=dict(method=m,areas=dict(computedM2=sum(area),faces=len(q),newGeometryFaces=0,newGeometryM2=0,jointExperimentalM2=0),
                originalBeadW=ref["beadAndJointsW"],experimentEquivalentW=exp["hWm2K"]*12*3.094e-6,experiment=exp,paths=outcomes,
                actualEnergy=ds["energy"],stability20=ds["relativeChange20"],window=ds["window"],plus20Window=ds["plus20Window"],
                unmatchedToMain=k["unmatchedToMain"],uiSeconds=receipt["elapsedSeconds"])
    out.mkdir(parents=True,exist_ok=False)
    (out/"comparison.json").write_bytes((json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    plot(result,out/"thermometry.png")
    stage=run/"steady01"
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        for n in ("request.json","completed.json","prepared/profile.json","surface/summary.json"):z.write(run/n,"case/"+n)
        for d in (stage,*(stage/n for n in ("reaction0","plus20","reaction20","step_control","step_observed"))):
            for n in ("case.cfg","history.csv","native.json","debugger_process.json"):z.write(d/n,str((d/n).relative_to(run)))
        for n in ("ui.json","debugger_process.json","solve/section.csv","archive/section.csv"):z.write(ui/n,"ui/"+n)


def plot(result,path):
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    m=result["method"];fig,(left,right)=plt.subplots(1,2,figsize=(9.2,3.5),constrained_layout=True)
    re=np.geomspace(1,1e3,160);r=m["ratio"]*(2/3+.3*np.sqrt(re))/(2/3+.3*np.sqrt(m["reynoldsL"]))
    left.semilogx(re,100*r,color="#425e75",label="附录式（2′）· 恒物性")
    g=m["figure12"];left.errorbar(g["reynoldsL"],100*g["value"],yerr=[[100*(g["value"]-g["range"][0])],[100*(g["range"][1]-g["value"])]],
                                  fmt="o",color="#ac622d",capsize=4,label="图 12 读数及墨迹范围")
    left.set(xlabel="长度雷诺数",ylabel="内部与表面温差 / 表面超温（%）",title="B 珠测温模型");left.legend(frameon=False,fontsize=8)
    exp=result["experiment"];interval=np.array(exp["graphReadRangeWm2K"])*3.094e-6*12*1000
    right.axvspan(*interval,color="#deddd8",label="图 7 印刷读图范围")
    for y,value,color in [(2,1000*result["originalBeadW"],"#425e75"),(1,1000*result["paths"]["solve"]["beadW"],"#287670"),(0,1000*result["experimentEquivalentW"],"#9c5e34")]:
        right.scatter(value,y,color=color,s=36);right.annotate(f"{value:.3f}",(value,y),xytext=(5,4),textcoords="offset points",fontsize=9)
    right.set(yticks=[0,1,2],yticklabels=["图 7 折算","按电阻测温","原给定表温"],xlabel="珠部向空气放热（mW）",title="3 m/s · 名义电阻超温 12 K",xlim=(7.65,10.8),ylim=(-.5,2.5))
    right.legend(frameon=False,fontsize=8,loc="center right",bbox_to_anchor=(.95,.25))
    for ax in (left,right):ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False);ax.grid(axis="x",alpha=.13)
    fig.savefig(path,dpi=170);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("ui","run","out"):p.add_argument(name,type=Path)
    a=p.parse_args();verify(a.ui.resolve(),a.run.resolve(),a.out.resolve())
