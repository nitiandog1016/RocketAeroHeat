"""用独立几何交叠、终态和同位置均温对照，核验长引线全场及工作台的定量出口。"""
import argparse
import csv
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v314"),str(HERE.parent/"v311")]
from audit_angular import independent_overlap
from verify_delivery import cut_values
from thermistorBead import read,sha
from thermistorReaction import terminal


def verify(work,matched,ui,out):
    s=read(work/"surface/summary.json");job=read(work/"completed.json");w=work/job["stage"];prep=work/"prepared"
    base=Path(read(prep/"profile.json")["base"]);model=read(prep/"case.json");body=read(prep/"body.json")
    ui_receipt=read(ui/"debugger_process.json");proof=read(ui/"ui.json")
    assert ui_receipt["status"]=="terminal" and ui_receipt["debuggerExitCode"]==0 and all(ui_receipt["inputsStable"].values())
    assert not(ui_receipt["stopReason"] or ui_receipt["supervisorError"] or ui_receipt["survivingProcessIDs"])
    assert proof["archiveSameWall"] and proof["changedBoundaryInvalidated"] and proof["bands"]==12 and proof["webglError"]==0 and proof["pixels"]>0
    saved=read(work.parent/"resume.json")
    assert saved["first"]["stage"]=="uniform" and saved["resume"]["stage"]=="uniform_control" and saved["resume"]["status"]=="paused"
    assert all(sha(prep/name)==digest for name,digest in saved["inputs"].items())
    assert proof["frame"]==s["sourceFrameSHA256"]==sha(w/"restart.csv")
    assert sha(work/"surface/wall.csv")==s["wallSHA256"] and sha(prep/"profile.dat")==s["boundarySHA256"]
    matched_stage=matched/w.name;u=read(matched_stage/"surface/summary.json")
    in_wire=sha(prep/"initial_flow.csv");in_uniform=sha(base/"restart.csv");steps=[]
    for dest in sorted(work.glob("steady[0-9][0-9]")):
        official=matched/dest.name;nr=read(dest/"native.json");ur=read(official/"field_audit.json")
        assert nr["steps"]==ur["steps"] and sha(dest/"initial_flow.csv")==in_wire and sha(official/"initial_flow.csv")==in_uniform
        in_wire=nr["outputs"]["restart.csv"];in_uniform=terminal(official)["outputSHA256"]["restart.csv"]
        steps.append(dict(steps=nr["steps"],wireFrame=in_wire,uniformFrame=in_uniform,wireWindow=nr["last20"],uniformWindow=ur["last20Residual"],
                          wireSeconds=nr["seconds"],uniformSeconds=ur["seconds"],wirePeakBytes=nr["familyPeak"]))
    assert in_wire==s["sourceFrameSHA256"] and in_uniform==u["sourceFrameSHA256"]
    # 全域起步时，只有真实壁面的强标量改变；六个原物理分量末窗按各自回执重新核对。
    with np.load(prep/"nodes.npz") as z:
        xyz,wall,near,pi,ci=(z[k] for k in ("xyz","wallMap","originalFirstMap","overlapParent","overlapChild"))
    wet=np.unique(wall);dry=np.ones(len(xyz),bool);dry[wet]=False
    ref_start=np.loadtxt(base/"restart.csv",delimiter=",",skiprows=1)
    wire_start=np.loadtxt(prep/"initial_flow.csv",delimiter=",",skiprows=1)
    np.testing.assert_array_equal(wire_start[:,:9],ref_start[:,:9]);np.testing.assert_array_equal(wire_start[dry,9],ref_start[dry,9])
    def close(actual,expected,tol):np.testing.assert_allclose(actual,expected,rtol=0,atol=tol)
    for record in (s["window"],s["plus20Window"],u["residualWindow"],u["plus20Window"]):assert len(record)==6 and max(record.values())<-8
    assert s["relativeChange20"]<.002 and u["relativeChange20"]<.002 and abs(s["energy"]["relativeDifference"])<.005 and abs(u["energy"]["relativeDifference"])<.005
    # 用实验式（13）独立复算双线根部热损，系数按原文有效数字，与式（6）至（8）相互核对。
    p=s["profile"];d=model["wireDiameterM"];visc=model["mu"]/model["rho"];wire_k=25.2*4184/3600
    root_kcal_h=1.332e-4*(model["k"]*3600/4184)**.5/visc**.18*model["speedMps"]**.18*12
    root_w=root_kcal_h*4184/3600;assert abs(root_w/(2*p["oneLeadRootW"])-1)<.0001
    exact_mu=p["oneLeadRootW"]/(12*wire_k*np.pi*d**2/4)
    tri=np.array(body["V"])[np.array(body["F"])-1];group=np.array(body["regionID"])
    root_x=abs(tri[group==6,:,0]).min();length=abs(tri[group==7,:,0]).max()-root_x
    wall_delta=12*np.exp(-exact_mu*np.maximum(0,np.abs(xyz[wet,0]*model["diameterM"])-root_x))
    close(wire_start[wet,9],wall_delta/24,1e-13);close(293.15+12*np.exp(-exact_mu*length),p["endK"],1e-12)
    def areas(triangles):return np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    small=xyz[wall]*model["diameterM"];sa=areas(small);area=areas(tri)
    weights=independent_overlap(tri,small,pi,ci)*sa[ci,None];lump=np.bincount(wall.ravel(),np.repeat(sa/3,3),minlength=len(xyz))
    def project(v):return np.bincount(pi,(weights*v[wall[ci]]).sum(axis=1),minlength=len(area))/area
    def recover(stage,reaction):
        with (reaction/"restart.csv").open(encoding="ascii") as f:
            keys=next(csv.reader(f));rhs=np.loadtxt(f,delimiter=",",usecols=keys.index("Residual_Species_0"))
        state=np.loadtxt(stage/"restart.csv",delimiter=",",skiprows=1,usecols=(4,9))
        flux=np.zeros(len(xyz));flux[wet]=rhs[wet]*model["rho"]*3*model["cp"]*24*model["diameterM"]**2/lump[wet]
        return project(flux),101325+project(1.18*9*state[:,0]),293.15+24*state[near,1].mean(axis=1),293.15+24*project(state[:,1])
    q,pressure,nearT,tw=recover(w,w/"reaction0")
    uq,up,ut,uniform_tw=recover(matched_stage,matched_stage/"reaction0");close(uniform_tw,305.15,2e-11)
    q20=recover(w/"plus20",w/"reaction20")[0];uq20=recover(matched_stage/"plus20",matched_stage/"reaction20")[0]
    meas=np.genfromtxt(work/"surface/wall.csv",names=True,delimiter=",");ref=np.genfromtxt(base/"surface/wall.csv",names=True,delimiter=",")
    np.testing.assert_array_equal(meas["faceID"],np.arange(len(area)))
    for name,values,tol in (("heatFluxWm2",q,1e-8),("pressurePa",pressure,1e-9),("wallInputK",tw,1e-10),("nearWallCellTemperatureK",nearT,1e-10)):
        close(meas[name],values,tol)
    close(ref["heatFluxWm2"],meas["uniformHeatFluxWm2"],0)
    control=np.genfromtxt(matched_stage/"surface/wall.csv",names=True,delimiter=",")
    close(uq,control["heatFluxWm2"],1e-8);close(pressure,up,.001)
    # 导出仅使用当前快照。按原三角面检查每个值、资格、温度身份、截线及七区与分带功率。
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip();head="window.ROCKET_CASE="
    assert raw.startswith(head) and raw.endswith(";");data=json.loads(raw[len(head):-1])
    m=data["thermistorWire"];fd={f["key"]:f for f in data["fields"]}
    for k,val in (("pressure",pressure/1000),("heatFlux",q/1000),("nearWallCellTemperature",nearT),("wallInput",tw),
                  ("wireHeatFluxChange",(q-ref["heatFluxWm2"])/1000)):close(fd[k]["values"],val,1e-9)
    assert data["resultCondition"]["wallInputType"]=="perFace" and data["resultCondition"]["wallUniformValueK"] is None
    assert m["boundarySHA256"]==s["boundarySHA256"] and m["jointExperimentVerifiedM2"]==0
    assert "wallSolved" not in fd or all(v is None for v in fd["wallSolved"]["values"])
    with (ui/"wall.csv").open(encoding="utf-8-sig",newline="") as f:
        exported=list(csv.DictReader(f))
    assert len(exported)==len(q)==10272
    for i,r in enumerate(exported):
        assert int(r["face_id"])==i+1 and r["source_frame_sha256"]==in_wire and r["uniform_reference_frame_sha256"]==s["baseFrame"]
        assert r["wall_boundary_sha256"]==s["boundarySHA256"] and r["experiment_doi"]==m["experiment"]["doi"]
        assert abs(float(r["pressure"])*1000-pressure[i])<1e-8 and abs(float(r["heat_flux"])*1000-q[i])<1e-7
        assert abs(float(r["wall_input_K"])-tw[i])<1e-8 and abs(float(r["change_from_uniform_kW_m2"])-(q[i]-ref["heatFluxWm2"][i])/1000)<1e-9
    receipt=read(ui/"temperature.json");assert receipt["csvSHA256"]==sha(ui/"temperature.csv") and receipt["thermistorWire"]==m
    with (ui/"temperature.csv").open(encoding="utf-8-sig",newline="") as f:temperature=list(csv.DictReader(f))
    assert len(temperature)==len(q)
    for i,r in enumerate(temperature):assert int(r["face_id"])==i+1 and r["unit"]=="K" and abs(float(r["value"])-tw[i])<1e-9
    segments=cut_values(ui/"section.csv",tri,proof["sectionXM"],(q-ref["heatFluxWm2"])/1000,"wireHeatFluxChange","kW/m²")
    x=np.abs(tri.mean(axis=1)[:,0])-root_x;bands=m["bands"]
    for j,band in enumerate(bands):
        a=length*j/12;b=length*(j+1)/12;use=(x>=a)&(x<b)
        close([band["fromM"],band["toM"]],[a,b],1e-14)
        loss=wire_k*np.pi*d**2/4*12*exact_mu*(np.exp(-exact_mu*a)-np.exp(-exact_mu*b))
        close(band["theorySideW"],loss,1e-13)
        for region,col in ((2,"backSideW"),(6,"frontSideW")):
            mask=use&(group==region);assert mask.sum()==96
            close(-q[mask]@area[mask],band[col],1e-13)
    bead=np.isin(group,[3,4,5]);deficit=-float((q[bead]-uq[bead])@area[bead])
    delta20=-float((q20[bead]-uq20[bead])@area[bead])
    regions=[]
    for i,region in enumerate(s["regions"],1):
        use=group==i;power=-q[use]@area[use];other=-uq[use]@area[use]
        close(power,region["outwardW"],1e-13)
        regions.append(dict(name=region["name"],wireW=power,matchedW=other,differenceW=power-other))
    h0=-float(uq[bead]@area[bead])/(12*model["sourceAreaM2"])
    primary=dict(stages=steps,reactionJointAreaM2=0,areaM2=area.sum(),faces=len(area),twoRootEquation13W=root_w,profile=p,
                 matchedFrame=in_uniform,wireFrame=in_wire,matchedMaxPressureDifferencePa=max(abs(pressure-up)),
                 matchedBeadChangeW=deficit,matchedBeadChangeAfter20W=delta20,uniformHWm2K=h0,wireHWm2K=m["beadHWm2K"],
                 hDifferenceWm2K=deficit/(12*model["sourceAreaM2"]),regions=regions,bands=bands,energy=s["energy"],
                 matchedEnergy=u["energy"],relativeChange20=s["relativeChange20"],peakFace=s["peakFaceID"],peakFluxWm2=min(q),
                 positiveIntoWallFaces=int(sum(q>0)),sectionSegments=len(segments),uiSeconds=ui_receipt["elapsedSeconds"])
    out.mkdir(parents=True,exist_ok=False)
    (out/"comparison.json").write_text(json.dumps(primary,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    with (out/"wire_bands.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(bands[0]));writer.writeheader();writer.writerows(bands)
    draw(primary,out/"wire.png")
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def add(folder,tag,names):
            for n in names:z.write(folder/n,tag+"/"+n)
        add(work,"wire",["completed.json","control.json","prepared/profile.json","surface/summary.json"])
        add(ui,"ui",["debugger_process.json","ui.json","section.csv"])
        paused=Path(proof["pausedNewJob"]);add(paused,"new_ui",["wire_job.json","progress.json"])
        for k in ("first","resume"):
            add(paused,"new_ui",["receipt-"+saved[k]["invocation"]+".json"])
        z.write(work.parent/"resume.json","new_ui/resume.json")
        for dest in [work/"uniform_control",*sorted(work.glob("steady[0-9][0-9]")),*(w/n for n in ("reaction0","plus20","reaction20","step_control","step_observed"))]:
            add(dest,"wire/"+str(dest.relative_to(work)).replace("\\","/"),["case.cfg","native.json","debugger_process.json","history.csv"])
        for dest in sorted(matched.glob("steady[0-9][0-9]")):
            add(dest,"matched/"+dest.name,["case.cfg","engine_process.json","field_audit.json","resources.json","history.csv"])
        for n in ("reaction0","plus20","reaction20"):add(matched_stage/n,"matched/"+n,["case.cfg","history.csv","engine_process.json","resources.json"])
        add(matched_stage,"matched",["surface/summary.json","surface/wall.csv","ledger/energy.json","ledger/observed/step_process.json","ledger/control/step_process.json"])
        build=Path(read(w/"debugger_process.json")["command"][0]).parent
        add(build,"build",["build.json","SU2_CFD_process.json","CIntegration_process.json","wire_process.json","wire_observer_process.json"])


def draw(data,path):
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    bands=data["bands"];xs=[1e6*b["fromM"] for b in bands]+[300];centers=(np.array(xs[:-1])+xs[1:])/2
    fig,(top,bottom)=plt.subplots(2,1,figsize=(7.2,5.8),sharex=True,constrained_layout=True)
    x=np.linspace(0,.0003,200);top.plot(x*1e6,293.15+12*np.exp(-data["profile"]["muPerM"]*x),color="#246382",label="长引线模型给定壁温")
    top.axhline(305.15,color="#777777",ls="--",label="均温控制");top.set_ylabel("壁温（K）");top.legend(frameon=False)
    for key,label,color,style in (("theorySideW","一维长线","#555555","--"),("backSideW","后引线三维线侧","#197470","-"),("frontSideW","前引线三维线侧","#ac622d","-")):
        bottom.plot(centers,[1e6*b[key] for b in bands],color=color,ls=style,marker=".",label=label)
    bottom.set(xlabel="距名义线根的轴向距离（µm）",ylabel="每 25 µm 线段放热（µW）")
    bottom.legend(frameon=False,ncol=3,fontsize=8);top.set_title("B 型热敏珠 · 3 m/s 侧吹 · 单根引线")
    for ax in (top,bottom):ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False);ax.grid(axis="y",alpha=.18)
    fig.savefig(path,dpi=170);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("work","matched","ui","out"):p.add_argument(n,type=Path)
    a=p.parse_args();verify(*(v.resolve() for v in (a.work,a.matched,a.ui,a.out)))
