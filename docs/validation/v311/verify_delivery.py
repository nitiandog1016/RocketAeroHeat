"""从全周原生节点、独立壁反力和三角面切线重算交付的云图、热量与周向曲线。"""
import argparse
import csv
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from readCM011TunnelSector import vtk_tokens
from thermistorBead import read, sha
from thermistorReaction import terminal


def cut_values(path,triangle,x,values,key,unit):
    """用真实面顶点独立切平面，逐端点核对浏览器截出的线段和两侧方位。"""
    lo=triangle[:,:,0].min(axis=1);hi=triangle[:,:,0].max(axis=1)
    expected=set(np.flatnonzero((lo<x)&(hi>x)));records=[]
    with path.open(encoding="utf-8-sig",newline="") as stream:
        for row in csv.DictReader(stream):
            idx=int(row["face_id"])-1;expected.remove(idx)
            assert row["field_key"]==key and row["unit"]==unit and float(row["section_x_m"])==x
            assert row["availability"]=="研究计算" and row["qualification"]=="未验证"
            ends=np.array([[float(row[v+"_y_m"]),float(row[v+"_z_m"])] for v in ("a","b")])
            pred=[];pts=triangle[idx]
            for a,b in zip(pts,np.roll(pts,-1,axis=0)):
                if (a[0]-x)*(b[0]-x)<0:pred.append((a+(b-a)*((x-a[0])/(b[0]-a[0])))[1:])
            assert len(pred)==2
            if np.linalg.norm(ends[0]-pred[1])<np.linalg.norm(ends[0]-pred[0]):pred.reverse()
            np.testing.assert_allclose(pred,ends,atol=1e-15,rtol=0)
            actual=float(row["value"]);assert abs(actual-values[idx])<1e-10
            mid=ends.mean(axis=0);theta=float(np.rad2deg(np.arctan2(mid[1],mid[0]))%360)
            records.append((theta,idx+1,actual))
    assert not expected and len(records)==96
    return sorted(records)


def plot(profile,out):
    """两幅图共用截面与角度，压力取表压，热流的负号表示试件向空气放热。"""
    plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","DejaVu Sans"]
    fig,axes=plt.subplots(2,1,figsize=(8.2,5.4),sharex=True,constrained_layout=True)
    fig.suptitle("带线热敏珠侧吹  ·  X = 0.172 mm",fontsize=15)
    for ax,rows,label,color in zip(axes,profile,("相对环境壁压（Pa）","入壁热流（kW/m²）"),("#355876","#a04b33")):
        theta,ids,val=np.asarray(rows).T
        ax.plot(theta,val,color=color,lw=1.25,marker=".",ms=3)
        ax.axvline(180,lw=.8,color="#969ba1",ls="--");ax.grid(axis="y",alpha=.22)
        ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False);ax.set_ylabel(label)
    axes[-1].set(xticks=[0,60,120,180,240,300,360],xlim=(0,360),xlabel="周向角（°）    0° 为迎风 +Y，180° 为背风 −Y")
    fig.savefig(out,dpi=165);plt.close(fig)


def verify(run,ui,out):
    run=run.resolve();ui=ui.resolve()
    rec=terminal(run);a=read(run/"field_audit.json");web=read(ui/"ui.json");u=read(ui/"debugger_process.json")
    assert u["status"]=="terminal" and u["debuggerExitCode"]==0 and all(u["inputsStable"].values())
    assert not (u["stopReason"] or u["supervisorError"] or u["survivingProcessIDs"]) and web["faces"]==10272 and web["webglError"]==0 and web["visibleSamples"]>0
    assert web["sourceFrame"]==a["lastFrame"]==sha(run/"restart.csv")
    case=read(run/"case.json");body=read(run/"body.json");s=read(run/"surface/summary.json");ledger=read(run/"ledger/energy.json")
    assert ledger==s["energy"] and s["sourceFrameSHA256"]==web["sourceFrame"] and s["wallSHA256"]==sha(run/"surface/wall.csv")
    arch=ROOT/"cases/thermistor_b/sidewind3.zip";assert sha(arch)==read(arch.with_suffix(".json"))["archiveSHA256"]
    with ZipFile(arch) as package:
        for p in (run/"body.json",run/"case.json",run/"surface/summary.json",run/"surface/wall.csv"):
            assert package.read(p.name)==p.read_bytes()
    f=np.asarray(body["F"])-1;tri=np.asarray(body["V"])[f];groups=np.asarray(body["regionID"])
    n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(n,axis=1)/2;n/=2*area[:,None]
    src=np.loadtxt(run/"restart.csv",delimiter=",",skiprows=1)
    assert terminal(run/"reaction0")["inputSHA256"][str(run/"reaction0/initial_flow.csv")]==a["lastFrame"]
    with (run/"reaction0/restart.csv").open(encoding="ascii") as stream:
        names=next(csv.reader(stream));rhs=np.loadtxt(stream,delimiter=",",usecols=(0,names.index("Residual_Species_0")))
    np.testing.assert_array_equal(rhs[:,0],np.arange(len(src)))
    with np.load(run/"nodes.npz") as z:xyz,wall,first=(z[k] for k in ("xyz","wallMap","firstMap"))
    speed=case["speedMps"];dt=case["wallK"]-case["airK"];scale=case["rho"]*speed*case["cp"]*dt*case["diameterM"]**2
    weight=np.zeros(len(xyz))
    for corner in wall.T:np.add.at(weight,corner,area/3)
    q=(rhs[wall,1]/weight[wall]).mean(axis=1)*scale
    gauge=src[wall,4].mean(axis=1)*case["rho"]*speed**2
    firstK=case["airK"]+dt*src[first,9].mean(axis=1)
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    d=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={v["key"]:v for v in d["fields"]}
    display=dict(d["thermistorSidewind"])
    for k in ("residualWindow","plus20Window"):
        assert display[k]=={n.replace("[","_").replace("]","_").replace("*","_"):v for n,v in s[k].items()}
        display[k]=s[k]
    assert display==s and d["resultCondition"]["physics"]["velocityMps"]==[0,-3,0]
    for k,v in (("pressure",(case["airPa"]+gauge)/1000),("pressureGauge",gauge),("nearWallCellTemperature",firstK),("heatFlux",q/1000),("wallInput",299.15)):
        np.testing.assert_allclose(fields[k]["values"],v,rtol=0,atol=1e-10)
        assert all(i==(2 if k=="wallInput" else 1) for i in fields[k]["status"])
    with (ui/"wall.csv").open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==len(f)==10272
    for j,row in enumerate(rows):
        assert int(row["face_id"])==j+1 and row["joint_verified"]=="0" and row["region"]==body["regionNames"][groups[j]-1]
        assert row["heat_flux_unit"]=="kW/m²" and row["pressure_unit"]=="kPa" and abs(float(row["area_m2"])-area[j])<1e-18
        assert abs(float(row["heat_flux"])*1000-q[j])<1e-8 and abs(float(row["pressure"])*1000-(case["airPa"]+gauge[j]))<1e-9
    release=-area*q;assert abs(sum(release)-ledger["wallReleaseW"])<1e-12
    regionW=[float(release[groups==i].sum()) for i in range(1,8)]
    np.testing.assert_allclose(regionW,[v["outwardW"] for v in s["regions"]],rtol=0,atol=1e-13)
    assert abs(ledger["relativeDifference"])<.005 and s["relativeChange20"]<.002 and max(s["residualWindow"].values())<-8
    p=cut_values(ui/"pressureGauge_section.csv",tri,web["sectionX"],gauge,"pressureGauge","Pa")
    h=cut_values(ui/"heatFlux_section.csv",tri,web["sectionX"],q/1000,"heatFlux","kW/m²")
    assert [v[:2] for v in p]==[v[:2] for v in h]
    sides=[]
    for angle in (0,180):
        indices=[i for i,v in enumerate(p) if v[0]>=angle]
        j=min(indices,key=lambda i:p[i][0]-angle)
        sides.append(dict(faceID=p[j][1],thetaDeg=p[j][0],gaugePa=p[j][2],intoWallWm2=1000*h[j][2]))
    # 同一周向偏角的迎风和背风两面须在同一轴向截面且方位相对。
    assert abs(sides[1]["thetaDeg"]-sides[0]["thetaDeg"]-180)<1e-10 and sides[0]["gaugePa"]>sides[1]["gaugePa"]
    assert sides[0]["intoWallWm2"]<sides[1]["intoWallWm2"]<0
    peak=int(np.argmin(q));assert web["peakFaceID"]==s["peakFaceID"]==peak+1
    mesh_nodes,flow=vtk_tokens(run/"surface.vtk");ids=np.unique(wall)
    dist,idx=cKDTree(mesh_nodes).query(xyz[ids]);assert len(np.unique(idx))==len(ids)==len(mesh_nodes) and dist.max()<2e-6
    cf=flow["Skin_Friction_Coefficient"][idx[np.searchsorted(ids,wall)]]
    first_height=case["diameterM"]*abs(np.einsum("fij,fj->fi",xyz[first]-xyz[wall],n)).mean(axis=1)
    yplus=first_height*case["rho"]*speed/case["mu"]*np.sqrt(np.linalg.norm(cf,axis=2).mean(axis=1)/2)
    assert yplus.min()>=0 and np.isfinite(yplus).all() and first_height.min()>0
    out.mkdir(parents=True,exist_ok=False)
    def limits(v):return [float(np.min(v)),float(np.max(v))]
    result=dict(faces=len(f),areaM2=float(np.sum(area)),regions=s["regions"],totalOutwardW=float(sum(release)),
                absolutePressurePa=limits(case["airPa"]+gauge),intoWallWm2=limits(q),firstLayerK=limits(firstK),
                firstHeightM=limits(first_height),estimatedFirstYPlus=limits(yplus),peakFaceID=peak+1,
                peakRegion=body["regionNames"][groups[peak]-1],peakCenterM=tri[peak].mean(axis=0).tolist(),
                sectionXM=web["sectionX"],segments=len(p),windwardAndLee=sides,jointExperimentVerifiedM2=0,
                ledger=ledger,relativeChange20=s["relativeChange20"],terminalSteps=a["cumulativeSteps"],
                residualWindow=s["residualWindow"],meshAuditSHA256=a["meshAuditSHA256"],sourceFrameSHA256=a["lastFrame"],
                guiSeconds=u["elapsedSeconds"],fields=sha(run/"surface/wall.csv"),fullExport=sha(ui/"wall.csv"),archive=sha(arch))
    (out/"verification.json").write_bytes((json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))
    plot((p,h),out/"circumference.png")
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        def add(folder,names,prefix):
            for name in names:z.write(folder/name,prefix+"/"+name)
        add(ui,("pressureGauge_section.csv","heatFlux_section.csv","ui.json","debugger_process.json"),"workbench")
        add(run.parent/"cold01",("inputs.json","mesh_audit.json","case.json","case.cfg"),"cold")
        for part in ("steady01","steady02"):
            add(run.parent/part,("engine_request.json","engine_process.json","field_audit.json","resources.json","history.csv","case.cfg"),part)
        for part in ("reaction0","reaction20","plus20"):
            add(run/part,("engine_process.json","resources.json","case.cfg"),part)
        for part in ("control","observed"):add(run/"ledger"/part,("step_process.json",),"ledger/"+part)
        add(run/"ledger",("energy.json",),"ledger")
        # 首段执行时，续算分支还未追加 STAGES.md 指纹；保存能复现两个真实启动指纹的源码。
        current=(HERE/"run_steady.py").read_bytes()
        import hashlib
        for part in ("steady01","steady02"):
            receipt=read(run.parent/part/"engine_process.json")
            path=str(HERE/"run_steady.py")
            version=current if part=="steady02" else current.replace(b'HERE/"STAGES.md",parent/"engine_process.json"',b'parent/"engine_process.json"')
            assert hashlib.sha256(version).hexdigest()==receipt["inputSHA256"][path]
            z.writestr(part+"/run_steady.py",version)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("run","ui","out"):parser.add_argument(key,type=Path)
    a=parser.parse_args();verify(a.run,a.ui,a.out)
