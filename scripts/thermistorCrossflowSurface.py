"""核对 B 型热敏珠侧吹的末窗、同态壁反力及独立热量账，再输出真实方位上的逐面结果。"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import sys
from zipfile import ZipFile

import numpy as np

from beadSurfaceResult import thermal_frame
from thermistorBead import ENGINE, ENGINE_SHA, ROOT, read, sha, write_json
from thermistorReaction import observe, terminal

sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def reaction(parent,probe,wall_ids,wall_value=1.0):
    """使用同一接受态的官方释放壁行读数，同时核实独立目录的精确配置。"""
    rec,_=thermal_frame(parent);terminal(probe)
    assert sha(probe/"initial_flow.csv")==rec["outputSHA256"]["restart.csv"]
    assert sha(probe/"mesh.su2")==sha(parent/"mesh.su2")
    cfg=re.sub(r"^ITER= \d+$","ITER= 1",(parent/"case.cfg").read_text("ascii"),flags=re.M)
    value=f"(wall, VALUE, {float(wall_value)})";assert cfg.count(value)==1
    cfg=cfg.replace(value,"(wall, FLUX, 0.0)").replace("SOLUTION, PRIMITIVE","SOLUTION, RESIDUAL")
    assert (probe/"case.cfg").read_text("ascii")==cfg+"WRT_RESTART_COMPACT= NO\n"
    with (probe/"restart.csv").open(encoding="ascii") as f:
        header=next(csv.reader(f))
        raw=np.loadtxt(f,delimiter=",",usecols=(0,header.index("Residual_Species_0")))
    np.testing.assert_array_equal(raw[:,0],np.arange(len(raw)))
    assert np.isfinite(raw).all() and np.all(raw[wall_ids,1]<0)
    return raw[:,1]


def observer_binary(work):
    """复用已发布指纹对应的观察器；独立安装可用项目内固定源码重建并比对构建记录。"""
    with ZipFile(ROOT/"docs/validation/v305/evidence/receipts.zip") as z:
        source=json.loads(z.read("build/build.json"))
    exe=Path(os.environ.get("RAH_THERMISTOR_LAYER_OBSERVER",ROOT/"outputs/v305_validation/build02/bead_layers.exe"))
    if exe.is_file():
        assert sha(exe)==source["executableSHA256"]
        return exe
    sys.path.insert(0,str(ROOT/"docs/validation/v305"))
    from observe_layer import build
    build(work)
    result=read(work/"build.json");exe=work/"bead_layers.exe"
    assert sha(exe)==result["executableSHA256"]
    for key in ("compilerSHA256","sourceArchiveSHA256","sourceSHA256","modifiedSHA256","cachedInputs"):
        assert result[key]==source[key]
    return exe


def energy(parent,out,source_theta,wall_ids,wall_rows,scale,exe,*,angular96=False,temperature_scale=6):
    """对比官方一步以及同阶段的壁面与远场收支。"""
    assert sha(ENGINE)==ENGINE_SHA
    cfg,n=re.subn(r"^ITER= \d+$","ITER= 1",(parent/"case.cfg").read_text("ascii"),flags=re.M);assert n==1
    out.mkdir(exist_ok=False)
    for name,binary in (("control",ENGINE),("observed",exe)):
        run=out/name;run.mkdir();os.link(parent/"mesh.su2",run/"mesh.su2");os.link(parent/"restart.csv",run/"initial_flow.csv")
        (run/"case.cfg").write_text(cfg,encoding="ascii")
        invoke("step",[binary,"case.cfg","--threads","1"],run,run,150 if angular96 else 90,
               4 if angular96 else 3,trace=name=="observed")
        rec=read(run/"step_process.json");assert rec["exitCode"]==0 and not (rec["stopReason"] or rec["livePids"])
        assert sum(p.stat().st_size for p in run.iterdir() if p.is_file())<(1536 if angular96 else 1024)*1024**2
    end=[np.loadtxt(out/name/"restart.csv",delimiter=",",skiprows=1) for name in ("control","observed")]
    delta=float(np.max(abs(end[0]-end[1])));assert delta<2e-7
    stages={}
    for name in ("start","beforeFar","afterFar","beforeWall","afterWall","end"):
        s=np.loadtxt(out/"observed"/(name+".csv"),delimiter=",",skiprows=1);s=s[np.argsort(s[:,0])]
        np.testing.assert_array_equal(s[:,0],np.arange(len(source_theta)));np.testing.assert_array_equal(s[:,1],source_theta)
        stages[name]=float(s[:,2].sum()*scale)
        if name=="beforeWall":error=float(np.max(abs(wall_rows[wall_ids]+s[wall_ids,2]))*scale)
    wall=stages["beforeWall"]-stages["afterWall"];far=stages["afterFar"]-stages["beforeFar"]
    assert error<1e-11 and wall>0 and far>0 and abs(stages["beforeFar"]-stages["start"])<1e-11
    assert abs(wall+wall_rows[wall_ids].sum()*scale)<1e-10
    return dict(wallReleaseW=wall,farOutW=far,relativeDifference=(wall-far)/wall,stageW=stages,
                maxWallNodeRowDifferenceW=error,maxOneUpdateDifference=delta,observerSHA256=sha(exe),
                inputFrameSHA256=sha(parent/"restart.csv"),temperatureScaleK=temperature_scale,domainCopies=1)


def collect(parent,exe=None):
    rec,residual=thermal_frame(parent)
    info=read(parent/"field_audit.json");assert info["allResidualsBelowMinus8"] and info["lastFrame"]==sha(parent/"restart.csv")
    c=read(parent/"case.json");body=read(parent/"body.json")
    assert c["caseKind"]=="thermistor-B-sidewind-two-leads" and c["velocityMps"]==[0,-3,0] and c["sectorDeg"]==360
    for name in ("case.json","nodes.npz","body.json","mesh.su2"):
        assert sha(parent/name)==rec["inputSHA256"][str(parent/name)]
    delta=c["wallK"]-c["airK"];temperature_scale=c.get("scalarReferenceDeltaK",delta)
    wall_value=delta/temperature_scale
    assert (delta,temperature_scale,wall_value) in ((6,6,1),(12,24,.5))
    scale=c["rho"]*c["speedMps"]*c["cp"]*temperature_scale*c["diameterM"]**2
    angular96=c["angularCells"]==96
    with np.load(parent/"nodes.npz") as z:
        wall,near,area,xyz=(z[k] for k in ("wallMap","firstMap","area","xyz"))
        if angular96:
            original_area,near,pi,ci,weights=(z[k] for k in ("originalArea","originalFirstMap","overlapParent","overlapChild","shapeAreaM2"))
            assert c["reportAngularCells"]==48 and len(original_area)==len(body["F"])
    wall_ids=np.unique(wall)
    budget=dict(maximum_seconds=180 if angular96 else 120,memory_bytes=(4 if angular96 else 3)*1024**3,
                output_bytes=(850 if angular96 else 450)*1024**2)
    for name,source,release in (("reaction0",parent,True),("plus20",parent,False),("reaction20",parent/"plus20",True)):
        dest=parent/name
        if not dest.exists():observe(source,dest,release=release,wall_value=wall_value,**budget)
        terminal(dest)
    plus=parent/"plus20"
    assert (plus/"case.cfg").read_text("ascii")==re.sub(r"^ITER= \d+$","ITER= 20",(parent/"case.cfg").read_text("ascii"),flags=re.M)
    assert sha(plus/"initial_flow.csv")==info["lastFrame"]
    raw=reaction(parent,parent/"reaction0",wall_ids,wall_value);other=reaction(plus,parent/"reaction20",wall_ids,wall_value)
    _,plus_res=thermal_frame(plus)
    for src in (plus,parent):
        state=np.loadtxt(src/"restart.csv",delimiter=",",skiprows=1)
        np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)))
        np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-13)
        assert np.isfinite(state).all() and state[:,9].min()>=-1e-12 and state[:,9].max()<=wall_value+1e-12
        assert np.all(state[wall_ids,9]==wall_value) and np.max(abs(state[wall_ids,5:8]))<1e-10 and np.min(c["airPa"]+c["rho"]*9*state[:,4])>0
    power=-scale*raw[wall_ids].sum();change=abs(other[wall_ids].sum()/raw[wall_ids].sum()-1)
    assert change<.002
    if exe is None:exe=observer_binary(parent/"observer_build")
    ledger=energy(parent,parent/"ledger",state[:,9],wall_ids,raw,scale,exe,angular96=angular96,temperature_scale=temperature_scale)
    write_json(parent/"ledger/energy.json",ledger);assert abs(ledger["relativeDifference"])<.005
    # 无扇区反射，每个真实壁节点按相邻三角面面积三分之一分摊功率。
    lump=np.bincount(wall.ravel(),weights=np.repeat(area/3,3),minlength=len(xyz))
    if angular96:
        # 先按真实细三角面积分摊节点反力，再积分细面形函数到原报告三角形。
        from thermistorAngularSpacing import project
        density=np.zeros(len(xyz));density[wall_ids]=scale*raw[wall_ids]/lump[wall_ids]
        q=project(density,wall,original_area,pi,ci,weights)
        gauge=project(c["rho"]*c["speedMps"]**2*state[:,4],wall,original_area,pi,ci,weights)
        area=original_area
    else:
        q=(scale*raw[wall]/lump[wall]).mean(axis=1)
        gauge=c["rho"]*c["speedMps"]**2*state[wall,4].mean(axis=1)
    temperature=c["airK"]+temperature_scale*state[near,9].mean(axis=1)
    watts=-area*q;assert abs(watts.sum()-power)<1e-12
    group=np.asarray(body["regionID"]);regions=[]
    for i,name in enumerate(body["regionNames"],1):
        use=group==i;regions.append(dict(name=name,faces=int(use.sum()),areaM2=float(area[use].sum()),outwardW=float(watts[use].sum())))
    out=parent/"surface";out.mkdir(exist_ok=False)
    table=out/"wall.csv"
    np.savetxt(table,np.column_stack((np.arange(len(q)),group,area,c["airPa"]+gauge,gauge,q,temperature)),delimiter=",",comments="",
               fmt=["%d","%d"]+["%.17g"]*5,header="faceID,regionID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,nearWallCellTemperatureK")
    summary=dict(faces=len(q),areaM2=float(area.sum()),regions=regions,totalOutwardW=float(power),beadAndJointsW=float(watts[(group>=3)&(group<=5)].sum()),
                 peakFaceID=int(np.argmin(q))+1,relativeChange20=float(change),residualWindow=residual,plus20Window=plus_res,
                 energy=ledger,jointExperimentVerifiedM2=0,meshSHA256=sha(parent/"mesh.su2"),bodySHA256=sha(parent/"body.json"),
                 caseSHA256=sha(parent/"case.json"),wallSHA256=sha(table),sourceFrameSHA256=info["lastFrame"],
                 sourceReactionSHA256=sha(parent/"reaction0/restart.csv"),cumulativeSteps=info["cumulativeSteps"])
    if angular96:
        sources=[Path(__file__).resolve(),ROOT/"scripts/thermistorAngularSpacing.py",ROOT/"scripts/thermistorReaction.py",
                 ROOT/"scripts/beadSurfaceResult.py",ROOT/"docs/validation/v272/observe_flux.py"]
        summary.update(nativeWallTriangles=len(wall),angularCells=96,reportAngularCells=48,
                       wallMapping="P1-overlap-integral",temperatureSampling="v313-common-first-nodes",
                       nodeMapSHA256=sha(parent/"nodes.npz"),plus20FrameSHA256=sha(plus/"restart.csv"),
                       reaction20SHA256=sha(parent/"reaction20/restart.csv"),
                       postprocessSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    write_json(out/"summary.json",summary)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("run",type=Path)
    args=p.parse_args();folder=args.run.resolve();folder.relative_to(ROOT/"outputs");collect(folder)
