"""用强等温标量的同态壁反力读取长引线算例，按原三角形交叠面积返回同一科研出口。"""
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from thermistorAngularSpacing import project
from thermistorBead import read, sha, write_json
from thermistorWireCase import native


def residual(path):
    with path.open(encoding="ascii") as f:
        header=next(csv.reader(f))
        values=np.loadtxt(f,delimiter=",",usecols=(0,header.index("Residual_Species_0")))
    np.testing.assert_array_equal(values[:,0],np.arange(len(values)))
    assert np.isfinite(values).all()
    return values[:,1]


def collect(prep,build,stage,work,checkpoint=None):
    profile=read(prep/"profile.json");c=read(prep/"case.json")
    base=Path(profile["base"]);body=read(prep/"body.json")
    frame=stage/"restart.csv";first=read(stage/"native.json")
    assert first["passed"] and sha(frame)==first["outputs"]["restart.csv"]
    reaction0=stage/"reaction0";plus=stage/"plus20";reaction20=stage/"reaction20"
    control=stage/"step_control";observed=stage/"step_observed"
    phases=((reaction0,frame,1,"release"),(plus,frame,20,"profile"),(reaction20,plus/"restart.csv",1,"release"),
            (control,frame,1,"profile"),(observed,frame,1,"profile"))
    for dest,start,steps,mode in phases:
        if checkpoint:checkpoint(dest.name)
        audit=native(prep,build,dest,start,steps,mode,trace=dest==observed)
        if dest==plus:
            extended=audit
            assert extended["passed"]
    if checkpoint:checkpoint("surface")
    a,b=[np.loadtxt(p/"restart.csv",delimiter=",",skiprows=1) for p in (control,observed)]
    step_delta=float(np.max(abs(a-b)));assert step_delta<2e-7
    del a,b
    with np.load(prep/"nodes.npz") as z:
        wall,near,area,original,pi,ci,weights=(z[k] for k in ("wallMap","originalFirstMap","area","originalArea","overlapParent","overlapChild","shapeAreaM2"))
    state=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=(4,9));wet=np.unique(wall)
    scale=c["rho"]*c["speedMps"]*c["cp"]*24*c["diameterM"]**2
    raw=residual(reaction0/"restart.csv");later=residual(reaction20/"restart.csv")
    power=-raw[wet].sum()*scale
    change=float(abs(later[wet].sum()/raw[wet].sum()-1));assert power>0 and change<.002
    # 同一更新之前的全壁组装行、释放反力和远场净流出，三者独立核对。
    stages={}
    for tag in ("start","beforeFar","afterFar","beforeWall","afterWall","end"):
        rows=np.loadtxt(observed/(tag+".csv"),delimiter=",",skiprows=1)
        rows=rows[np.argsort(rows[:,0])]
        np.testing.assert_array_equal(rows[:,0],np.arange(len(state)));np.testing.assert_array_equal(rows[:,1],state[:,1])
        stages[tag]=float(rows[:,2].sum()*scale)
        if tag=="beforeWall":row_error=float(np.max(abs(rows[wet,2]+raw[wet]))*scale)
    released=stages["beforeWall"]-stages["afterWall"];far=stages["afterFar"]-stages["beforeFar"]
    assert row_error<1e-11 and abs(stages["beforeFar"]-stages["start"])<1e-11 and abs(released-power)<1e-10
    imbalance=(released-far)/released;assert far>0 and abs(imbalance)<.005
    def average(nodes):
        return project(nodes,wall,original,pi,ci,weights)
    lump=np.bincount(wall.ravel(),weights=np.repeat(area/3,3),minlength=len(state))
    def heat(rhs):
        density=np.zeros(len(state));density[wet]=scale*rhs[wet]/lump[wet]
        return average(density)
    q=heat(raw);q20=heat(later);gauge=average(c["rho"]*9*state[:,0])
    temperatures=c["airK"]+24*average(state[:,1])
    adjacent=c["airK"]+24*state[near,1].mean(axis=1)
    # 交叠积分精度内允许常值面的舍入偏差，保留原始计算值。
    assert min(temperatures)>=profile["endK"]-1e-10 and max(temperatures)<=c["wallK"]+1e-10
    watts=-original*q;assert abs(power-watts.sum())<1e-12
    group=np.array(body["regionID"]);reference=np.genfromtxt(base/"surface/wall.csv",names=True,delimiter=",")
    np.testing.assert_array_equal(reference["regionID"],group);np.testing.assert_array_equal(reference["areaM2"],original)
    base_summary=read(base/"surface/summary.json")
    regions=[]
    for i,name in enumerate(body["regionNames"],1):
        sel=group==i
        regions.append(dict(name=name,faces=int(sel.sum()),areaM2=float(original[sel].sum()),outwardW=float(watts[sel].sum()),
                            controlW=base_summary["regions"][i-1]["outwardW"],wallRangeK=[min(temperatures[sel]),max(temperatures[sel])]))
    bead=np.isin(group,[3,4,5]);bead_w=float(watts[bead].sum());base_bead=base_summary["beadAndJointsW"]
    method=read(prep/"profile.json").copy();method.pop("inputs");method.pop("base")
    energy=dict(sourceFrameSHA256=sha(frame),stageW=stages,wallReleaseW=released,farOutW=far,relativeDifference=imbalance,
                maxNodeDifferenceW=row_error,maxOneUpdateDifference=step_delta,observerSHA256=sha(build/"wire_observer.exe"))
    summary=dict(profile=method,sourceFrameSHA256=sha(frame),reactionSHA256=sha(reaction0/"restart.csv"),baseFrame=profile["baseFrame"],
        meshSHA256=profile["meshSHA256"],nodeMapSHA256=sha(prep/"nodes.npz"),bodySHA256=sha(prep/"body.json"),caseSHA256=sha(prep/"case.json"),
        boundarySHA256=sha(prep/"profile.dat"),solverSHA256=sha(build/"wire.exe"),faces=len(q),areaM2=float(original.sum()),
        totalOutwardW=power,beadAndJointsW=bead_w,beadHWm2K=bead_w/(c["sourceAreaM2"]*12),beadChangeW=bead_w-base_bead,
        controlBeadW=base_bead,controlTotalW=base_summary["totalOutwardW"],regions=regions,peakFaceID=int(np.argmin(q))+1,
        energy=energy,relativeChange20=change,beadChange20W=float(-original[bead]@(q20[bead]-q[bead])),
        plus20FrameSHA256=sha(plus/"restart.csv"),reaction20SHA256=sha(reaction20/"restart.csv"),
        window=first["last20"],plus20Window=extended["last20"],jointExperimentVerifiedM2=0)
    names="faceID,regionID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,nearWallCellTemperatureK,wallInputK,uniformHeatFluxWm2"
    with TemporaryDirectory(dir=work) as tmp:
        out=Path(tmp)/"surface";out.mkdir()
        np.savetxt(out/"wall.csv",np.column_stack((np.arange(len(q)),group,original,c["airPa"]+gauge,gauge,q,adjacent,temperatures,reference["heatFluxWm2"])),
                   fmt=["%d","%d"]+["%.17g"]*7,delimiter=",",header=names,comments="")
        summary["wallSHA256"]=sha(out/"wall.csv");write_json(out/"summary.json",summary)
        target=work/"surface"
        if not target.exists():out.rename(target)
        else:assert sha(target/"wall.csv")==summary["wallSHA256"] and read(target/"summary.json")==summary
    completed=dict(stage=stage.name,frame=summary["sourceFrameSHA256"],wallSHA256=summary["wallSHA256"],summarySHA256=sha(target/"summary.json"))
    if (work/"completed.json").exists():assert read(work/"completed.json")==completed
    else:write_json(work/"completed.json",completed)
