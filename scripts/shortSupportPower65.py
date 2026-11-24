"""以外露银球功率为边界，分别求解两个括区工况和反求壁温处的独立原生场。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[variable]="1"
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"docs/validation/v329"))
from run_shared import advance, check, read, sha, table, write_json, DRIVER
from map_surface import CASE, DATA, mapped, project_surface, residual
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive
from shortSupportUniformDriver import build as build_uniform

PLAN=ROOT/"docs/validation/v340/PLAN.md"
BASE=ROOT/"outputs/v332_validation"
PREP_HASH="22ef3bf27abf1a093fae6bfb03145fac1b9f401229f11a3d12042be19cd9a04c"
BASE_FRAME="c23c9ca1beb024105c600642f3fb016c9b7077599c8cecff227eb6d3ab0dbe63"
TRANSFER=DATA/"supportedAngularTransfer65.npz"


def sources(job):
    """核对稳态父场及驱动。新安装从现有项目入口建立相同试件的本地依赖。"""
    prep=BASE/"prepared01";driver=BASE/"build01";steady=BASE/"run01";release=DRIVER/"wire.exe"
    legacy=all(p.is_file() for p in (prep/"prepared.json",driver/"build.json",driver/"shared_process.json",steady/"thermal.json",release))
    if legacy:
        legacy=Path(read(driver/"shared_process.json")["args"][2])==driver/"shared.exe"
    if legacy:
        assert sha(driver/"build.json")=="0bbcde57b7bc0aa2889bf1426613a002c498da980f61e06b67b8ae2db7e94c19"
        assert sha(steady/"thermal.json")=="746119a21a0c9f6902e9cbe331bd22c18dbc2d6a3ea3060906a63afe35e6657d"
        assert sha(release)=="9721e25124a866652e8d004b531854ec88d2f856e99ca7b80c0ab34d2bbebb55"
        assert read(steady/"converged.json")["frame"]==BASE_FRAME
        # 稀疏静态库引用独立目标文件；链接前将这些输入和原编译参数与归档核对。
        ref=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
        files=ref["cachedBuildInputs"]
        old=next(Path(k) for k in files if Path(k).as_posix().endswith("/source/SU2_CFD/src/integration/CIntegration.cpp")).parents[4]
        sdk=ROOT/"tools/runtime/su2-thermal-build"
        assert all(sha(sdk/Path(k).relative_to(old))==digest for k,digest in files.items())
        with ZipFile(ROOT/"docs/validation/v332/evidence/receipts.zip") as records:
            for file in ("build.json","CSpeciesSolver_process.json","CIncEulerSolver_process.json","shared_process.json"):
                assert records.read("build01/"+file)==(driver/file).read_bytes()
    else:
        from shortSupportFirstRun import inputs
        prep,native,steady,_=inputs(job);driver=native/"steady";release=driver/"wire.exe"
        assert sha(release)==read(driver/"build.json")["files"]["wire.exe"]
    build=read(driver/"build.json");engine=driver/"shared.exe"
    assert sha(prep/"prepared.json")==build["preparedSHA256"]==PREP_HASH and sha(engine)==build["exe"]
    spec=read(prep/"prepared.json")
    assert all(sha(prep/name)==digest for name,digest in spec["files"].items()) and sha(CASE)==spec["caseSHA256"]
    assert all(sha(driver/name)==digest for name,digest in build["files"].items())
    conv=read(steady/"converged.json");last=steady/conv["stage"];state=last/"restart_000199.csv"
    thermal=read(steady/"thermal.json");rec=read(last/"record.json");proc=read(last/"debugger_process.json")
    assert thermal["passed"] and sha(state)==conv["frame"]==thermal["early"]["sourceFrame"]==rec["end"]
    assert max(rec["window"]) < -8 and rec["rawMax20"]<1e-8 and rec["exe"]==sha(engine)
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["survivingProcessIDs"] or proc["supervisorError"])
    assert sha(TRANSFER)==read(DATA/"supportedAngular65.json")["transferSHA256"]
    write_json(job/"source.json",dict(prepared=str(prep),preparedSHA256=PREP_HASH,frame=str(state),frameSHA256=sha(state),
        buildSHA256=sha(driver/"build.json"),engineSHA256=sha(engine),releaseSHA256=sha(release),caseSHA256=sha(CASE),
        transferSHA256=sha(TRANSFER),thermalSHA256=sha(steady/"thermal.json"),planSHA256=sha(PLAN),
        uniformContractSHA256=sha(PLAN.with_name("NATIVE.md"))))
    return prep,state,engine,release


def prepare(directory,base,frame,excess):
    """为本次试算写独立温界、面映射中的壁温和全体初值；原网和几何权重只读复用。"""
    directory.mkdir()
    spec=read(base/"prepared.json");condition=spec["condition"]
    theta=excess/condition["scalarUnitK"];wall=condition["airK"]+excess
    for name in ("mesh.su2","transfer.npz"):os.link(base/name,directory/name)
    with np.load(base/"mapping.npz") as grid:data=dict(grid)
    wet=data["wet"];previous=data["theta"];data["theta"]=np.full_like(previous,theta)
    np.savez_compressed(directory/"mapping.npz",**data)
    for name,scalar in (("profile.dat",theta),("release.dat",0.)):
        np.savetxt(directory/name,np.column_stack((wet,np.full(len(wet),1+theta),np.full(len(wet),scalar))),fmt=["%d","%.17g","%.17g"])
    values=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=range(10))
    np.testing.assert_array_equal(values[:,1:4],data["xyz"]);np.testing.assert_allclose(values[wet,9],previous,rtol=0,atol=2e-14)
    values[:,9]*=theta/max(previous)
    values[wet,9]=theta
    values[wet,8]+=theta-previous
    assert np.isfinite(values).all() and values[:,9].min()>=-1e-12 and values[:,9].max()<=theta+1e-12
    with frame.open("r",encoding="ascii") as f:header=f.readline().strip()
    np.savetxt(directory/"initial.csv",values,fmt=["%d"]+["%.17g"]*9,delimiter=",",header=header,comments="")
    spec.update(equatorK=wall,temperatureBoundary="uniform-calorimetry-trial",basePreparedSHA256=PREP_HASH,
                parentFrame=sha(frame),uniformExcessK=excess)
    spec["files"]={name:sha(directory/name) for name in spec["files"]}
    write_json(directory/"prepared.json",spec)


def second_surface(work,prep,dest,thermal):
    """在追加二十次更新后的真实端点重新积分五区功率和峰值。"""
    with np.load(prep/"mapping.npz") as m:wet=m["wet"]
    n=np.fromfile(work/"late/nodes.bin","<f8").reshape(-1,16);n=n[np.argsort(n[:,0])]
    state=table(work/"late","initial_flow.csv")
    assert sha(work/"late/initial_flow.csv")==thermal["late"]["sourceFrame"]
    np.testing.assert_array_equal(n[:,1],state[:,9])
    row=residual(work/"releaseLate/restart.csv")
    nodeW=(row[wet]+n[wet,9]-n[wet,10])*thermal["scale"]
    return project_surface(dest,prep,TRANSFER,sha(TRANSFER),state,nodeW,thermal["late"]["wallW"])


def trial(job,name,excess,inputs):
    base,frame,engine,release=inputs
    directory=job/name;directory.mkdir();prep=directory/"prepared";prepare(prep,base,frame,excess)
    engine=build_uniform(directory/"driver",engine.parent,excess/read(prep/"prepared.json")["condition"]["scalarUnitK"])
    work=directory/"run";work.mkdir();seed=prep/"initial.csv";trend=[]
    options=dict(prepared=prep,plan=PLAN)
    for index in range(1,6):
        stage=work/f"steady{index:02d}"
        atomic_write_json(job/"progress.json",dict(stage=name,message=f"{name}，均匀超温 {excess:.6f} K，第 {index} 段"))
        print(name,stage.name,excess,flush=True)
        r=advance(stage,seed,engine,200,**options)
        worst=max(max(r["window"]),np.log10(r["rawMax20"]));trend.append(worst);seed=stage/"restart.csv"
        if worst < -8:break
        if len(trend)>=3 and trend[-1]>=trend[-2]>=trend[-3]:raise RuntimeError("两段连续不改善")
    else:raise RuntimeError("本候选已用尽稳态预算")
    write_json(work/"converged.json",dict(stage=stage.name,frame=sha(seed),trend=trend))
    for label,start,steps,free,capture in (("plus20",seed,20,False,False),("control",seed,1,False,False),
            ("early",seed,1,False,True),("releaseEarly",seed,1,True,False),
            ("late",work/"plus20/restart.csv",1,False,True),("releaseLate",work/"plus20/restart.csv",1,True,False)):
        advance(work/label,start,release if free else engine,steps,free,capture,**options)
    thermal=check(work,prep)
    early=mapped(work,directory/"surface",prepared=prep,expected_mesh=thermal["meshSHA256"],transfer_file=TRANSFER,transfer_sha=sha(TRANSFER))
    late=second_surface(work,prep,directory/"lateSurface",thermal)
    result=dict(excessK=excess,temperatureK=early["condition"]["airK"]+excess,preparedSHA256=sha(prep/"prepared.json"),
                initialSHA256=sha(prep/"initial.csv"),earlyFrame=early["frameSHA256"],lateFrame=thermal["late"]["sourceFrame"],
                earlyRegionW=[r["outW"] for r in early["regions"]],lateRegionW=[r["outW"] for r in late["regions"]],
                earlyPeak=early["peakOutFace"],latePeak=late["peakOutFace"],earlyPeakWm2=early["peakOutWm2"],latePeakWm2=late["peakOutWm2"],
                thermalSHA256=sha(work/"thermal.json"),surfaceSHA256=sha(directory/"surface/surface.zip"))
    write_json(directory/"trial.json",result)
    print(name,"sphere W",result["earlyRegionW"][0],flush=True)
    return result


def solve(job):
    job=job.resolve();job.relative_to(ROOT/"outputs")
    if shutil.disk_usage(ROOT).free<12*1024**3:raise OSError("此工况需要 12 GiB 可用磁盘。")
    job.mkdir(parents=True,exist_ok=False)
    with exclusive(job):
        inputs=sources(job);case=read(CASE)
        target=case["observedBtuS"]*1055.05585262/case["supportAreaFactor"]
        high=trial(job,"high",40.,inputs);low=trial(job,"low",20.,inputs)
        qlo=low["earlyRegionW"][0];qhi=high["earlyRegionW"][0]
        assert 0<qlo<target<qhi and abs(qhi/40/(qlo/20)-1)<.002
        lo=table(job/"low/surface","wall.csv");hi=table(job/"high/surface","wall.csv")
        np.testing.assert_array_equal(lo[:,:3],hi[:,:3])
        air=read(inputs[0]/"prepared.json")["condition"]["airK"]
        shape=dict(maxHeatPerKDifference=float(max(abs(hi[:,5]/40-lo[:,5]/20))),
                   maxFirstRatioDifference=float(max(abs((hi[:,7]-air)/40-(lo[:,7]-air)/20))))
        assert shape["maxHeatPerKDifference"]<.002*max(abs(hi[:,5]/40)) and shape["maxFirstRatioDifference"]<.002
        candidate=20+(target-qlo)*20/(qhi-qlo)
        power=trial(job,"balance",candidate,inputs)
        for key in ("earlyRegionW","lateRegionW"):
            assert abs(power[key][0]/target-1)<.0002
        out=job/"balance/surface";archive=job/"power65.zip"
        with ZipFile(archive,"x") as z:
            for name in ("summary.json","wall.csv","stations.csv"):
                z.writestr(ZipInfo(name),(out/name).read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
        ack=dict(schema="short65-sphere-power-v1",targetSphereW=target,temperatureK=power["temperatureK"],
                 source=read(job/"source.json"),high=high,low=low,balance=power,normalizedFields=shape,
                 archiveSHA256=sha(archive),faces=len(hi),areaM2=sum(hi[:,2]),experimentalJointAreaM2=0.)
        atomic_write_json(job/"completed.json",ack)
        atomic_write_json(job/"progress.json",dict(stage="complete",message="球面功率、壁温和五区原面场已核对"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("job",type=Path)
    solve(parser.parse_args().job)
