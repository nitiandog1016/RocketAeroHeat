"""在同一 B 珠体网施加长引线的一维温度分布，组织实际的三维原生推进和来源核对。"""
import csv
import os
from pathlib import Path

import numpy as np

from runSU2NativeDebugger import observe
from thermistorBead import ROOT, read, sha, write_json
from thermistorMeasuredCold import GRID_KEY, HEADER, RESIDUALS, config
from thermistorReaction import terminal


def prepare(job, out):
    """从已审核的同一 12 K 冷启动读取初猜，只改变引线的给定热边界。"""
    done=read(job/"completed.json")
    assert done["jobSHA256"]==sha(job/"job.json") and done["stage"] in [f"steady{i:02d}" for i in range(1,9)]
    base=job/done["stage"];end=terminal(base);summary=read(base/"surface/summary.json")
    assert done["frameSHA256"]==summary["sourceFrameSHA256"]==end["outputSHA256"]["restart.csv"]
    assert done["summarySHA256"]==sha(base/"surface/summary.json") and done["wallSHA256"]==sha(base/"surface/wall.csv")
    key=read(GRID_KEY)
    for name,field in (("mesh.su2","meshSHA256"),("body.json","bodySHA256"),("nodes.npz","nodeMapSHA256")):
        assert sha(base/name)==key[field]==done[field]
    c=read(base/"case.json");body=read(base/"body.json")
    assert c["speedMps"]==3 and c["velocityMps"]==[0,-3,0] and c["airK"]==293.15 and c["wallK"]==305.15
    assert c["wireDiameterM"]==54e-6 and c["exposedStraightLengthM"]==.000300 and c["scalarReferenceDeltaK"]==24
    assert sha(base/"case.json")==done["caseSHA256"] and max(summary["residualWindow"].values())<-8
    out.mkdir(parents=True,exist_ok=False)
    for name in ("mesh.su2","body.json","nodes.npz","case.json"):
        os.link(base/name,out/name)
    # 论文式（6）至（8），热导率单位由 kcal/(m·h·°C) 换为 W/(m·K)，空气物性沿用控制工况。
    d=c["wireDiameterM"];length=c["exposedStraightLengthM"];wire_k=25.2*4184/3600
    re=c["rho"]*c["speedMps"]*d/c["mu"];nu=.841*re**.36;h=nu*c["k"]/d
    mu=np.sqrt(4*h/(wire_k*d));cross=np.pi*d**2/4;delta=c["wallK"]-c["airK"]
    root_w=wire_k*cross*mu*delta
    with np.load(out/"nodes.npz") as z:
        xyz=z["xyz"];wet=np.unique(z["wallMap"])
    vertices=np.array(body["V"]);faces=np.array(body["F"])-1;group=np.array(body["regionID"])
    roots=[min(abs(vertices[faces[group==i],0]).ravel()) for i in (2,6)]
    assert abs(roots[0]-roots[1])<1e-15
    x=xyz[wet,0]*c["diameterM"];axis_distance=np.maximum(abs(x)-roots[0],0.)
    assert abs(max(axis_distance)-length)<1e-15
    theta=delta/24*np.exp(-mu*axis_distance)
    # 表中同时保存弱焓的原生温度值和强标量值。壁反力读数释放强行时仅将通量置零。
    for name,temp,scalar in (("uniform",np.full(len(wet),1.5),np.full(len(wet),.5)),
                             ("profile",1+theta,theta),("release",1+theta,np.zeros(len(wet)))):
        np.savetxt(out/(name+".dat"),np.column_stack((wet,temp,scalar)),fmt=["%d","%.17g","%.17g"])
    state=np.loadtxt(base/"restart.csv",delimiter=",",skiprows=1)
    np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_allclose(state[:,1:4],xyz,atol=1e-12,rtol=0)
    assert np.all(state[wet,9]==.5);state[wet,9]=theta
    np.savetxt(out/"initial_flow.csv",state,delimiter=",",header=HEADER,comments="",fmt=["%d"]+["%.17g"]*9)
    spec=dict(base=str(base),baseFrame=done["frameSHA256"],baseWallSHA256=done["wallSHA256"],meshSHA256=key["meshSHA256"],
              source="Tsubouchi and Sato (1959), pp. 812-813, Eqs. 6-8, s infinite, long-wire limit",
              wireK=wire_k,reynoldsD=re,nusseltD=nu,hWm2K=h,muPerM=mu,rootXM=roots[0],lengthM=length,
              wallNodes=len(wet),temperatureUnitK=24,endK=c["airK"]+delta*np.exp(-mu*length),
              oneLeadRootW=root_w,oneLeadThroughCutW=root_w*np.exp(-mu*length),oneLeadSideW=root_w*(-np.expm1(-mu*length)),
              inputs={n:sha(out/n) for n in ("initial_flow.csv","uniform.dat","profile.dat","release.dat","mesh.su2","body.json","nodes.npz","case.json")})
    write_json(out/"profile.json",spec)


def native(prep, build, dest, start, steps, mode, trace=False):
    """在独立目录按时间、进程族内存和磁盘限额运行。接续仅使用正常结束且输入输出一致的运行。"""
    assert mode in ("uniform","profile","release") and 1<=steps<=320
    exe=build/("wire_observer.exe" if trace else "wire.exe")
    assert read(build/"build.json")["programs"][exe.stem]==sha(exe)
    spec=read(prep/"profile.json");cfg=config(steps)+"MARKER_PYTHON_CUSTOM= (wall)\n"
    if mode=="release":
        assert steps==1 and not trace
        cfg=cfg.replace("(wall, VALUE, 0.5)","(wall, FLUX, 0.0)").replace("SOLUTION, PRIMITIVE","SOLUTION, RESIDUAL")+"WRT_RESTART_COMPACT= NO\n"
    created=not dest.exists()
    if created:
        dest.mkdir()
        for src,name in ((prep/"mesh.su2","mesh.su2"),(prep/(mode+".dat"),"wire_boundary.dat"),(start,"initial_flow.csv")):
            os.link(src,dest/name)
        (dest/"case.cfg").write_text(cfg,encoding="ascii")
    assert (dest/"case.cfg").read_text("ascii")==cfg and sha(dest/"mesh.su2")==spec["meshSHA256"]
    assert sha(dest/"wire_boundary.dat")==spec["inputs"][mode+".dat"] and sha(dest/"initial_flow.csv")==sha(start)
    locks={str(p):sha(p) for p in [exe,build/"build.json",prep/"profile.json",Path(__file__).resolve(),
           *(dest/n for n in ("mesh.su2","wire_boundary.dat","initial_flow.csv","case.cfg"))]}
    threads=4 if steps>20 else 1
    command=[str(exe),"case.cfg","--threads",str(threads)]
    if created:
        # 独立命令行进程内的环境仅作用于本次调用。
        if trace:os.environ["RAH_CM011_LOW_SPEED_LEDGER"]="1"
        else:os.environ.pop("RAH_CM011_LOW_SPEED_LEDGER",None)
        kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build"))
        bin_path=str(kit/"toolchain/w64devkit/bin")
        if bin_path not in os.environ["PATH"].split(os.pathsep):os.environ["PATH"]=bin_path+os.pathsep+os.environ["PATH"]
        observe(command,dest,locks,720 if threads==4 else (180 if steps==20 or mode=="release" else 150),
                (1536 if trace or threads==4 else 850)*1024**2)
    r=read(dest/"debugger_process.json")
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and r["command"]==command and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"]) and r["inputSHA256"]==locks
    assert f'Applied thermistor wall nodes: {spec["wallNodes"]}' in (dest/"debugger.log").read_text("ascii",errors="replace")
    with (dest/"history.csv").open(encoding="ascii") as stream:
        header=[k.strip() for k in next(csv.reader(stream,skipinitialspace=True))]
        history=np.loadtxt(stream,delimiter=",",ndmin=2)
    assert np.isfinite(history).all() and 1<=len(history)<=steps
    np.testing.assert_array_equal(history[:,header.index("Inner_Iter")],np.arange(len(history)))
    window={k:float(history[-20:,header.index(k)].max()) for k in RESIDUALS}
    assert len(history)==steps or (len(history)>=20 and max(window.values())<-8)
    last=len(history)-1
    files=["history.csv","restart.csv","surface.vtk"]
    if trace:files.extend(k+".csv" for k in ("start","beforeFar","afterFar","beforeWall","afterWall","end"))
    with np.load(prep/"nodes.npz") as z:xyz=z["xyz"];wet=np.unique(z["wallMap"])
    targets=np.loadtxt(prep/(mode+".dat"))
    np.testing.assert_array_equal(targets[:,0],wet)
    for i in sorted(set(range(100,last+1,100))|{last}):
        name=f"restart_{i:06d}.csv";files.extend((name,f"surface_{i:06d}.vtk"))
        if mode!="release":
            state=np.loadtxt(dest/name,delimiter=",",skiprows=1)
            assert state.shape==(len(xyz),10) and np.isfinite(state).all()
            np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-12)
            assert np.max(abs(state[wet,5:8]))<1e-10 and np.min(101325+1.18*9*state[:,4])>0
            assert state[:,9].min()>=-1e-12 and state[:,9].max()<=.5+1e-12
            np.testing.assert_allclose(state[wet,9],targets[:,2],rtol=0,atol=2e-14)
    assert sha(dest/name)==sha(dest/"restart.csv")
    result=dict(steps=len(history),last20=window,passed=len(history)>=20 and max(window.values())<-8,
                seconds=r["elapsedSeconds"],familyPeak=r["peakFamilyRSSBytes"],outputs={n:sha(dest/n) for n in files})
    seal=dest/"native.json"
    if seal.exists():assert read(seal)==result
    else:write_json(seal,result)
    return result
