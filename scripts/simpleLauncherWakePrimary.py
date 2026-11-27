"""固定低速风洞的两半步细热网计算。工作台主结果只读本次作业的终帧。"""
import json
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from runSimpleLauncherWakeCase import CASE, DOC, ROOT, flow_available, flow_sources, matches, momentum_step, read, saved_record
sys.path.insert(0,str(DOC.parent/"v355"))
from refined import launch_step, save
from compare import stored_power
from finalize import read_step
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe
from simpleLauncherTunnel import geometry
import numpy as np


def fine_seed(job,upstream,prepared,allow_cached):
    """只复用有稳态基流、几何、温度初场和网格身份的输入；新工况开始时尚无瞬态热场。"""
    cache=ROOT/"outputs/v355_validation/run01"
    expected=read(DOC.parent/"v355/evidence.json")["initial"]
    proof=saved_record(cache/"initial/results.json") if allow_cached else {}
    summary=saved_record(cache/"steady/summary.json") if allow_cached else {}
    if (proof==expected and proof["accepted"] and all(proof["gates"].values()) and
        summary.get("sourceRestartSHA256")==sha(upstream/"coarse/restart.csv")==proof["sourceFlowSHA256"] and
        summary.get("surfaceSHA256")==read(CASE/"geometry.json")["fileSHA256"] and summary.get("heatSpacingD")==.002 and
        matches(cache/"steady/temperature.npz",proof["sourceHeatMeshSHA256"]) and
        matches(cache/"initial/fields.npz",proof["fieldsSHA256"])):
        return cache,True
    # 初场缺失时在独立进程重建，超时或中断保留原始日志，并在进入第一个热步前停止。
    folder=job/"seed_preparation"
    folder.mkdir()
    script=DOC.parent/"v355/refined.py"
    locks=[script,CASE/"tunnel.json",CASE/"geometry.json",upstream/"coarse/restart.csv",
           upstream/"coarse/condition.json",upstream/"heat_coarse/temperature.npz"]
    command=[sys.executable,str(script),str(upstream),str(job/"flow"),str(prepared),str(job),"prepare"]
    rec=observe(command,folder,{str(p):sha(p) for p in locks},180,300*1024**2)
    if rec["debuggerExitCode"] or rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]:
        raise RuntimeError(f"细热网初场重建未完成，详见 {folder/'debugger.log'}")
    proof=read(job/"initial/results.json")
    assert proof["accepted"] and all(proof["gates"].values()) and proof["sourceFlowSHA256"]==sha(upstream/"coarse/restart.csv")
    assert matches(job/"steady/temperature.npz",proof["sourceHeatMeshSHA256"]) and matches(job/"initial/fields.npz",proof["fieldsSHA256"])
    return job,False


def solve_fine(job,allow_cached=True):
    job.relative_to(ROOT/"outputs")
    job.mkdir(parents=True,exist_ok=False)
    upstream,prepared=flow_sources(job,allow_cached)
    seed,seed_reused=fine_seed(job,upstream,prepared,allow_cached)
    condition=read(CASE/"wake/condition.json")
    boundary=read(CASE/"tunnel.json")
    assert all(condition[k]==v for k,v in boundary.items() if k not in ("scope","caseLabel"))
    definition,tri,diameter,_=geometry()
    assert definition["fileSHA256"]==condition["surfaceSHA256"]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    index=read(CASE/"wake/flow_index.json")
    names=["disturbed_half_1","disturbed_half_2"]
    flow=(ROOT/"outputs/v354_validation/flow01").resolve()
    flow_reused=allow_cached and flow_available(flow,{name:index[name] for name in names},prepared)
    if not flow_reused:
        flow=job/"flow"
    initial=seed/"initial/fields.npz"
    mesh_file=seed/"steady/temperature.npz"
    proof=read(seed/"initial/results.json")
    with np.load(mesh_file) as source:
        mesh={k:source[k] for k in ("xy","triangles")}
    scale=(2*np.pi*diameter**3*condition["densityKgM3"]*condition["specificHeatJkgK"]*
           (condition["wallTemperatureK"]-condition["ambientTemperatureK"]))
    steps=[]
    dt=condition["physicalStepSizeS"]
    for number,name in enumerate(names,1):
        if not momentum_step(flow,"disturbed","half",number,prepared):
            raise RuntimeError(f"{name} 动量终帧未通过")
        launch_step(upstream,flow,prepared,job,name,seed_directory=seed,registered_frame=flow_reused,cached_modes=5)
        table,record=read_step(flow,job/"heat","disturbed","half",number,area,definition)
        f,h=record["flow"],record["heat"]
        assert abs(f["dtDU"]*diameter/.05-dt)<1e-13 and abs(f["elapsedPhysicalS"]-number*dt)<1e-13
        assert h["heatMeshSHA256"]==proof["sourceHeatMeshSHA256"] and h["initialThetaSHA256"]==sha(initial)
        if steps:
            assert f["seedSHA256"]==steps[-1]["flow"]["frameSHA256"]
        final=job/"heat"/name/"thermal.npz"
        storage=stored_power(initial,final,mesh,dt,scale)
        assert abs(storage-h["storageW"])<1e-14 and abs(h["outwardW"]-h["farW"]-storage)<1e-12
        assert record["heatPeakRSSBytes"]<4*1024**3 and record["heatSeconds"]<240
        record["independentStorageW"]=storage
        steps.append(record)
        initial=final
    assert h["elapsedPhysicalS"]==condition["elapsedPhysicalTimeS"]
    assert np.all((table[:,4]>=condition["ambientTemperatureK"]-1e-12) &
                  (table[:,4]<=condition["wallTemperatureK"]+1e-12))
    summary=dict(accepted=True,primaryHeatGrid="fine",condition=condition,surfaceSHA256=definition["fileSHA256"],
                 originalFaces=len(area),areaM2=float(area.sum()),wallSHA256=h["wallSHA256"],
                 outwardW=h["outwardW"],regionW=h["regionW"],sourceFrameSHA256=f["frameSHA256"],
                 thermalFrameSHA256=h["fieldSHA256"],initialFlowSHA256=steps[0]["flow"]["seedSHA256"],
                 initialThermalSHA256=proof["fieldsSHA256"],flowMeshSHA256=f["meshSHA256"],
                 heatMeshSHA256=h["heatMeshSHA256"],storageW=h["storageW"],farW=h["farW"],
                 pressureRangePa=[table[:,2].min(),table[:,2].max()],
                 heatFluxRangeWm2=[table[:,3].min(),table[:,3].max()],verifiedByExperimentM2=0)
    dest=job/"result"
    dest.mkdir()
    bundle=dest/"wind_tunnel.zip"
    with ZipFile(bundle,"x",ZIP_DEFLATED) as z:
        z.writestr("summary.json",json.dumps(summary,ensure_ascii=False,allow_nan=False,indent=2))
        z.write(job/"heat"/names[-1]/"wall.csv","wall.csv")
    save(dest/"receipt.json",dict(sha256=sha(bundle),faces=len(area),outwardW=h["outwardW"]))
    save(job/"completed.json",dict(flowReusedAfterVerification=flow_reused,initialReusedAfterVerification=seed_reused,
         flowDirectory=str(flow),seedDirectory=str(seed),steps=steps,summary=summary))
    print("细热网两半步完成",h["wallSHA256"],flush=True)
