"""默认五步工作台的 0.05 m/s 风洞尾迹快照：内部准备初场、推进流动、复算显热。"""
import argparse
import json
import os
from pathlib import Path
import sys
import uuid

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[variable]="1"

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs/validation/v354"
sys.path[:0]=[str(DOC),str(DOC.parent/"v353"),str(DOC.parent/"v352")]
from finalize import finalize
from prepare_wake import prepare as prepare_full
from run_control import run as check_initial_heat
from run_flow import run as momentum_step
from runSU2NativeDebugger import observe
from runEmbeddedSU2 import sha
from simpleLauncherTemperature import evaluate
from simpleLauncherTunnel import CASE,prepare
from thermistorBead import solve


def read(p):
    return json.loads(p.read_bytes())


def matches(path,digest):
    return path.is_file() and sha(path)==digest


def saved_record(path):
    """损坏或不完整的缓存记录按未命中处理，原目录留存。"""
    try:
        record=read(path)
        return record if isinstance(record,dict) else {}
    except (OSError,ValueError):
        return {}


def sources(job,allow_cached=True):
    """只复用指纹已知的本机基流与初场。缺少时由冻结 STL 和风洞边界重新生成。"""
    ref=read(DOC.parent/"v352/baseline.json")["coarse"]
    old=read(DOC.parent/"v350/evidence.json")["fresh"]["runDirectory"]
    upstream=ROOT/"outputs"/Path(old.replace("\\","/")).name
    receipt=saved_record(upstream/"coarse/engine_process.json")
    condition_file=upstream/"coarse/condition.json"
    expected_wall=read(DOC.parent/"v353/evidence.json")["coarse"]["sourceFaceSHA256"]
    if not (allow_cached and receipt.get("status")=="terminal" and receipt.get("exitCode")==0 and
            not (receipt.get("cancelled") or receipt.get("timedOut")) and
            matches(condition_file,receipt.get("inputSHA256",{}).get(str(condition_file))) and
            matches(upstream/"coarse/restart.csv",ref["flowSHA256"]) and
            matches(upstream/"heat_coarse/temperature.npz",ref["heatMeshSHA256"]) and
            matches(upstream/"heat_coarse/wall.csv",expected_wall)):
        upstream=job/"axisymmetric"
        flow=upstream/"coarse";prepare(flow,"coarse")
        request=read(flow/"engine_request.json");request["jobToken"]=str(uuid.uuid4())
        (flow/"engine_request.json").write_text(json.dumps(request),encoding="utf8")
        solve(flow,memory_bytes=4*1024**3,output_bytes=2*1024**3)
        evaluate(flow,upstream/"heat_coarse",.004)
    assert matches(upstream/"coarse/restart.csv",ref["flowSHA256"])
    assert matches(upstream/"heat_coarse/temperature.npz",ref["heatMeshSHA256"])
    prepared=ROOT/"outputs/v352_validation/prepared01"
    expected=read(DOC.parent/"v352/evidence.json")["geometry"]["files"]
    if not (allow_cached and saved_record(prepared/"prepared.json").get("files")==expected and
            all(matches(prepared/name,digest) for name,digest in expected.items())):
        prepared=job/"prepared"
        prepare_full(prepared,upstream)
    initial=ROOT/"outputs/v353_validation/coarse01/fields.npz"
    seed_record=saved_record(initial.parent/"results.json")
    if not (allow_cached and seed_record.get("accepted") and seed_record.get("sourceFlowSHA256")==ref["flowSHA256"] and
            matches(initial,read(DOC.parent/"v353/evidence.json")["coarse"]["fieldsSHA256"])):
        initial=job/"initial_heat/fields.npz"
        assert check_initial_heat(upstream,"coarse",initial.parent)
    return upstream,prepared,initial


def flow_available(folder,index,prepared):
    """检查六个实际流动终帧、初始文件链、残差历史与本模型网格。"""
    for name,expected in index.items():
        f=folder/name
        if not all((f/k).is_file() for k in ("step.json","debugger_process.json","debugger.log","case.cfg")):
            return False
        rec=saved_record(f/"step.json");n=rec.get("step")
        run=saved_record(f/"debugger_process.json")
        if not (rec.get("accepted") and isinstance(n,int) and run.get("status")=="terminal" and
                run.get("debuggerExitCode")==0 and not (run.get("stopReason") or run.get("survivingProcessIDs")) and
                all(rec.get(k)==v for k,v in expected.items())):return False
        checks={f"restart_{n:05d}.csv":"frameSHA256",f"seed_{n-1:05d}.csv":"seedSHA256",
                f"history_{n:05d}.csv":"historySHA256","mesh.su2":"meshSHA256"}
        if not all(matches(f/key,rec[identity]) for key,identity in checks.items()):return False
        if rec["meshSHA256"]!=read(prepared/"prepared.json")["files"]["mesh.su2"]:return False
    return len(index)==6


def fresh(job,allow_cached=True):
    job.relative_to(ROOT/"outputs");job.mkdir(parents=True,exist_ok=False)
    upstream,prepared,initial=sources(job,allow_cached)
    index=read(CASE/"wake/flow_index.json")
    flow=ROOT/"outputs/v354_validation/flow01"
    reused=allow_cached and flow_available(flow,index,prepared)
    if not reused:flow=job/"flow"
    heat=job/"heat"
    for size,count in (("full",1),("half",2)):
        for branch in ("control","disturbed"):
            for step in range(1,count+1):
                name=f"{branch}_{size}_{step}"
                if not momentum_step(flow,branch,size,step,prepared):
                    raise ValueError(f"{name} 的动量时间步未通过")
                frame=flow/name/f"restart_{step:05d}.csv";out=heat/name;out.mkdir(parents=True)
                script=DOC/"heat_step.py"
                paths=[script,DOC/"flux_transfer.py",DOC/"PLAN.md",frame,flow/name/"step.json",
                       initial,upstream/"coarse/restart.csv",upstream/"heat_coarse/temperature.npz",
                       prepared/"stations.npz",prepared/"control.csv",ROOT/"scripts/simpleLauncherAngularHeat.py"]
                if step==2:paths.append(heat/f"{branch}_half_1/thermal.npz")
                rec=observe([sys.executable,str(script),str(upstream),str(flow),str(prepared),str(heat),name,
                             "--initial-field",str(initial)],out,{str(p):sha(p) for p in paths},240,300*1024**2)
                if rec["debuggerExitCode"] or rec["stopReason"] or rec["survivingProcessIDs"]:
                    raise RuntimeError(f"{name} 的三维显热未通过，详见 {out/'debugger.log'}")
    summary=finalize(flow,heat,job/"result")
    record=dict(flowReusedAfterVerification=reused,flowDirectory=str(flow),initialField=str(initial),
                surfaceSHA256=summary["surfaceSHA256"],wallSHA256=summary["wallSHA256"],
                thermalFrameSHA256=summary["thermalFrameSHA256"])
    (job/"completed.json").write_text(json.dumps(record,indent=2),encoding="utf8")
    print(summary["wallSHA256"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("job",type=Path)
    parser.add_argument("--from-start",action="store_true",help="从模型与风洞边界重新生成初场并推进全部动量步")
    args=parser.parse_args();fresh(args.job.resolve(),not args.from_start)
