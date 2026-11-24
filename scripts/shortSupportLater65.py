"""接续 Short Test 65 银球 14D 的两条已完成 BDF1 时间链，在 tU/D=1 导出同刻结果。"""
import argparse
import os
from pathlib import Path
import shutil
import sys

for var in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[var]="1"

from shortSupportTime65 import ROOT, PREP, read, sha, write_json
from physical_step import DT, advance, inspect_step
from energy import endpoint
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive

PLAN=ROOT/"docs/validation/v339/PLAN.md"
TARGET=1


def checked_source(job):
    """核对两个出发帧的原生终态、独立释放、时间索引与实际前一层。"""
    job=job.resolve();job.relative_to(ROOT/"outputs")
    receipt=read(job/"completed.json");source=Path(receipt["source"]).resolve();source.relative_to(ROOT/"outputs")
    pair=read(job/"export/pair.json");descriptor=read(source/"job.json")
    assert sha(job/"export/pair.json")==receipt["pairSHA256"] and pair["elapsedDU"]==.25
    assert receipt["fullFrame"]==pair["fullFrame"] and receipt["halfFrame"]==pair["halfFrame"]
    assert pair["faces"]==receipt["faces"]==32480 and pair["meshSHA256"]=="a56a4c3dee8fd8bb74bc23ae4e85b71eaf7d5552dd99fccf8a6c8cb6a40789cd"
    assert pair["seed"]==sha(source/"seed.csv")==descriptor["seed"] and descriptor["baseFrame"]==pair["steadyFrame"]
    prep=PREP;steady=ROOT/"cases/short_sphere/supportedAngularSurface65.zip"
    if "constructionSHA256" in receipt:
        assert sha(job/"construction.json")==receipt["constructionSHA256"]
        c=read(job/"construction.json");prep=(ROOT/c["prepared"]).resolve();prep.relative_to(ROOT/"outputs")
        steady=(ROOT/c["steady"]).resolve().parent/"surface/surface.zip";steady.relative_to(ROOT/"outputs")
        assert c["steadyFrame"]==pair["steadyFrame"] and c["meshSHA256"]==pair["meshSHA256"]
    assert sha(prep/"prepared.json")==descriptor["preparedSHA256"] and pair["meshSHA256"]==sha(prep/"mesh.su2")
    for branch,step,field in (("full",1,"fullEnergy"),("half",1,"firstHalfEnergy"),("half",2,"halfEnergy")):
        stage=source/branch/f"step{step}";e=read(stage/"energy.json")
        assert e==pair[field] and max(e["window"])< -8 and e["rawMax"]<1e-8 and e["heatChange"]<.002 and abs(e["relativeBalance"])<.005
        previous=source/"seed.csv" if step==1 else source/"half/step1/restart_00001.csv"
        frame=inspect_step(source,branch,step,previous,False,prepared=prep)
        inspect_step(source,branch,step,frame,True,prepared=prep)
        assert sha(frame)==e["frame"] and e["old"]==sha(previous) and e["dt"]==DT[branch] and e["step"]==step
        assert all(sha(stage/name)==digest for name,digest in e["files"].items())
    return source,prep,steady,pair


def stage_range():
    for branch,dt in DT.items():
        for step in range(int(.25/dt)+1,int(TARGET/dt)+1):
            yield branch,step


def progress(job,code,message):
    atomic_write_json(job/"progress.json",dict(stage=code,message=message))
    print(message,flush=True)


def prepare(job,source_job):
    source,prep,steady,pair=checked_source(source_job)
    work=job/"run";work.mkdir()
    for name in ("job.json","seed.csv"):os.link(source/name,work/name)
    for branch,step in (("full",1),("half",2)):
        stage=source/branch/f"step{step}";dest=work/branch/stage.name;dest.mkdir(parents=True)
        for name in ("native.json",f"restart_{step:05d}.csv"):os.link(stage/name,dest/name)
    record=dict(schema="short65-time-continuation-v1",targetDU=TARGET,sourceJob=str(source_job.resolve()),
                sourceReceipt=sha(source_job/"completed.json"),sourcePair=sha(source_job/"export/pair.json"),sourceRun=str(source),
                prepared=str(prep),preparedSHA256=sha(prep/"prepared.json"),steadySurface=str(steady),steadySurfaceSHA256=sha(steady),
                earlyFull=pair["fullFrame"],earlyHalf=pair["halfFrame"],seed=pair["seed"],baseFrame=pair["steadyFrame"])
    write_json(work/"continuation.json",record)
    return work


def advance_chain(job,work):
    c=read(work/"continuation.json");source_job=Path(c["sourceJob"]);prep=Path(c["prepared"])
    assert c["sourceReceipt"]==sha(source_job/"completed.json") and c["sourcePair"]==sha(source_job/"export/pair.json")
    assert sha(prep/"prepared.json")==c["preparedSHA256"] and sha(Path(c["steadySurface"]))==c["steadySurfaceSHA256"]
    journal=job/"acceptedSteps.json"
    saved=read(journal)["steps"] if journal.exists() else []
    stages=list(stage_range())
    if [(row["branch"],row["step"]) for row in saved]!=stages[:len(saved)]:
        raise ValueError("已有接续回执的步次与计划不符。")
    ledger=[]
    for index,(branch,step) in enumerate(stages):
        stem=work/branch;stage=stem/f"step{step}";prev=step-1
        previous=stem/f"step{prev}"/f"restart_{prev:05d}.csv"
        if index<len(saved):
            for file in ("native.json","release/native.json","energy.json"):
                if not (stage/file).is_file():
                    raise ValueError("已登记的物理步缺少原始回执："+str(stage/file))
        if step==int(.25/DT[branch])+1:
            assert sha(previous)==c["early"+branch.title()]
        progress(job,"physical",f"{branch} 第 {step} 步：从真实前一时间层推进至 tU/D = {DT[branch]*step:g}")
        if stage.exists():frame=inspect_step(work,branch,step,previous,False,prepared=prep)
        else:frame=advance(work,branch,step,previous,prepared=prep,plan=PLAN)
        if (stage/"release").exists():inspect_step(work,branch,step,frame,True,prepared=prep)
        else:advance(work,branch,step,frame,True,prepared=prep,plan=PLAN)
        e=endpoint(work,branch,step,prepared=prep)
        assert e["step"]==step and e["elapsedDU"]==step*DT[branch] and e["old"]==sha(previous) and e["frame"]==sha(frame)
        ledger.append(dict(branch=branch,step=step,energySHA256=sha(stage/"energy.json"),frame=e["frame"],old=e["old"],
                           elapsedS=e["elapsedS"],wallW=e["wallW"],storageW=e["storageW"],relativeBalance=e["relativeBalance"]))
        # 恢复期间逐步复核原回执。完整的旧记录保留至确认下一步可以追加。
        if index<len(saved):
            if ledger[-1]!=saved[index]:
                raise ValueError("已接受物理步的来源或能量账发生变化："+str(stage))
        else:
            atomic_write_json(journal,dict(steps=ledger))
    return ledger


def launch(job,source_job=None,resume=False):
    job=job.resolve();job.relative_to(ROOT/"outputs");job.mkdir(parents=True,exist_ok=resume)
    with exclusive(job):
        assert not (job/"completed.json").exists()
        if resume:
            work=job/"run"
            if source_job:assert source_job.resolve()==Path(read(work/"continuation.json")["sourceJob"])
        else:
            if shutil.disk_usage(job).free<7*1024**3:raise OSError("接续计算需要 7 GiB 可用磁盘。")
            if source_job is None:
                from shortSupportTime65 import main
                source_job=job/"early"
                progress(job,"early","先生成两种步长的 2.565681 ms 出发帧")
                main(source_job)
            progress(job,"bind","核对两条早期时间链、原生回执、体网及温界")
            work=prepare(job,source_job)
        ledger=advance_chain(job,work)
        progress(job,"export","核对 tU/D = 1 的两种步长、五区热量、近壁测站与周向结果")
        from shortSupportLaterExport import finish
        finish(job,work,ledger)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("output",type=Path)
    parser.add_argument("--from-job",type=Path);parser.add_argument("--resume",action="store_true")
    args=parser.parse_args();launch(args.output,args.from_job,args.resume)
