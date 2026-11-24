"""逐层核对接续帧、储热及释放反力，再从真实体场重建两套终点和同面时间差。"""
import argparse
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import sys
from zipfile import ZIP_DEFLATED, ZipFile

for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[var]="1"
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CASE=ROOT/"cases/short_sphere"
sys.path.insert(0,str(HERE.parent/"v333"))
from verify_time import verify,table
from build_time import read,sha

def validation_module(name,path):
    # 各版本有同名 pack.py，按登记版本的完整路径载入独立核对程序。
    spec=importlib.util.spec_from_file_location(name,HERE.parent/path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

check_ledger=validation_module("short65_time_ledger","v333/pack.py").check
circumference=validation_module("short65_circumference","v338/verify.py")
STEPS=[("full",i) for i in range(2,5)]+[("half",i) for i in range(3,9)]


def copy_archive(job,dest):
    """只发布可直接从工作台读取的两份表面解、小型诊断和源摘要。"""
    files={"supportedTime65.json":"pair.json","supportedTime65Full.zip":"full/surface.zip","supportedTime65Half.zip":"half/surface.zip",
           "stationGradients65.json":"stationGradients.json","stationAzimuth.zip":"stationAzimuth.zip",
           "stationAzimuthFull.zip":"stationAzimuthFull.zip","timeDifference.zip":"timeDifference.zip"}
    receipt=read(job/"completed.json");assert sha(job/"export/pair.json")==receipt["pairSHA256"]
    dest.mkdir(parents=True,exist_ok=True)
    for name,rel in files.items():
        old=job/"export"/rel;new=dest/name
        if new.exists():assert sha(new)==sha(old)
        else:shutil.copyfile(old,new)


def compare_reread(original,repeated):
    """重读只用同一源任务，两个时刻的回执和科学小包均应逐字节一致。"""
    files=["completed.json","export/pair.json","export/full/surface.zip","export/half/surface.zip",
           "export/stationGradients.json","export/stationAzimuth.zip"]
    receipt=read(original/"completed.json")
    if "continuationSHA256" in receipt:
        files.extend(["continuation.json","export/stationAzimuthFull.zip","export/timeDifference.zip"])
    if "constructionSHA256" in receipt:
        files.append("construction.json")
    hashes={}
    for name in files:
        hashes[name]=sha(original/name)
        assert hashes[name]==sha(repeated/name),name
    return hashes


def evidence(job,dest):
    work=job/"run";members={};large={}
    with ZipFile(dest,"x",compression=ZIP_DEFLATED,compresslevel=8) as z:
        def add(name,data):
            members[name]=hashlib.sha256(data).hexdigest();z.writestr(name,data)
        # 回执锁定的是 Windows 实际读取的字节，也保存其原换行，便于从 Git 归档核对。
        locked=(HERE/"PLAN.md",HERE.parent/"v333/physical_step.py")
        for path in locked:
            add("locked/"+path.name,path.read_bytes())
        for f in (job/"completed.json",job/"acceptedSteps.json",job/"continuation.json",work/"job.json",*job.glob("mapping_*/export_process.json")):
            add(f.relative_to(job).as_posix(),f.read_bytes())
        for branch,step in STEPS:
            stage=work/branch/f"step{step}";e=read(stage/"energy.json")
            for folder in (stage,stage/"release"):
                process=read(folder/"debugger_process.json")
                for path in locked:
                    assert members["locked/"+path.name]==process["inputSHA256"][str(path)]
                for f in folder.iterdir():
                    if not f.is_file():continue
                    if f.name in ("case.cfg","native.json","observed.json","energy.json","debugger_process.json","wall_budget.npz",
                                  "debugger.log","far_pairs.csv","stages.csv") or f.name.startswith("history_"):
                        add(f.relative_to(job).as_posix(),f.read_bytes())
            for name,h in e["files"].items():
                if name.endswith(".bin"):large[(stage/name).relative_to(job).as_posix()]=h
            for name,key in ((f"seed_{step-1:05d}.csv","old"),(f"restart_{step:05d}.csv","frame"),(f"surface_{step:05d}.vtk","vtk")):
                assert sha(stage/name)==e[key];large[(stage/name).relative_to(job).as_posix()]=e[key]
            old=np.loadtxt(stage/f"seed_{step-1:05d}.csv",delimiter=",",skiprows=1,usecols=9)
            new=np.loadtxt(stage/f"restart_{step:05d}.csv",delimiter=",",skiprows=1,usecols=9)
            nodes=np.fromfile(stage/"nodes.bin","<f8").reshape(-1,18);ix=np.argsort(nodes[:,0]);volume=nodes[ix,5]
            np.testing.assert_array_equal(nodes[ix,0],np.arange(len(new)))
            before=np.fromfile(stage/"before_time.bin","<f8")[ix];after=np.fromfile(stage/"after_time.bin","<f8")[ix]
            begin=np.arange(0,len(new),4096)
            sums=np.column_stack((begin,np.minimum(begin+4096,len(new)),
                    np.add.reduceat(volume*(new-old)*e["scale"]/e["dt"],begin),np.add.reduceat((after-before)*e["scale"],begin)))
            buf=io.StringIO();np.savetxt(buf,sums,delimiter=",",comments="",fmt=["%d","%d","%.17g","%.17g"],
                  header="firstID,endExclusive,storageFromFramesW,storageFromNativeW")
            add(f"storage/{branch}{step}.csv",buf.getvalue().encode())
        z.writestr("manifest.json",json.dumps(dict(files=members,largeLocal=large),indent=2))
    check_ledger(dest,prefix="run",steps=STEPS)


def independent(job,archive):
    work=job/"run";c=read(work/"continuation.json");prep=Path(c["prepared"]);early=Path(c["sourceJob"])
    steady=ROOT/"outputs/v332_validation/run01/steady05/restart_000199.csv"
    if (early/"construction.json").is_file():
        path=ROOT/read(early/"construction.json")["steady"];steady=path/read(path/"converged.json")["stage"]/"restart_000199.csv"
    pair=read(archive/"supportedTime65.json");previous=read(early/"export/pair.json")
    proof=dict(raw=verify(archive,work,base_frame=c["baseFrame"],steady=steady))
    with ZipFile(archive/"timeDifference.zip") as diffs:
        for branch,step in (("full",4),("half",8)):
            data=pair["history"][branch];dt=1/step
            prior=pair["seed"]
            for k,r in enumerate(data,1):
                assert r["old"]==prior and r["step"]==k and r["dt"]==dt and r["elapsedDU"]==k*dt
                folder=(Path(c["sourceRun"]) if k*dt<=.25 else work)/branch/f"step{k}"
                assert r["energySHA256"]==sha(folder/"energy.json") and r["frame"]==sha(folder/f"restart_{k:05d}.csv")
                prior=r["frame"]
            assert prior==pair[branch+"Frame"] and data[round(.25/dt)-1]["frame"]==previous[branch+"Frame"]
            with ZipFile(archive/("supportedTime65"+branch.title()+".zip")) as z, ZipFile(early/"export"/branch/"surface.zip") as q:
                w=table(z,"wall.csv");old=table(q,"wall.csv");s=json.loads(z.read("summary.json"));first=json.loads(q.read("summary.json"))
            d=table(diffs,branch+".csv");np.testing.assert_array_equal(d["referenceFaceID"],old["faceID"])
            np.testing.assert_array_equal(w["faceID"],old["faceID"]);np.testing.assert_array_equal(w["wallInputK"],old["wallInputK"])
            for col,k in (("deltaGaugePa","gaugePa"),("deltaFluxInWm2","fluxInWm2"),("deltaFirstK","firstK")):
                np.testing.assert_array_equal(d[col],w[k]-old[k]);stats=pair["timeComparison"][branch]["fields"][k]
                np.testing.assert_allclose(np.sqrt(np.average(d[col]**2,weights=w["areaM2"])),stats["rms"],rtol=1e-12)
            region=[-w["areaM2"][w["regionID"]==k]@d["deltaFluxInWm2"][w["regionID"]==k] for k in range(1,6)]
            np.testing.assert_allclose(region,pair["timeComparison"][branch]["regionChangeW"],atol=1e-12,rtol=0)
            assert abs(sum(region)-(s["energy"]["wallW"]-first["energy"]["wallW"]))<1e-10
            bundle=archive/("stationAzimuth.zip" if branch=="half" else "stationAzimuthFull.zip")
            stage=work/branch/f"step{step}"
            circ=circumference.verify(bundle,True,pair=pair,stage=stage,prepared=prep,branch=branch)
            a=np.loadtxt(stage/f"restart_{step:05d}.csv",delimiter=",",skiprows=1,usecols=(1,2,3,5,6,7))
            n=round(.25/dt);base=Path(c["sourceRun"])/branch/f"step{n}"/f"restart_{n:05d}.csv"
            old_u=np.loadtxt(base,delimiter=",",skiprows=1,usecols=(5,6,7))
            wake=((a[:,0]+1.15)/.55)**2+(a[:,1]**2+a[:,2]**2)/.45**2<1
            mod=np.linalg.norm(a[wake,3:]-old_u[wake],axis=1)
            nodes=np.fromfile(stage/"nodes.bin","<f8").reshape(-1,18);weights=nodes[np.argsort(nodes[:,0]),5][wake]
            fields={}
            for key in ("absolutePa","gaugePa","fluxInWm2","firstK","wallInputK"):
                values=w[key]
                assert np.isfinite(values).all()
                low=int(np.argmin(values));high=int(np.argmax(values))
                fields[key]=dict(minimum=values[low],minimumFace=int(w["faceID"][low]),
                                 maximum=values[high],maximumFace=int(w["faceID"][high]))
            proof[branch]=dict(azimuth=circ,regionalDeltaW=region,surfaceRanges=fields,
                               maxWakeVelocityChangeOverU=max(mod),wakeRMSChangeOverU=np.sqrt(np.average(mod**2,weights=weights)))
    return proof


def verify_release(archive,receipts):
    """只用版本内的小文件核对物理链、计算读数及发布包，让无体场的安装也能检查。"""
    pair=read(archive/"supportedTime65.json")
    early=read(CASE/"supportedTime65.json")
    assert pair["schema"]=="short65-supported-time-continuation-v1" and pair["seed"]==early["seed"]
    assert pair["elapsedS"]==4*early["elapsedS"] and pair["elapsedDU"]==1
    check_ledger(receipts,prefix="run",steps=STEPS)
    verify(archive)
    assert sha(archive/"stationGradients65.json")==pair["stationGradientsSHA256"]
    changes=archive/"timeDifference.zip"
    assert sha(changes)==pair["timeComparison"]["differenceSHA256"]
    with ZipFile(receipts) as proof, ZipFile(changes) as differences:
        completed=json.loads(proof.read("completed.json"))
        assert completed["pairSHA256"]==sha(archive/"supportedTime65.json")
        assert completed["continuationSHA256"]==hashlib.sha256(proof.read("continuation.json")).hexdigest()
        steps=json.loads(proof.read("acceptedSteps.json"))["steps"]
        assert [(s["branch"],s["step"]) for s in steps]==STEPS
        for branch in ("full","half"):
            dt=pair[branch+"Energy"]["dt"]
            chain=pair["history"][branch]
            assert len(chain)*dt==1
            previous=pair["seed"]
            for k,row in enumerate(chain,1):
                assert row["old"]==previous and row["step"]==k and row["dt"]==dt and row["elapsedDU"]==k*dt
                assert row["elapsedS"]==k*dt*pair["condition"]["diameterM"]/pair["condition"]["speedMps"]
                if k*dt>.25:
                    record=proof.read(f"run/{branch}/step{k}/energy.json")
                    e=json.loads(record)
                    assert row["energySHA256"]==hashlib.sha256(record).hexdigest()
                    assert all(row[key]==e[key] for key in ("step","old","frame","wallW","storageW","elapsedS","relativeBalance"))
                previous=row["frame"]
            assert previous==pair[branch+"Frame"]==completed[branch+"Frame"]
            assert chain[round(.25/dt)-1]["frame"]==early[branch+"Frame"]
            filename="supportedTime65"+branch.title()+".zip"
            with ZipFile(archive/filename) as last, ZipFile(CASE/filename) as first:
                w=table(last,"wall.csv")
                old=table(first,"wall.csv")
            diff=table(differences,branch+".csv")
            np.testing.assert_array_equal(diff["referenceFaceID"],w["faceID"])
            np.testing.assert_array_equal(old["faceID"],w["faceID"])
            for key,column in (("gaugePa","deltaGaugePa"),("fluxInWm2","deltaFluxInWm2"),("firstK","deltaFirstK")):
                np.testing.assert_array_equal(diff[column],w[key]-old[key])
            suffix="Full" if branch=="full" else ""
            curve=archive/("stationAzimuth"+suffix+".zip")
            assert sha(curve)==pair["stationAzimuth"+suffix+"SHA256"]
            circumference.verify(curve,False,pair=pair,branch=branch)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--job",type=Path);p.add_argument("--publish",action="store_true")
    p.add_argument("--pack",action="store_true");p.add_argument("--raw",action="store_true")
    p.add_argument("--compare-reread",type=Path);a=p.parse_args()
    archive=CASE/"later1";dest=HERE/"evidence";dest.mkdir(exist_ok=True)
    if a.publish:copy_archive(a.job.resolve(),archive)
    if a.pack:evidence(a.job.resolve(),dest/"receipts.zip")
    if a.raw:
        result=independent(a.job.resolve(),archive)
        with (dest/"comparison.json").open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2);f.write("\n")
    if a.compare_reread:
        print(json.dumps(compare_reread(a.job,a.compare_reread),indent=2))
    if not(a.publish or a.pack or a.raw or a.compare_reread):
        verify_release(archive,dest/"receipts.zip")
        print("Verified both time chains, nine energy ledgers, surface fields, temporal differences and azimuth curves.")
