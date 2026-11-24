"""审核银球 10.263 ms 的完整物理时间链，生成全表面场、两时刻同面差及来源回执。"""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

for var in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[var]="1"
import numpy as np
from shortSupportTime65 import ROOT, read, sha
from physical_step import DT, inspect_step
from energy import endpoint
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive
from time_surface import export
from shortSupportStationGradient import save_pair
from shortSupportAzimuth import compute as azimuth


def history(work,c,audited,recheck):
    source=Path(c["sourceRun"]);manifest=read(work/"job.json")
    assert sha(work/"job.json")==sha(source/"job.json") and sha(work/"seed.csv")==c["seed"]
    assert len(audited)==9
    chains={}
    for branch,dt in DT.items():
        chain=[];prior=source/"seed.csv"
        for step in range(1,int(1/dt)+1):
            late=step>int(.25/dt);run=work if late else source
            folder=run/branch/f"step{step}";file=folder/"energy.json";e=read(file)
            frame=folder/f"restart_{step:05d}.csv"
            assert e["step"]==step and e["dt"]==dt and e["elapsedDU"]==dt*step and e["old"]==sha(prior)
            assert e["elapsedS"]==dt*step*manifest["condition"]["diameterM"]/manifest["condition"]["speedMps"]
            assert e["frame"]==sha(frame) and e["old"]==sha(folder/f"seed_{step-1:05d}.csv")
            assert e["rawMax"]<1e-8 and max(e["window"])< -8 and e["heatChange"]<.002 and abs(e["relativeBalance"])<.005
            for place in (folder,folder/"release"):
                proc=read(place/"debugger_process.json")
                assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
                assert not(proc["survivingProcessIDs"] or proc["supervisorError"] or proc["stopReason"])
            if late:
                receipt,=[r for r in audited if r["branch"]==branch and r["step"]==step]
                assert receipt["energySHA256"]==sha(file) and receipt["old"]==e["old"] and receipt["frame"]==e["frame"]
                if recheck:
                    prep=Path(c["prepared"])
                    inspect_step(work,branch,step,prior,False,prepared=prep);inspect_step(work,branch,step,frame,True,prepared=prep)
                    assert endpoint(work,branch,step,prepared=prep)==e
            chain.append(dict(step=step,dt=dt,elapsedDU=dt*step,elapsedS=e["elapsedS"],old=e["old"],frame=e["frame"],
                              energySHA256=sha(file),wallW=e["wallW"],storageW=e["storageW"],relativeBalance=e["relativeBalance"]))
            prior=frame
        assert chain[int(.25/dt)-1]["frame"]==c["early"+branch.title()]
        chains[branch]=chain
    return chains


def read_wall(file):
    with ZipFile(file) as z:
        summary=json.loads(z.read("summary.json"))
        wall=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    return summary,wall


def time_difference(early,out,pair,early_pair):
    comparison=dict(earlyTimeS=early_pair["elapsedS"],laterTimeS=pair["elapsedS"])
    assert pair["elapsedS"]==4*early_pair["elapsedS"] and pair["seed"]==early_pair["seed"]
    with ZipFile(out/"timeDifference.zip","x") as z:
        for branch in DT:
            first=early/"export"/branch/"surface.zip";last=out/branch/"surface.zip"
            assert sha(first)==early_pair[branch+"Archive"] and sha(last)==pair[branch+"Archive"]
            a,w=read_wall(first);b,v=read_wall(last)
            assert a["dtDU"]==b["dtDU"] and a["condition"]==b["condition"] and a["meshSHA256"]==b["meshSHA256"]
            for key in ("faceID","regionID","areaM2","wallInputK"):np.testing.assert_array_equal(w[key],v[key])
            area=w["areaM2"];columns=[];stats={}
            for field in ("gaugePa","fluxInWm2","firstK"):
                d=v[field]-w[field];columns.append(d);i=np.argmax(abs(d))
                stats[field]=dict(rms=np.sqrt(np.average(d*d,weights=area)),maxAbs=abs(d[i]),atFace=int(w["faceID"][i]),signed=d[i])
            power=b["energy"]["wallW"]-a["energy"]["wallW"]
            assert abs(area@columns[1]+power)<1e-10
            groups=[y["outW"]-x["outW"] for x,y in zip(a["regions"],b["regions"])]
            assert abs(sum(groups)-power)<1e-10
            comparison[branch]=dict(earlyFrame=a["frameSHA256"],laterFrame=b["frameSHA256"],wallChangeW=power,
                regionChangeW=groups,earlyRegionW=[x["outW"] for x in a["regions"]],fields=stats,
                earlyPeak=a["peakOutFace"],laterPeak=b["peakOutFace"],earlyPeakWm2=a["peak"]["peakWm2"],laterPeakWm2=b["peak"]["peakWm2"],
                earlyNu=early_pair[branch+"Nu"],laterNu=pair[branch+"Nu"],earlyStationRMSE=a["stationRMSE"],laterStationRMSE=b["stationRMSE"],
                earlyCalorimetryError=a["equivalentRelativeError"],laterCalorimetryError=b["equivalentRelativeError"])
            buf=io.StringIO();np.savetxt(buf,np.column_stack([w["faceID"],*columns]),fmt=["%d","%.17g","%.17g","%.17g"],
                delimiter=",",header="referenceFaceID,deltaGaugePa,deltaFluxInWm2,deltaFirstK",comments="")
            z.writestr(ZipInfo(branch+".csv"),buf.getvalue().encode("ascii"),compress_type=ZIP_DEFLATED,compresslevel=9)
    comparison["differenceSHA256"]=sha(out/"timeDifference.zip")
    return comparison


def publish(job,work,recheck=False):
    work.relative_to(ROOT/"outputs");job.relative_to(ROOT/"outputs")
    c=read(work/"continuation.json");early=Path(c["sourceJob"]);receipt=read(early/"completed.json")
    assert sha(early/"completed.json")==c["sourceReceipt"] and sha(early/"export/pair.json")==c["sourcePair"]==receipt["pairSHA256"]
    original=read(early/"export/pair.json");steps=read(work.parent/"acceptedSteps.json")["steps"]
    chains=history(work,c,steps,recheck);prep=Path(c["prepared"])
    tmp=Path(tempfile.mkdtemp(prefix="export_",dir=str(job)))
    out=tmp/"export"
    pair=export(work,out,prepared=prep,prepared_sha=c["preparedSHA256"],base_frame=c["baseFrame"],
                steady_surface=Path(c["steadySurface"]),continuation=chains)
    gradients=save_pair(prep,work,out);curve=azimuth(prep,work,out);full_curve=azimuth(prep,work,out,branch="full")
    pair["timeComparison"]=time_difference(early,out,pair,original)
    pair.update(stationGradientsSHA256=gradients,stationAzimuthSHA256=curve,stationAzimuthFullSHA256=full_curve)
    for branch,file in (("full","stationAzimuthFull.zip"),("half","stationAzimuth.zip")):
        with ZipFile(out/file) as z:pair["timeComparison"][branch]["rear"]=json.loads(z.read("summary.json"))["rear"]
    with (out/"pair.json").open("w",encoding="utf8",newline="\n") as f:json.dump(pair,f,indent=2,allow_nan=False);f.write("\n")
    ack=dict(source=str(work),pairSHA256=sha(out/"pair.json"),fullFrame=pair["fullFrame"],halfFrame=pair["halfFrame"],
             stationGradientsSHA256=gradients,stationAzimuthSHA256=curve,stationAzimuthFullSHA256=full_curve,
             continuationSHA256=sha(work/"continuation.json"),faces=pair["faces"],areaM2=pair["areaM2"],elapsedS=pair["elapsedS"])
    files={"continuation.json":work/"continuation.json"}
    if "constructionSHA256" in receipt:
        assert sha(early/"construction.json")==receipt["constructionSHA256"]
        files["construction.json"]=early/"construction.json";ack["constructionSHA256"]=receipt["constructionSHA256"]
    for name,file in files.items():
        if (job/name).exists():assert sha(job/name)==sha(file)
        else:shutil.copyfile(file,job/name)
    if (job/"export").exists():
        for name in ("pair.json","full/surface.zip","half/surface.zip","timeDifference.zip","stationGradients.json","stationAzimuth.zip","stationAzimuthFull.zip"):
            assert sha(job/"export"/name)==sha(out/name)
        tmp.resolve().relative_to(job.resolve())
        assert tmp.parent==job and tmp.name.startswith("export_")
        shutil.rmtree(tmp)
    else:
        out.replace(job/"export");tmp.rmdir()
    atomic_write_json(job/"completed.json",ack)
    atomic_write_json(job/"progress.json",dict(stage="complete",message="两条时间链和同一时刻的五区面表已保存"))


def finish(job,work,ledger):
    # 映射与原生求解分别监督，读回不会改写已接受的体场。
    assert len(ledger)==9
    sys.path.insert(0,str(ROOT/"docs/validation/v272"))
    from observe_flux import invoke
    log=Path(tempfile.mkdtemp(prefix="mapping_",dir=str(job)))
    invoke("export",[sys.executable,Path(__file__),job,"--run",work],ROOT,log,240,2,bin_directory=Path(sys.executable).parent)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path);p.add_argument("--run",required=True,type=Path)
    p.add_argument("--reread",action="store_true",help="在新目录从已保存物理层、原生装配和释放壁行独立复读")
    a=p.parse_args();job=a.output.resolve();work=a.run.resolve()
    if a.reread:
        job.relative_to(ROOT/"outputs");job.mkdir(parents=True,exist_ok=False)
        with exclusive(job):publish(job,work,True)
    else:publish(job,work)
