"""在同一物理时刻、同一原面编号下比较热敏珠的一整步与两半步，并生成轻量对照包。"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import zipfile

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE, ROOT, RUN, input_file, observe, read, sha
from thermistorStartupSurface import columns,write


def difference(delta,area,scale):
    peak=int(np.argmax(abs(delta)))
    rms=float(np.sqrt(np.average(delta**2,weights=area)))
    return dict(weightedRMS=rms,maximumAbsolute=float(abs(delta[peak])),maxDifferenceFaceID=peak,
                signedAtMax=float(delta[peak]),scale=scale,rmsFraction=rms/scale,
                maximumFraction=float(abs(delta[peak])/scale))


def compare(first,second,reference,out):
    mid=read(first);end=read(second);fine=end["surface"]
    assert mid["accepted"] and end["accepted"] and mid["prefixRows"]==240
    assert sha(reference)=="25bdbfa8f11b6a12a6488b81b41002e7c3bc9808f8490f4b63987ee782c54b61"
    with zipfile.ZipFile(reference) as z:
        coarse=json.loads(z.read("summary.json"))
        coarse_bytes=z.read("wall.csv")
        body=z.read("body.json")
    assert sha(input_file("body.json"))==hashlib.sha256(body).hexdigest()==coarse["bodySHA256"]==fine["bodySHA256"]
    assert sha(input_file("mesh.su2"))==coarse["meshSHA256"]==fine["meshSHA256"]
    base=coarse["energy"];one=mid["energy"];two=end["energy"]
    assert base["oldFrameSHA256"]==one["oldFrameSHA256"]==sha(input_file("seed_00000.csv"))
    assert one["acceptedFrameSHA256"]==two["oldFrameSHA256"]
    assert one["dtSeconds"]==two["dtSeconds"] and 2*two["dtSeconds"]==base["dtSeconds"]==end["elapsedPhysicalTimeS"]
    observe_dir=second.parent.parent/"second_observe"
    n=read(second.parent.parent/"second_state/verified.json")["innerRows"]
    calls=(observe_dir/"readout_updates.csv").read_text("ascii").splitlines()
    assert calls==[f"2,{n},flow",f"2,{n},scalar"]
    sample=columns(observe_dir/"before_wall.csv",["id","p","u","v","w","theta"])
    ordered=sample[np.argsort(sample[:,0])]
    scratch=columns(observe_dir/"restart_00002.csv",["PointID","Pressure","Velocity_x","Velocity_y","Velocity_z","Species_0"])
    np.testing.assert_array_equal(scratch[:,0],ordered[:,0])
    scratch_error=np.max(abs(scratch[:,1:]-ordered[:,1:]),axis=0)
    assert scratch_error.max()<2e-12
    assert abs(base["relativeImbalance"])<.005 and abs(one["relativeImbalance"])<.005 and abs(two["relativeImbalance"])<.005
    assert hashlib.sha256(coarse_bytes).hexdigest()==coarse["wallTableSHA256"]
    wall=second.parent/"wall.csv"
    assert sha(wall)==fine["wallTableSHA256"]
    a=np.genfromtxt(io.BytesIO(coarse_bytes),delimiter=",",names=True)
    b=np.genfromtxt(str(wall),delimiter=",",names=True)
    assert a.dtype.names==b.dtype.names and np.array_equal(a["faceID"],np.arange(5216))
    assert np.array_equal(a["faceID"],b["faceID"]) and np.array_equal(a["areaM2"],b["areaM2"])
    area=a["areaM2"];assert np.all(area>0)
    geom=read(input_file("geometry.json"))
    dyn=geom["rho"]*geom["speedMps"]**2/2
    p=difference(b["pressureGaugePa"]-a["pressureGaugePa"],area,dyn)
    q=difference(b["heatFluxWm2"]-a["heatFluxWm2"],area,float(np.max(abs(a["heatFluxWm2"]))))
    first_layer=difference(b["nearWallCellTemperatureK"]-a["nearWallCellTemperatureK"],area,geom["nominalWallK"]-geom["airK"])
    for rows,summary in ((a,coarse),(b,fine)):
        assert abs(np.dot(area,rows["heatFluxWm2"])/summary["powerIntoWallW"]-1)<1e-11
        assert np.argmin(rows["heatFluxWm2"])==summary["peakOutwardFaceID"]
    rel_power=(two["wallReleaseW"]-base["wallReleaseW"])/base["wallReleaseW"]
    passed=p["rmsFraction"]<.01 and q["rmsFraction"]<.01 and max(p["maximumFraction"],q["maximumFraction"])<.05 and abs(rel_power)<.01
    summary=dict(elapsedPhysicalTimeS=end["elapsedPhysicalTimeS"],originalFaces=len(area),areaM2=float(area.sum()),
                 initialSHA256=one["oldFrameSHA256"],referenceFrameSHA256=base["acceptedFrameSHA256"],
                 intermediateFrameSHA256=one["acceptedFrameSHA256"],fineFrameSHA256=two["acceptedFrameSHA256"],
                 fineWallTableSHA256=sha(wall),referenceWallTableSHA256=coarse["wallTableSHA256"],bodySHA256=fine["bodySHA256"],
                 meshSHA256=fine["meshSHA256"],firstHalfEnergy=one,secondHalfEnergy=two,pressurePa=p,heatFluxWm2=q,
                 firstLayerK=first_layer,referenceWallReleaseW=base["wallReleaseW"],relativePowerDifference=rel_power,
                 referencePeakFaceID=coarse["peakOutwardFaceID"],finePeakFaceID=fine["peakOutwardFaceID"],
                 referencePeakFluxWm2=float(a["heatFluxWm2"].min()),finePeakFluxWm2=float(b["heatFluxWm2"].min()),
                 referencePeakCentroidM=coarse["peakFaceCentroidM"],finePeakCentroidM=fine["peakFaceCentroidM"],
                 peakDisplacementM=float(np.linalg.norm(np.array(fine["peakFaceCentroidM"])-coarse["peakFaceCentroidM"])),
                 withinComparisonThresholds=bool(passed),jointExperimentVerifiedM2=0,
                 appendedUpdatesSkipped=calls,scratchFrameMaxDifference=scratch_error.tolist())
    write(out/"summary.json",summary)
    with zipfile.ZipFile(out/"startup_6ms_halfsteps.zip","x") as z:
        for file in (out/"summary.json",wall):
            info=zipfile.ZipInfo(file.name,date_time=(2026,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,file.read_bytes())


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--root",type=Path,default=RUN)
    cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN)
    first=area/"first_energy/accepted.json";second=area/"second_energy/accepted.json"
    reference=ROOT/"cases/thermistor_b/startup_6ms.zip";out=area/"compare"
    if a.worker:
        compare(first,second,reference,out)
    else:
        out.mkdir(parents=True,exist_ok=False)
        paths=[HERE/"PLAN.md",Path(__file__).resolve(),HERE/"run_step.py",first,second,second.parent/"wall.csv",
               reference,input_file("body.json"),input_file("geometry.json"),input_file("mesh.su2"),input_file("seed_00000.csv"),
               area/"second_state/verified.json",area/"second_observe/readout_updates.csv",
               area/"second_observe/restart_00002.csv",area/"second_observe/before_wall.csv"]
        rec=observe([sys.executable,str(Path(__file__).resolve()),"--root",str(area),"--worker"],out,
                    {str(p):sha(p) for p in paths},120,100*1024**2)
        assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
