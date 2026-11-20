"""在相同物理终点按原三角面比较 48 层与 56 层，并单独核对相同位置的流体温度。"""
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
from advance import GRID,HERE,ROOT,RUN,grid_info,observe,read,sha
from thermistorStartupSurface import write

BASE=ROOT/"cases/thermistor_b/startup_6ms_halfsteps.zip"
BASE_SHA="9a1ce5800c51a19b845877b83c735cd56d082a873c04024194f1f666eb1f84f3"


def metrics(delta,area,scale):
    peak=int(np.argmax(abs(delta)))
    rms=float(np.sqrt(np.average(delta**2,weights=area)))
    return dict(weightedRMS=rms,maximumAbsolute=float(abs(delta[peak])),maxDifferenceFaceID=peak,
                signedAtMax=float(delta[peak]),scale=scale,rmsFraction=rms/scale,maximumFraction=float(abs(delta[peak])/scale))


def compare(area,out):
    g=grid_info();first=read(area/"first_energy/accepted.json");end=read(area/"second_energy/accepted.json")
    assert first["accepted"] and end["accepted"] and first["verifiedPrefixRows"]==560
    one,two=first["energy"],end["energy"];fine=end["surface"]
    assert sha(BASE)==BASE_SHA
    with zipfile.ZipFile(BASE) as z:
        reference=json.loads(z.read("summary.json"));raw=z.read("wall.csv")
    assert hashlib.sha256(raw).hexdigest()==reference["fineWallTableSHA256"]
    assert reference["bodySHA256"]==fine["bodySHA256"] and g["gridSHA256"]==fine["meshSHA256"]
    assert reference["meshSHA256"]!=fine["meshSHA256"] and reference["initialSHA256"]==g["sourceSeedSHA256"]
    assert one["oldFrameSHA256"]==g["seedSHA256"] and two["oldFrameSHA256"]==one["acceptedFrameSHA256"]
    assert one["dtSeconds"]==two["dtSeconds"]==reference["secondHalfEnergy"]["dtSeconds"]
    assert end["elapsedPhysicalTimeS"]==2*one["dtSeconds"]==reference["elapsedPhysicalTimeS"]
    wall=area/"second_energy/wall.csv";temperature=area/"second_energy/matched_temperature.csv"
    assert sha(wall)==fine["wallTableSHA256"] and sha(temperature)==end["matchedTableSHA256"]
    c=np.genfromtxt(io.BytesIO(raw),names=True,delimiter=",")
    f=np.genfromtxt(str(wall),names=True,delimiter=",")
    t=np.genfromtxt(str(temperature),names=True,delimiter=",")
    assert c.dtype.names==f.dtype.names and np.array_equal(c["faceID"],np.arange(5216))
    assert np.array_equal(c["faceID"],f["faceID"]) and np.array_equal(t["faceID"],f["faceID"])
    a=f["areaM2"];assert np.array_equal(a,c["areaM2"]) and np.all(a>0)
    base_energy=reference["secondHalfEnergy"]
    for rows,e in ((c,base_energy),(f,two)):
        assert abs(a@rows["heatFluxWm2"]/e["wallReleaseW"]+1)<1e-11
    pressure=metrics(f["pressureGaugePa"]-c["pressureGaugePa"],a,21.24)
    heat=metrics(f["heatFluxWm2"]-c["heatFluxWm2"],a,float(max(abs(c["heatFluxWm2"]))))
    matched=metrics(t["temperatureAtOriginalFirstLayerK"]-c["nearWallCellTemperatureK"],a,10.)
    power=two["wallReleaseW"]/base_energy["wallReleaseW"]-1
    assert int(np.argmin(c["heatFluxWm2"]))==reference["finePeakFaceID"]
    assert int(np.argmin(f["heatFluxWm2"]))==fine["peakOutwardFaceID"]
    result=dict(originalFaces=len(a),areaM2=float(a.sum()),elapsedPhysicalTimeS=end["elapsedPhysicalTimeS"],
        bodySHA256=fine["bodySHA256"],fineSurface=fine,firstHalfEnergy=one,referenceMeshSHA256=reference["meshSHA256"],
        referenceFrameSHA256=reference["fineFrameSHA256"],referenceWallTableSHA256=reference["fineWallTableSHA256"],
        referenceInitialSHA256=reference["initialSHA256"],gridTransfer=g,matchedTableSHA256=sha(temperature),
        pressurePa=pressure,heatFluxWm2=heat,matchedTemperatureK=matched,referenceWallReleaseW=base_energy["wallReleaseW"],
        relativePowerDifference=power,referencePeakFaceID=reference["finePeakFaceID"],referencePeakCentroidM=reference["finePeakCentroidM"],
        peakDisplacementM=float(np.linalg.norm(np.array(fine["peakFaceCentroidM"])-reference["finePeakCentroidM"])),
        withinComparisonThresholds=bool(max(pressure["rmsFraction"],heat["rmsFraction"])<.01 and
            max(pressure["maximumFraction"],heat["maximumFraction"])<.05 and abs(power)<.01),jointExperimentVerifiedM2=0)
    summary=out/"summary.json";write(summary,result)
    with zipfile.ZipFile(out/"startup_6ms_nearwall56.zip","x") as archive:
        for p in (summary,wall,temperature):
            item=zipfile.ZipInfo(p.name,date_time=(2026,1,1,0,0,0));item.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(item,p.read_bytes())


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("--root",type=Path,default=RUN)
    cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);out=area/"compare"
    if a.worker:compare(area,out)
    else:
        out.mkdir(exist_ok=False)
        files=[HERE/"PLAN.md",Path(__file__).resolve(),HERE/"advance.py",BASE,GRID/"grid.json",
               area/"first_energy/accepted.json",area/"second_energy/accepted.json",area/"second_energy/wall.csv",
               area/"second_energy/matched_temperature.csv"]
        rec=observe([sys.executable,str(Path(__file__).resolve()),"--root",str(area),"--worker"],out,
                    {str(p):sha(p) for p in files},180,200*1024**2)
        assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
