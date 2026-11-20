"""将两个官方独立壁反力读数之差保守映回热敏珠原面，核对微小热流变化。"""
import io
import os
from pathlib import Path
import sys
import zipfile

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE,ROOT,RUN,input_file,observe,read,sha
from beadSurfaceResult import projection
from thermistorStartupSurface import columns,write


def check(coarse,fine,package,wall,out):
    ref=read(HERE.parent/"v290/evidence/surface.json")
    trial=read(wall.parent/"accepted.json")
    assert sha(coarse)==ref["energy"]["reactionCheckFrameSHA256"]
    assert sha(fine)==trial["energy"]["reactionCheckFrameSHA256"]
    assert sha(wall)==trial["surface"]["wallTableSHA256"]
    with zipfile.ZipFile(package) as z:
        base=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True)
    new=np.genfromtxt(str(wall),delimiter=",",names=True)
    assert np.array_equal(new["areaM2"],base["areaM2"]) and np.array_equal(new["faceID"],base["faceID"])
    a=new["areaM2"];published=new["heatFluxWm2"]-base["heatFluxWm2"]
    body=read(input_file("body.json"));g=read(input_file("geometry.json"))
    with np.load(input_file("nodes.npz")) as grid:
        weights,ids=projection(np.array(body["V"]),np.array(body["F"])-1,grid["ids"][:,:65],np.linspace(0,2*np.pi,49))
    fields=["PointID","Residual_Species_0"]
    c=columns(coarse,fields);f=columns(fine,fields)
    assert np.array_equal(f[:,0],np.arange(len(f))) and np.array_equal(c[:,0],f[:,0])
    scale=g["rho"]*g["cp"]*g["speedMps"]*g["transverseDiameterM"]**2*(g["nominalWallK"]-g["airK"])
    nodal=(f[ids,1]-c[ids,1])*scale
    independent=weights@(nodal/(weights.T@a))
    assert abs(a@independent-nodal.sum())<1e-16
    err=independent-published
    write(out/"difference.json",dict(faces=len(a),wallNodes=len(ids),units="W/m2",independentMaxAbs=float(max(abs(independent))),
          independentWeightedRMS=float(np.sqrt(np.average(independent**2,weights=a))),
          maxAbsDifferenceFromPublishedDelta=float(max(abs(err))),weightedRMSDifferenceFromPublishedDelta=float(np.sqrt(np.average(err**2,weights=a))),
          publishedReleasePowerChangeW=float(-a@published),independentReleasePowerChangeW=float(-a@independent),
          coarseIndependentFrameSHA256=sha(coarse),fineIndependentFrameSHA256=sha(fine)))


if __name__=="__main__":
    area=RUN/"pair";out=area/"independent_delta"
    coarse=ROOT/"outputs/v290_validation/release01/restart_00001.csv";fine=area/"second_release/restart_00002.csv"
    package=ROOT/"cases/thermistor_b/startup_6ms.zip";wall=area/"second_energy/wall.csv"
    if "--worker" in sys.argv:
        check(coarse,fine,package,wall,out)
    else:
        out.mkdir(exist_ok=False)
        files=[HERE/"PLAN.md",Path(__file__).resolve(),coarse,fine,package,wall,wall.parent/"accepted.json",
               HERE.parent/"v290/evidence/surface.json",ROOT/"scripts/beadSurfaceResult.py",
               input_file("body.json"),input_file("geometry.json"),input_file("nodes.npz")]
        rec=observe([sys.executable,str(Path(__file__).resolve()),"--worker"],out,{str(p):sha(p) for p in files},120,100*1024**2)
        assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
