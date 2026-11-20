"""诊断已拒绝装配的时间层，只读既有全场；保持原验收状态。"""
import argparse
import csv
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE, ROOT, RUN, input_file, observe, read, sha
from thermistorStartupSurface import checked, columns, write


def inspect(area,out):
    native=area/"first01";obs=area/"first_observe";ind=area/"first_release"
    proof=read(area/"first_state/verified.json");assert proof["passed"] and proof["prefixRows"]==240
    for folder in (native,obs,ind):checked(folder)
    sample=columns(obs/"before_wall.csv",["id","rhs","theta","theta_old","volume","p","u","v","w","dt"])
    sample=sample[np.argsort(sample[:,0])]
    coord=columns(input_file("seed_00000.csv"),["x","y","z","Species_0"])
    ids=np.arange(len(coord))
    np.testing.assert_array_equal(sample[:,0],ids)
    assert np.all(sample[:,9]==.025)
    with np.load(input_file("nodes.npz")) as z:
        wall=np.unique(z["ids"][:,:65])
    fields=["PointID","Pressure","Velocity_x","Velocity_y","Velocity_z","Species_0","Pressure_Coefficient"]
    result={}
    for label,folder in (("official260",native),("observer261",obs),("officialReleaseOne",ind)):
        new=columns(folder/"restart_00001.csv",fields)
        np.testing.assert_array_equal(new[:,0],ids)
        dp=sample[:,5]-new[:,1]
        peak=int(np.argmax(abs(dp)))
        result[label]=dict(maxAbsPressure=float(abs(dp[peak])),minPressureDifference=float(dp.min()),
            maxPressureDifference=float(dp.max()),medianPressureDifference=float(np.median(dp)),
            maxWallPressureDifference=float(max(abs(dp[wall]))),nodesOverPressureLimit=int(np.sum(abs(dp)>1e-8)),
            wallNodesOverPressureLimit=int(np.sum(abs(dp[wall])>1e-8)),peakNode=peak,peakXYZOverD=coord[peak,:3].tolist(),
            maxVelocityDifference=float(np.max(abs(sample[:,6:9]-new[:,2:5]))),
            maxCurrentScalarDifference=float(np.max(abs(sample[:,2]-new[:,5]))),
            maxPrimitivePressureVsHalfCp=float(np.max(abs(sample[:,5]-.5*new[:,6]))))
    after=columns(obs/"after_wall.csv",["id","rhs"]);after=after[np.argsort(after[:,0])]
    np.testing.assert_array_equal(after[:,0],ids)
    reaction=sample[:,1]-after[:,1]
    assert np.all(after[wall,1]==0) and np.count_nonzero(reaction)==np.count_nonzero(reaction[wall])
    independent=columns(ind/"restart_00001.csv",["PointID","Residual_Species_0"])
    np.testing.assert_array_equal(independent[:,0],ids)
    e=read(input_file("geometry.json"));scale=e["rho"]*e["cp"]*e["speedMps"]*e["transverseDiameterM"]**2*(e["nominalWallK"]-e["airK"])
    stages={key:float(value) for i,key,value in csv.reader((obs/"stages.csv").read_text("ascii").splitlines()) if int(i)==proof["innerRows"]}
    q=(stages["beforeWall"]-stages["afterWall"])*scale
    far=(stages["afterFar"]-stages["beforeFar"])*scale;storage=(stages["afterTime"]-stages["beforeTime"])*scale
    output=dict(acceptedForSameFrameHeat=False,threshold=1e-8,nodes=len(ids),wallNodes=len(wall),
                oldScalarMaxDifference=float(np.max(abs(sample[:,3]-coord[:,3]))),snapshots=result,
                observerOnly=dict(wallReleaseW=q,farOutW=far,storageW=storage,relativeImbalance=(q-far-storage)/q,
                                  independentWallReactionErrorOverPeak=float(np.max(abs(reaction[wall]+independent[wall,1]))/np.max(abs(reaction[wall])))),
                frameHashes={k:sha(v/"restart_00001.csv") for k,v in (("official260",native),("observer261",obs),("officialReleaseOne",ind))})
    assert result["official260"]["maxAbsPressure"]>output["threshold"]
    write(out/"assembly.json",output)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--root",type=Path,default=RUN);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);out=area/"assembly_check"
    if a.worker:
        inspect(area,out)
    else:
        out.mkdir(parents=True,exist_ok=False)
        paths=[Path(__file__).resolve(),HERE/"PLAN.md",area/"first_state/verified.json",input_file("seed_00000.csv"),input_file("nodes.npz"),
               input_file("geometry.json"),ROOT/"scripts/thermistorStartupSurface.py"]
        paths.extend(area/d/n for d in ("first01","first_observe","first_release") for n in ("debugger_process.json","restart_00001.csv"))
        paths.extend(area/"first_observe"/n for n in ("before_wall.csv","after_wall.csv","stages.csv"))
        rec=observe([sys.executable,str(Path(__file__).resolve()),"--root",str(area),"--worker"],out,
                    {str(p):sha(p) for p in paths},120,100*1024**2)
        assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
