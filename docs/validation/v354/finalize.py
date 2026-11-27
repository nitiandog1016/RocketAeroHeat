"""从六个独立动量和显热终态核算同刻误差，打包完整低速研究快照。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZIP_DEFLATED,ZipFile

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from simpleLauncherTunnel import CASE,geometry


def read(path):
    return json.loads(path.read_bytes())


def finalize(flow,heat,dest):
    definition,tri,d,_=geometry()
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    tables={};records={};comparisons={};files={}
    for branch in ("control","disturbed"):
        for size,step in (("full",1),("half",1),("half",2)):
            name=f"{branch}_{size}_{step}";f=flow/name;t=heat/name
            a=read(f/"step.json");b=read(t/"heat.json")
            for path in (f,t):
                rec=read(path/"debugger_process.json")
                assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
                assert not (rec["stopReason"] or rec["survivingProcessIDs"] or rec["supervisorError"])
            assert a["accepted"] and max(a["last20Maximum"].values())< -8 and b["accepted"]
            for field in ("branch","stepSize","step","elapsedPhysicalS"):assert a[field]==b[field]
            assert (branch,size,step)==(a["branch"],a["stepSize"],a["step"])
            assert sha(f/f"restart_{step:05d}.csv")==b["sourceFrameSHA256"]==a["frameSHA256"]
            assert sha(f/f"history_{step:05d}.csv")==a["historySHA256"]
            assert sha(t/"thermal.npz")==b["fieldSHA256"] and sha(t/"wall.csv")==b["wallSHA256"]
            assert b["relativeEnergyImbalance"]<.005 and b["originalResidual"]<1e-8 and b["normalizedDivergence"]<1e-8
            assert b["relativeProjection"]<.05 and b["minTheta"]>=-1e-12 and b["maxTheta"]<=1+1e-12
            table=np.loadtxt(t/"wall.csv",delimiter=",",skiprows=1)
            np.testing.assert_array_equal(table[:,0],np.arange(1,len(area)+1))
            np.testing.assert_allclose(table[:,1],area,rtol=0,atol=1e-15)
            assert np.isfinite(table).all() and table[:,2].min()>0
            assert abs(-area@table[:,3]-b["outwardW"])<1e-13
            if step==2:
                prior=records[f"{branch}_half_1"]
                assert prior["flow"]["frameSHA256"]==a["seedSHA256"] and prior["heat"]["fieldSHA256"]==b["initialThetaSHA256"]
            elif size=="half":
                prior=records[f"{branch}_full_1"]
                assert prior["flow"]["seedSHA256"]==a["seedSHA256"] and prior["heat"]["initialThetaSHA256"]==b["initialThetaSHA256"]
            for group in definition["regions"]:
                s=slice(group["firstFace"]-1,group["lastFace"])
                assert abs(-area[s]@table[s,3]-b["regionW"][group["name"]])<1e-13
            files[name]=t/"wall.csv";tables[name]=table
            records[name]=dict(flow=a,heat={k:v for k,v in b.items() if k!="baseTapGaugeMilliPa"},
                heatSeconds=rec["elapsedSeconds"],heatPeakRSSBytes=rec["peakFamilyRSSBytes"])
        full=tables[f"{branch}_full_1"];half=tables[f"{branch}_half_2"]
        r=records[f"{branch}_full_1"];s=records[f"{branch}_half_2"]
        assert r["flow"]["elapsedPhysicalS"]==s["flow"]["elapsedPhysicalS"]
        dp=full[:,2]-half[:,2];dq=full[:,3]-half[:,3]
        comparison=dict(pressureRMSOverDynamic=float(np.sqrt(np.average(dp**2,weights=area))/(.5*1.177*.05**2)),
            heatFluxRelativeRMS=float(np.sqrt((area@dq**2)/(area@half[:,3]**2))),
            outwardDifference=float(abs(r["heat"]["outwardW"]/s["heat"]["outwardW"]-1)),
            regionPowerDifferences={k:abs(v/s["heat"]["regionW"][k]-1) for k,v in r["heat"]["regionW"].items()})
        comparisons[branch]=comparison
        assert comparison["pressureRMSOverDynamic"]<.02 and comparison["heatFluxRelativeRMS"]<.05,comparison
        assert comparison["outwardDifference"]<.02 and max(comparison["regionPowerDifferences"].values())<.05,comparison
    condition=read(CASE/"tunnel.json")
    condition.update(caseLabel="简化箭体 · 低速风洞尾迹扰动",
        elapsedPhysicalTimeS=.125*d/.05,physicalStepSizeS=.0625*d/.05,
        timeCaption="距构造初场 0.600 s",initialVelocityPerturbationFraction=.01,
        scope="0.05 m/s 恒物性无浮力风洞，14D 粗网、全周启动瞬态，给定等温壁。实际尾迹和实验精度待验证。")
    ref=records["disturbed_half_2"];b=ref["heat"];table=tables["disturbed_half_2"]
    summary=dict(accepted=True,condition=condition,surfaceSHA256=definition["fileSHA256"],originalFaces=len(area),
        areaM2=float(area.sum()),wallSHA256=b["wallSHA256"],outwardW=b["outwardW"],regionW=b["regionW"],
        sourceFrameSHA256=ref["flow"]["frameSHA256"],thermalFrameSHA256=b["fieldSHA256"],
        initialFlowSHA256=records["disturbed_half_1"]["flow"]["seedSHA256"],
        initialThermalSHA256=records["disturbed_half_1"]["heat"]["initialThetaSHA256"],
        flowMeshSHA256=ref["flow"]["meshSHA256"],storageW=b["storageW"],farW=b["farW"],
        pressureRangePa=[float(table[:,2].min()),float(table[:,2].max())],
        heatFluxRangeWm2=[float(table[:,3].min()),float(table[:,3].max())],
        comparison=comparisons,verifiedByExperimentM2=0)
    dest.mkdir(parents=True,exist_ok=False)
    (dest/"condition.json").write_bytes((json.dumps(condition,ensure_ascii=False,indent=2)+"\n").encode())
    zpath=dest/"wind_tunnel.zip"
    with ZipFile(zpath,"x",ZIP_DEFLATED) as z:
        z.writestr("summary.json",json.dumps(summary,ensure_ascii=False,allow_nan=False,indent=2))
        z.write(files["disturbed_half_2"],"wall.csv")
    (dest/"receipt.json").write_bytes((json.dumps(dict(sha256=sha(zpath),faces=len(area),outwardW=b["outwardW"]),indent=2)+"\n").encode())
    (dest/"study.json").write_bytes((json.dumps(dict(records=records,comparisons=comparisons),indent=2,allow_nan=False)+"\n").encode())
    return summary


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("flow","heat","dest"):ap.add_argument(key,type=Path)
    a=ap.parse_args();result=finalize(a.flow.resolve(),a.heat.resolve(),a.dest.resolve())
    print(result["outwardW"],result["comparison"])
