"""核对同一热网、同一初温、同一时刻的原面面积、守恒及动量网格敏感性。"""
import argparse
import json
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED,ZipFile

from refine_flow import BASE,DOC,OLD,ROOT,SEED,UP,geometry,read,save,sha
import numpy as np
from finalize import read_step
from flux_transfer import build_model
from heat_step import surface_pressure
sys.path.insert(0,str(DOC.parent/"v355"))
from compare import stored_power


def compare(work,dest):
    geo,tri,diameter,_=geometry()
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    folder=ROOT/"cases/simple_launcher/wake";seal=read(folder/"receipt.json")
    assert sha(folder/"fine_primary.json")==seal["fineSummarySHA256"]
    assert sha(folder/"wake_heat_grid.zip")==seal["wakeHeatGridSHA256"]
    base=read(folder/"fine_primary.json");condition=base["condition"]
    assert condition==read(folder/"condition.json") and base["originalFaces"]==len(area)
    with ZipFile(folder/"wake_heat_grid.zip") as z:
        import hashlib
        assert hashlib.sha256(z.read("fine_wall.csv")).hexdigest()==base["wallSHA256"]
        old=np.loadtxt(z.open("fine_wall.csv"),delimiter=",",skiprows=1)
    dt=condition["physicalStepSizeS"]
    scale=(2*np.pi*diameter**3*condition["densityKgM3"]*condition["specificHeatJkgK"]*
           (condition["wallTemperatureK"]-condition["ambientTemperatureK"]))
    initial=SEED/"initial/fields.npz"
    assert sha(initial)==base["initialThermalSHA256"] and sha(SEED/"steady/temperature.npz")==base["heatMeshSHA256"]
    with np.load(SEED/"steady/temperature.npz") as src:
        mesh={key:src[key] for key in ("xy","triangles")}
    records=[]
    for n in (1,2):
        table,record=read_step(work/"flow",work/"heat","disturbed","half",n,area,geo)
        h,f=record["heat"],record["flow"]
        assert h["heatMeshSHA256"]==base["heatMeshSHA256"] and sha(initial)==h["initialThetaSHA256"]
        assert abs(f["elapsedPhysicalS"]-n*dt)<1e-12
        if n==2:assert f["seedSHA256"]==records[0]["flow"]["frameSHA256"]
        new=work/"heat"/f"disturbed_half_{n}"/"thermal.npz"
        storage=stored_power(initial,new,mesh,dt,scale)
        assert abs(storage-h["storageW"])<1e-14 and abs(h["outwardW"]-h["farW"]-storage)<1e-12
        record["independentStorageW"]=storage;records.append(record);initial=new
    assert condition["elapsedPhysicalTimeS"]==h["elapsedPhysicalS"]
    np.testing.assert_array_equal(table[:,:2],old[:,:2])
    assert np.isfinite(table).all() and ((table[:,4]>=300-1e-12)&(table[:,4]<=300.001+1e-12)).all()
    dp=table[:,2]-old[:,2];dq=table[:,3]-old[:,3];qinf=.5*condition["densityKgM3"]*.05**2
    def p_rms(x):
        return float(np.sqrt((area@x**2)/area.sum())/qinf)
    region={}
    for group in geo["regions"]:
        name=group["name"];s=slice(group["firstFace"]-1,group["lastFace"])
        p0=float(-area[s]@old[s,3]);p1=float(-area[s]@table[s,3])
        assert abs(p0-base["regionW"][name])<1e-13 and abs(p1-h["regionW"][name])<1e-13
        region[name]=dict(originalW=p0,refinedW=p1,relativeDifference=abs(p0/p1-1),
                          pressureRMSOverDynamic=float(np.sqrt((area[s]@dp[s]**2)/area[s].sum())/qinf),
                          heatFluxRelativeRMS=float(np.sqrt((area[s]@dq[s]**2)/(area[s]@table[s,3]**2))),
                          pressureSquaredDifferenceIntegral=float(area[s]@dp[s]**2),
                          fluxSquaredDifferenceIntegral=float(area[s]@dq[s]**2))
        if name=="sphericalNose":
            q0=-old[s,3];q1=-table[s,3]
            peak=dict(originalWm2=float(q0.max()),refinedWm2=float(q1.max()),relativeDifference=float(abs(q0.max()/q1.max()-1)),
                      originalFace=int(old[s,0][q0.argmax()]),refinedFace=int(table[s,0][q1.argmax()]))
    # 同一加密压力帧另取原有节点，重做原三顶点/底盘积分，单列原面取样细化的影响。
    with np.load(work/"prepared/meridian.npz") as src:old_ids=src["oldToNew"]
    native=np.loadtxt(work/"flow/disturbed_half_2/restart_00002.csv",delimiter=",",skiprows=1)
    m=build_model(UP,OLD,SEED/"steady")
    restricted,_=surface_pressure(m,native[old_ids])
    errors=dict(pressureRMSOverDynamic=p_rms(dp),sameNodesPressureRMSOverDynamic=p_rms(restricted-old[:,2]),
                samplingOnlyPressureRMSOverDynamic=p_rms(table[:,2]-restricted),
                fluxRelativeRMS=float(np.sqrt((area@dq**2)/(area@table[:,3]**2))),
                outwardDifference=float(abs(base["outwardW"]/h["outwardW"]-1)),region=region,nosePeak=peak,
                maxPressureChangePa=float(dp[np.argmax(abs(dp))]),maxPressureChangeFace=int(table[np.argmax(abs(dp)),0]),
                maxFluxChangeWm2=float(dq[np.argmax(abs(dq))]),maxFluxChangeFace=int(table[np.argmax(abs(dq)),0]))
    gates=dict(pressure=errors["pressureRMSOverDynamic"]<.05,heatFlux=errors["fluxRelativeRMS"]<.05,
               power=errors["outwardDifference"]<.05,sixRegions=max(v["relativeDifference"] for v in region.values())<.10,
               nosePeak=peak["relativeDifference"]<.10)
    g=read(work/"prepared/refinement.json")
    summary=dict(condition=condition,faces=len(area),areaM2=float(area.sum()),surfaceSHA256=geo["fileSHA256"],
                 heatMeshSHA256=h["heatMeshSHA256"],initialThermalSHA256=base["initialThermalSHA256"],
                 originalFlowSHA256=base["sourceFrameSHA256"],originalThermalSHA256=base["thermalFrameSHA256"],
                 originalWallSHA256=base["wallSHA256"],originalFlowMeshSHA256=base["flowMeshSHA256"],
                 refinedFlowSHA256=f["frameSHA256"],refinedThermalSHA256=h["fieldSHA256"],refinedWallSHA256=h["wallSHA256"],
                 refinedFlowMeshSHA256=f["meshSHA256"],refinedSeedSHA256=records[0]["flow"]["seedSHA256"],
                 refinedNodes=g["nodes"],originalFirstNormalMedianMm=g["oldResolution"]["firstNormalMedianMm"],
                 refinedFirstNormalMedianMm=g["refinedResolution"]["firstNormalMedianMm"],
                 outwardW=h["outwardW"],farW=h["farW"],storageW=h["storageW"],
                 pressureRangePa=[table[:,2].min(),table[:,2].max()],fluxRangeWm2=[table[:,3].min(),table[:,3].max()],
                 errors=errors,gates=gates,allSpatialGatesPassed=all(gates.values()),jointExperimentalAreaM2=0.)
    dest.mkdir(parents=True,exist_ok=False);bundle=dest/"momentum_grid.zip"
    with ZipFile(bundle,"x",ZIP_DEFLATED) as z:
        z.writestr("comparison.json",json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False))
        z.write(work/"heat/disturbed_half_2/wall.csv","wall.csv")
    save(dest/"fields.json",dict(summary=summary,steps=records,sampling=read(work/"sampling.json"),bundleSHA256=sha(bundle)))
    print(json.dumps(errors,indent=2))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work",type=Path);ap.add_argument("dest",type=Path);a=ap.parse_args()
    work=a.work.resolve();work.relative_to(BASE);compare(work,a.dest)
