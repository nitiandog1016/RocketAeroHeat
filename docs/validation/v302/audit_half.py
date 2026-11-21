"""对独立半温差场逐面核对壁热、压力、首层温度、七区功率及反求壁温。"""
import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

from half_case import BASE, HERE, ROOT, read, sha, thermal_frame
from thermistorBead import CP,D,DT,P0,RHO,T0,U
from thermistorBeadSector import COPIES
sys.path.insert(0,str(HERE.parent/"v285"))
from audit_temperature import reaction


def save(path,data):
    with path.open("x",encoding="utf8") as f:json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)


def audit(folder,out):
    rec,res=thermal_frame(folder);plus,res20=thermal_frame(folder/"plus20");source=read(folder/"inputs.json")
    assert all(sha(folder/k)==h for k,h in source["hashes"].items())
    base=read(HERE.parent/"v300/evidence/summary.json")
    assert sha(BASE/"surface/wall.csv")==base["wallSHA256"]==source["referenceTable"]
    assert sha(BASE/"restart.csv")==base["sourceFrameSHA256"]==source["referenceFrame"]
    assert sha(folder/"case.json")==rec["inputSHA256"][str(folder/"case.json")] and read(folder/"case.json")["wallK"]==299.15
    cfg=(folder/"case.cfg").read_text("ascii")
    assert cfg==(BASE/"case.cfg").read_text("ascii").replace("(wall, VALUE, 1.0)","(wall, VALUE, 0.5)").replace("MARKER_ISOTHERMAL= (wall, 2)","MARKER_ISOTHERMAL= (wall, 1.5)")
    assert (folder/"plus20/case.cfg").read_text("ascii")==cfg.replace("ITER= 7000\n","ITER= 20\n")
    assert sha(folder/"plus20/initial_flow.csv")==sha(folder/"restart.csv")
    with np.load(folder/"nodes.npz") as z:xyz,wall,near,area=(z[k] for k in ("xyz","wallMap","firstMap","area"))
    walls=np.unique(wall);initial=np.loadtxt(BASE/"initial_flow.csv",delimiter=",",skiprows=1)
    initial[walls,8]-=.5;initial[walls,9]*=.5
    np.testing.assert_array_equal(initial,np.loadtxt(folder/"initial_flow.csv",delimiter=",",skiprows=1))
    history=np.loadtxt(folder/"history.csv",delimiter=",",skiprows=1,usecols=2,dtype=int);last=int(history[-1])
    steps=[*range(500,last,500),last];frames={}
    for name in ["initial_flow.csv",*(f"restart_{i:06d}.csv" for i in steps),"plus20/restart.csv"]:
        path=folder/name
        state=np.loadtxt(path,delimiter=",",skiprows=1)
        assert state.shape==(len(xyz),10) and np.isfinite(state).all()
        np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)))
        np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-13)
        theta=state[:,9];assert theta.min()>=-1e-12 and theta.max()<=.5+1e-12
        assert np.all(theta[walls]==.5) and np.max(abs(state[walls,5:8]))<1e-10 and (P0+RHO*U**2*state[:,4]).min()>0
        frames[name]=dict(sha256=sha(path),thetaRange=[float(theta.min()),float(theta.max())])
        if name.startswith("restart_"):
            vtk=name.replace("restart_","surface_").replace(".csv",".vtk")
            assert frames[name]["sha256"]==rec["outputSHA256"][name] and sha(folder/vtk)==rec["outputSHA256"][vtk]
            frames[name]["surfaceSHA256"]=rec["outputSHA256"][vtk]
    end_state=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1)
    rhs=reaction(folder,folder/"reaction0")*(COPIES*D**2*RHO*U*CP*DT)
    rhs20=reaction(folder/"plus20",folder/"reaction20")*(COPIES*D**2*RHO*U*CP*DT)
    assert np.all(rhs[walls]<0) and np.all(rhs20[walls]<0)
    power=float(-rhs[walls].sum());change=float(abs(rhs20[walls].sum()/rhs[walls].sum()-1))
    energy=read(folder/"ledger/energy.json")
    assert change<.002 and abs(energy["relativeDifference"])<.005 and abs(energy["wallReleaseW"]-power)<1e-10
    assert energy["inputFrameSHA256"]==sha(folder/"restart.csv") and energy["sampledThetaExact"] and energy["maxWallNodeReactionDifferenceW"]<1e-11
    assert energy["observerSHA256"]==sha(folder/"ledger/thermistor_observe.exe")==source["observerSHA256"]
    weights=np.bincount(wall.ravel(),weights=np.repeat(area/3,3),minlength=len(xyz))
    q=np.mean(rhs[wall]/weights[wall],axis=1);p=P0+RHO*U**2*end_state[wall,4].mean(axis=1)
    t=T0+DT*end_state[near,9].mean(axis=1);assert abs(float(-area@q)-power)<1e-13
    ref=np.genfromtxt(str(BASE/"surface/wall.csv"),delimiter=",",names=True);groups=ref["regionID"].astype(int)
    body=read(folder/"body.json");assert np.array_equal(groups,body["regionID"]) and np.array_equal(ref["faceID"],np.arange(len(q)))
    np.testing.assert_array_equal(ref["areaM2"],area)
    errq=float(np.max(abs(q-.5*ref["heatFluxWm2"])));relq=float(np.max(abs(q/(.5*ref["heatFluxWm2"])-1)))
    errp=float(np.max(abs(p-ref["pressurePa"])));errt=float(np.max(abs(t-(T0+.5*(ref["nearWallCellTemperatureK"]-T0)))))
    inverse=T0+DT*power/base["totalOutwardW"]
    assert relq<.0005 and abs(power/(.5*base["totalOutwardW"])-1)<.0005
    assert errp<.001 and errt<.001 and abs(inverse-299.15)<.001
    regions=[]
    for j,label in enumerate(body["regionNames"],1):
        mask=groups==j
        regions.append(dict(name=label,faces=int(mask.sum()),halfNativeW=float(-area[mask]@q[mask]),
                       linearPredictionW=float(-area[mask]@(.5*ref["heatFluxWm2"][mask]))))
    out.mkdir(parents=True,exist_ok=False)
    np.savetxt(out/"direct_wall.csv",np.column_stack((ref["faceID"],groups,p,q,t)),delimiter=",",comments="",
               fmt=["%d","%d"]+["%.17g"]*3,header="faceID,regionID,pressurePa,heatFluxWm2,firstLayerK")
    save(out/"half.json",dict(faces=len(q),areaM2=float(area.sum()),regions=regions,wallReleaseW=power,scale=.5,referenceW=base["totalOutwardW"],
         solvedFromNativePowerK=inverse,maxAbsoluteFluxErrorWm2=errq,maxRelativeFluxError=relq,
         maxPressureErrorPa=errp,maxFirstLayerErrorK=errt,relativeHeatChange20=change,energy=energy,residualWindow=res,continuedWindow=res20,
         frames=frames,lastIteration=last,nativeSeconds=rec["elapsedSeconds"],continuationSeconds=plus["elapsedSeconds"],
         directory=folder.relative_to(ROOT).as_posix(),source=source,wallSHA256=sha(out/"direct_wall.csv")))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for path in folder.rglob("*.json"):
            if path.name!="body.json":z.write(path,path.relative_to(folder).as_posix())
        for name in ("case.cfg","history.csv","engine_stdout.log"):z.write(folder/name,name)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("folder",type=Path);parser.add_argument("output",type=Path)
    a=parser.parse_args();folder=a.folder.resolve();out=a.output.resolve()
    for path in (folder,out):path.relative_to(ROOT/"outputs")
    audit(folder,out)
