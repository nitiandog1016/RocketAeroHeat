"""读取两网原始初场的第一次方程装配，定位近壁细化引入的离散残差。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE,ROOT,RUN,LIMIT,OFFICIAL,ENGINE_SHA,BIN,inputs,configuration,replace,read,sha,observe,bound_source
from thermistorStartupSurface import columns,checked,write

REF=ROOT/"outputs/v292_validation/first01"
NAMES=["Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"]
RMS=["rms["+s+"]" for s in ("P","U","V","W","h","rho*Y_0")]


def probe(fine,out):
    decision=read(RUN/"first_state/state.json");assert not decision["passed"]
    src=RUN/"first01" if fine else REF
    checked(src);assert sha(OFFICIAL)==ENGINE_SHA
    if fine:
        assert sha(src/"seed_00000.csv")==inputs(RUN/"grid01")["seedSHA256"]
        assert (src/"case.cfg").read_text("ascii")==configuration(1)
    else:
        assert sha(src/"seed_00000.csv")==sha(bound_source("seed_00000.csv"))
        assert sha(src/"case.cfg")==read(HERE.parent/"v292/evidence/first_state.json")["configSHA256"]
    config=replace((src/"case.cfg").read_text("ascii"),f"INNER_ITER= {LIMIT if fine else 260}\n","INNER_ITER= 1\n")
    (out/"case.cfg").write_text(config,encoding="ascii")
    for name in ("mesh.su2","seed_00000.csv"):os.link(src/name,out/name)
    files=[HERE/"INITIAL_RHS.md",Path(__file__).resolve(),RUN/"first_state/state.json",OFFICIAL,src/"case.cfg",
           out/"case.cfg",out/"mesh.su2",out/"seed_00000.csv",ROOT/"scripts/runSU2NativeDebugger.py"]
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    observe([str(OFFICIAL),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in files},100,200*1024**2)
    checked(out)


def audit(out):
    arrays=[];summaries={};layers=[];volumes=[]
    for kind in ("coarse","fine"):
        small=RUN/("initial_"+kind);checked(small)
        fine=kind=="fine";src=RUN/"first01" if fine else REF
        history=(small/"history_00001.csv").read_bytes()
        assert len(history.splitlines())==2 and (src/"history_00001.csv").read_bytes().startswith(history)
        with np.load(RUN/"grid01/nodes.npz" if fine else bound_source("nodes.npz")) as z:
            xyz=z["xyz"];ids=z["ids"].reshape(49,-1,65)
            if fine:volume=z["dualVolume"]
        if not fine:
            v=columns(ROOT/"outputs/v292_validation/endpoint_observe/before_wall.csv",["id","volume"])
            v=v[np.argsort(v[:,0])];np.testing.assert_array_equal(v[:,0],np.arange(len(xyz)));volume=v[:,1]
        a=columns(small/"restart_00001.csv",["PointID","x","y","z"]+["Residual_"+s for s in NAMES])
        np.testing.assert_array_equal(a[:,0],np.arange(len(xyz)));np.testing.assert_allclose(a[:,1:4],xyz,rtol=1e-14,atol=2e-14)
        r=a[:,4:];assert np.isfinite(r).all()
        log=np.log10(np.sqrt(np.mean(r**2,axis=0)));h=columns(small/"history_00001.csv",RMS)
        error=float(np.max(abs(log-h[0])));assert error<1e-7
        layer=np.empty(len(xyz),int)
        for j in range(ids.shape[1]):layer[ids[:,j]]=j
        masks={"refinedBand":layer<=(16 if fine else 8),"wall":layer==0,"outerBeyondOldLayer10":layer>(18 if fine else 10)}
        if fine:masks["newMidpoints"]=np.arange(len(xyz))>=len(arrays[0])
        squared=r**2
        rate=np.sqrt(np.sum(squared/volume[:,None],axis=0)/volume.sum())
        summaries[kind]=dict(nodes=len(xyz),rmsLog=dict(zip(RMS,log)),nodeRMSMaxLogError=error,
            historyPrefixRows=1,volumeWeightedRMSLocalRate=dict(zip(NAMES,rate)),
            squaredResidualShares={name:dict(zip(NAMES,squared[mask].sum(axis=0)/squared.sum(axis=0))) for name,mask in masks.items()},
            nativeSeconds=read(small/"debugger_process.json")["elapsedSeconds"],frameSHA256=sha(small/"restart_00001.csv"))
        arrays.append(r);layers.append(layer);volumes.append(volume)
    a,b=arrays;shared=b[:len(a)];outside=layers[0]>10;inside=layers[0]<=8
    delta=shared-a;d2=delta**2;total=d2.sum(axis=0)
    common=dict(outerMaxRawDifference=dict(zip(NAMES,abs(delta[outside]).max(axis=0))),
        outerMaxLocalRateDifference=dict(zip(NAMES,abs(shared[outside]/volumes[1][:len(a)][outside,None]-a[outside]/volumes[0][outside,None]).max(axis=0))),
        refinedBandFractionOfSquaredSharedNodeDifference=dict(zip(NAMES,d2[inside].sum(axis=0)/total)))
    write(out/"initial_rhs.json",dict(probes=summaries,commonNodes=common,acceptedPhysicalSteps=0))


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("stage",choices=["coarse","fine","audit"]);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();out=RUN/("initial_"+a.stage)
    if a.worker:
        assert a.stage=="audit";audit(out)
    else:
        out.mkdir(exist_ok=False)
        if a.stage!="audit":probe(a.stage=="fine",out)
        else:
            files=[HERE/"INITIAL_RHS.md",Path(__file__).resolve(),RUN/"grid01/nodes.npz",bound_source("nodes.npz"),
                   ROOT/"outputs/v292_validation/endpoint_observe/before_wall.csv",RUN/"first_state/state.json"]
            files.extend(RUN/("initial_"+k)/n for k in ("coarse","fine") for n in ("debugger_process.json","restart_00001.csv","history_00001.csv"))
            rec=observe([sys.executable,str(Path(__file__).resolve()),"audit","--worker"],out,{str(p):sha(p) for p in files},180,200*1024**2)
            assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
