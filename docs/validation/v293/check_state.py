"""从加密网的全部实际节点和原生历史核对半步终点，并定位残差的壁法向分布。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE, ROOT, RUN, OFFICIAL, ENGINE_SHA, LIMIT, configuration, inputs, observe, read, sha, bound_source
from thermistorStartupSurface import columns,checked,write
sys.path[:0]=[str(HERE.parent/v) for v in ("v287","v288","v289")]
from audit_step import RMS, frame


def audit(index,grid,folder,out):
    info=inputs(grid);rec=checked(folder);end=folder/f"restart_{index:05}.csv"
    hist=folder/f"history_{index:05}.csv";seed=folder/f"seed_{index-1:05}.csv"
    assert (folder/"case.cfg").read_text("ascii")==configuration(index)
    assert rec["command"]==[str(OFFICIAL),"case.cfg","--threads","1"] and sha(OFFICIAL)==ENGINE_SHA
    assert sha(folder/"mesh.su2")==info["gridSHA256"] and end.name in (folder/"debugger.log").read_text(errors="replace")
    if index==1:expected=info["seedSHA256"]
    else:
        parent=read(folder.parent/"first_energy/accepted.json");assert parent["accepted"]
        expected=parent["energy"]["acceptedFrameSHA256"]
    assert sha(seed)==expected
    with np.load(grid/"nodes.npz") as z:xyz=z["xyz"];ids=z["ids"].reshape(49,57,65)
    assert len(xyz)==info["nodes"] and np.array_equal(ids[0],ids[-1])
    before=frame(seed,xyz);new=frame(end,xyz)
    assert np.any(before[:,4:]!=new[:,4:])
    h=columns(hist,["Time_Iter","Inner_Iter"]+RMS)
    assert np.isfinite(h).all() and 20<=len(h)<=LIMIT and np.all(h[:,0]==index)
    np.testing.assert_array_equal(h[:,1],np.arange(len(h)))
    rhs=columns(end,["Residual_"+s for s in ("Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0")])
    assert rhs.shape==(len(xyz),6) and np.isfinite(rhs).all()
    log=np.log10(np.sqrt(np.mean(rhs**2,axis=0)));rhs_error=float(max(abs(log-h[-1,2:])))
    layer=np.empty(len(xyz),int)
    for j in range(57):layer[ids[:,j]]=j
    eq=int(np.argmax(h[-20:,2:].max(axis=0)))
    residual=rhs[:,eq]**2;peak=int(np.argmax(residual))
    portions={k:float(residual[mask].sum()/residual.sum()) for k,mask in
              (("wall",layer==0),("first16FineIntervals",layer<=16),("newMidpoints",np.arange(len(xyz))>=info["oldNodes"]),
               ("axis",np.hypot(xyz[:,1],xyz[:,2])==0),("far",layer==56))}
    # 对原始方程右端作诊断，仍按全部原六分量作正式验收。
    pressure_rhs=rhs[:,0];enthalpy_rhs=rhs[:,4]
    coefficient=float((pressure_rhs@enthalpy_rhs)/(pressure_rhs@pressure_rhs))
    remainder=enthalpy_rhs-coefficient*pressure_rhs
    residual_alignment=dict(leastSquaresCoefficient=coefficient,
        cosine=float((pressure_rhs@enthalpy_rhs)/(np.linalg.norm(pressure_rhs)*np.linalg.norm(enthalpy_rhs))),
        unexplainedSquaredFraction=float((remainder@remainder)/(enthalpy_rhs@enthalpy_rhs)))
    angular=rhs[ids[:-1,:,1:-1],eq]
    modes=(abs(np.fft.rfft(angular,axis=0))**2).sum(axis=(1,2))/48
    modes[1:-1]*=2
    modes[0]+=residual[np.unique(ids[0,:,[0,-1]])].sum()
    assert abs(modes.sum()/residual.sum()-1)<1e-10
    worst_layers=np.bincount(layer,weights=residual)/residual.sum()
    g=read(bound_source("geometry.json"));p=g["airPa"]+g["rho"]*g["speedMps"]**2*new[:,4]
    t=new[:,9];w=np.unique(ids[:,0]);slip=float(np.max(abs(new[w,5:8])));wall_error=float(np.max(abs(t[w]-1)))
    window=h[-20:,2:].max(axis=0);valid=(h[:,2:]<-8).all(axis=1)
    where=np.flatnonzero(np.convolve(valid.astype(int),np.ones(20,int),"valid")==20)
    passed=bool(max(window)<-8 and rhs_error<1e-7 and min(p)>0 and min(t)>=-1e-12 and max(t)<=1+1e-12 and slip<=1e-12 and wall_error<=1e-12)
    write(out/"state.json",dict(step=index,passed=passed,innerRows=len(h),dtSeconds=g["dtSeconds"]/2,
        elapsedPhysicalTimeS=index*g["dtSeconds"]/2,firstCompleteWindowEndIndex=int(where[0]+19) if len(where) else None,
        trailingPassingRows=int(len(h)-1-np.flatnonzero(~valid)[-1]) if not valid.all() else len(h),
        firstLogRMS=dict(zip(RMS,h[0,2:])),last20MaxLogRMS=dict(zip(RMS,window)),finalLogRMS=dict(zip(RMS,h[-1,2:])),
        nodeRHSMaxLogError=rhs_error,limitingEquation=RMS[eq],squaredResidualShares=portions,
        enthalpyContinuityAlignment=residual_alignment,limitingResidualByLayer=worst_layers.tolist(),
        limitingResidualByAngularMode=(modes/residual.sum()).tolist(),
        peakResidualNode=peak,peakRadialLayer=int(layer[peak]),peakXYZOverD=xyz[peak].tolist(),
        pressureRangePa=[float(min(p)),float(max(p))],scalarRange=[float(min(t)),float(max(t))],maxWallScalarError=wall_error,
        maxWallVelocityOverU=slip,frameSHA256=sha(end),historySHA256=sha(hist),oldFrameSHA256=expected,
        gridSHA256=info["gridSHA256"],normalNativeTermination=True,nativeSeconds=rec["elapsedSeconds"]))


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("step",type=int,choices=[1,2]);cli.add_argument("--root",type=Path,default=RUN)
    cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);label="first" if a.step==1 else "second"
    grid=area/"grid01";src=area/(label+"01");out=area/(label+"_state")
    if a.worker:audit(a.step,grid,src,out)
    else:
        out.mkdir(exist_ok=False)
        files=[HERE/"PLAN.md",HERE/"run_step.py",Path(__file__).resolve(),ROOT/"scripts/thermistorStartupSurface.py",
               HERE.parent/"v287/audit_step.py",grid/"nodes.npz",grid/"grid.json",src/"case.cfg",src/"debugger_process.json",
               src/f"restart_{a.step:05}.csv",src/f"seed_{a.step-1:05}.csv",src/f"history_{a.step:05}.csv"]
        if a.step==2:files.append(area/"first_energy/accepted.json")
        rec=observe([sys.executable,str(Path(__file__).resolve()),str(a.step),"--root",str(area),"--worker"],out,
                    {str(p):sha(p) for p in files},180,200*1024**2)
        assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
