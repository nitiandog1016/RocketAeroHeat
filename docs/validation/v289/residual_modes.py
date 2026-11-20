"""按完整周向网格分解标量残差；轴端节点只计一次，保留节点级原生平方和。"""
import argparse
from pathlib import Path
import numpy as np
from check_fields import columns
from run_relaxation import ASSETS, FULL, LAST, frozen, read, save, sha


def spectrum(rhs,ids,xyz):
    axis_mask=np.hypot(xyz[:,1],xyz[:,2])==0
    ring_ids=ids[:-1,~axis_mask[ids[0]]]
    np.testing.assert_array_equal(np.sort(ring_ids.ravel()),np.flatnonzero(~axis_mask))
    rays=len(ring_ids)
    harmonics=np.fft.rfft(rhs[ring_ids],axis=0)
    energy=np.sum(np.abs(harmonics)**2,axis=1)/rays
    energy[1:-1]*=2
    norm2=np.sum(rhs**2)
    axis_energy=np.sum(rhs[axis_mask]**2)
    energy[0]+=axis_energy
    assert np.isclose(energy.sum(),norm2,rtol=1e-11,atol=0) and norm2>0
    j=int(np.argmax(np.abs(rhs)))
    return dict(nodeRMS=float(np.sqrt(norm2/len(rhs))),m0=float(energy[0]/norm2),
                m1=float(energy[1]/norm2),higher=float(energy[2:].sum()/norm2),
                axisFraction=float(axis_energy/norm2),peakNode=j,peakXYZoverD=xyz[j].tolist())


def compare(run,out):
    forced=frozen(LAST/"probe_mirror01/restart_00001.csv",ASSETS["probe_mirror01/restart_00001.csv"]["sha256"])
    calm=frozen(run/"probe01/restart_00001.csv",read(run/"control01/control.json")["frameSHA256"])
    names=["Residual_Pressure","Residual_Enthalpy"]
    a,b=columns(forced,names),columns(calm,names)
    with np.load(FULL/"nodes.npz") as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    result={}
    for j,key in enumerate(("mass","enthalpy")):
        result[key]={name:spectrum(v,ids,xyz) for name,v in
                     (("original",a[:,j]),("noImpulse",b[:,j]),("difference",a[:,j]-b[:,j]))}
        result[key]["differenceOverOriginalNorm"]=float(np.linalg.norm(a[:,j]-b[:,j])/np.linalg.norm(a[:,j]))
    if (run/"step01/restart_00001.csv").is_file():
        audit=read(run/"finish01/first_step.json")
        endpoint=frozen(run/"step01/restart_00001.csv",audit["frameSHA256"])
        final=columns(endpoint,names)
        for j,key in enumerate(("mass","enthalpy")):
            these=spectrum(final[:,j],ids,xyz)
            native=audit["lastLogRMS"][["rms[P]","rms[h]"][j]]
            assert abs(np.log10(these["nodeRMS"])-native)<1e-7
            result[key]["physicalStepFinalAssembly"]=these
    if (run/"window_audit01/window.json").is_file():
        window=read(run/"window_audit01/window.json")
        endpoint=frozen(run/"window01/restart_00001.csv",window["sha256"]["restart_00001.csv"])
        final=columns(endpoint,names)
        for j,key in enumerate(("mass","enthalpy")):
            these=spectrum(final[:,j],ids,xyz)
            assert abs(np.log10(these["nodeRMS"])-window["finalLogRMS"][["rms[P]","rms[h]"][j]])<1e-7
            result[key]["windowCheckFinalAssembly"]=these
    save(out,result)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("run",type=Path)
    cli.add_argument("output",type=Path)
    args=cli.parse_args()
    compare(args.run,args.output)
