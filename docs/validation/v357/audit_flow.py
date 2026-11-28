"""逐文件核对固定低速风洞的动量网格、首步终态及四个原生残差历史。"""
import argparse
import csv
from pathlib import Path

from refine_flow import BASE,DOC,OLD,ROOT,SEED,read,save,sha
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from run_flow import EXE,HEADER

FIELDS=("rms[P]","rms[U]","rms[V]","rms[W]")


def checked_history(folder,prepared,n=1):
    p=read(prepared/"prepared.json")
    h=folder/f"history_{n:05d}.csv";frame=folder/f"restart_{n:05d}.csv"
    step=read(folder/"step.json");run=read(folder/"debugger_process.json")
    seed=folder/f"seed_{n-1:05d}.csv"
    assert run["status"]=="terminal" and run["inputSHA256"][str(EXE)]==sha(EXE)
    assert all(run["inputsStable"].values()) and not (run["survivingProcessIDs"] or run["supervisorError"])
    assert sha(folder/"case.cfg")==run["inputSHA256"][str(folder/"case.cfg")]
    assert sha(seed)==step["seedSHA256"]==run["inputSHA256"][str(seed)]
    assert p["files"]["mesh.su2"]==sha(folder/"mesh.su2")==step["meshSHA256"]
    common=dict(branch=step["branch"],stepSize=step["stepSize"],elapsedPhysicalS=step["elapsedPhysicalS"],
                meshSHA256=step["meshSHA256"],seedSHA256=step["seedSHA256"],cfgSHA256=sha(folder/"case.cfg"),
                status=run["status"],exitCode=run["debuggerExitCode"],stopReason=run["stopReason"],
                elapsedSeconds=run["elapsedSeconds"],peakRSSBytes=run["peakFamilyRSSBytes"],
                peakOutputBytes=run["peakOutputBytes"],inputsStable=True,survivingProcesses=run["survivingProcessIDs"])
    if run["stopReason"] or run["debuggerExitCode"]!=0:
        assert not step["accepted"]
        return dict(**common,accepted=False,normalNativeCompletion=False),None
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")
    with h.open() as stream:
        reader=csv.DictReader(stream,skipinitialspace=True);reader.fieldnames=[k.strip() for k in reader.fieldnames]
        rows=list(reader)
    assert len(rows)>=20 and all(int(r["Time_Iter"])==n and int(r["Inner_Iter"])==i for i,r in enumerate(rows))
    arr=np.array([[float(r[key]) for key in FIELDS] for r in rows])
    limit=arr[-20:].max(axis=0)
    np.testing.assert_array_equal(limit,[step["last20Maximum"][k] for k in FIELDS])
    assert sha(frame)==step["frameSHA256"] and sha(h)==step["historySHA256"]
    with frame.open() as stream:assert next(csv.reader(stream))==HEADER
    raw=np.loadtxt(frame,delimiter=",",skiprows=1)
    assert raw.shape==(p["nodes"],8) and np.isfinite(raw).all()
    np.testing.assert_array_equal(raw[:,0],np.arange(p["nodes"]))
    with np.load(prepared/"stations.npz") as st:
        np.testing.assert_allclose(raw[:,1:4],st["xyz"],atol=1e-13,rtol=0)
        vmax=float(abs(raw[st["wall"],5:8]).max())
    pmin=float((101325+1.177*.05**2*raw[:,4]).min())
    assert pmin==step["minimumAbsolutePressurePa"] and vmax==step["maxWallSpeedOverU"]
    passed=bool(np.isfinite(arr).all() and (limit < -8).all() and pmin>0 and vmax<1e-9)
    assert passed==step["accepted"]
    return dict(**common,accepted=passed,normalNativeCompletion=True,rows=len(arr),
                last20Maximum=dict(zip(FIELDS,limit)),lastSample=dict(zip(FIELDS,arr[-1])),
                last100Drop=dict(zip(FIELDS,arr[-100]-arr[-1])) if len(arr)>=100 else {},
                frameSHA256=sha(frame),historySHA256=sha(h),minimumAbsolutePressurePa=pmin,maxWallSpeedOverU=vmax),arr


def audit(work,output):
    candidate=work/"prepared";g=read(candidate/"refinement.json");model=read(candidate/"prepared.json")
    assert sha(candidate/"meridian.npz")==g["meridianSHA256"] and all(sha(candidate/k)==v for k,v in g["newFiles"].items())
    for path,key in ((SEED/"initial/fields.npz","savedInitial"),(SEED/"steady/temperature.npz","heatMesh")):
        assert sha(path)==g[key]
    old=read(OLD/"prepared.json")
    assert g["oldNodesChecked"]==old["nodes"] and model["originalFaces"]==old["originalFaces"]==9216
    with np.load(candidate/"meridian.npz") as layout,np.load(OLD/"stations.npz") as original,np.load(candidate/"stations.npz") as new:
        to_new=layout["oldToNew"]
        np.testing.assert_allclose(new["xyz"][to_new],original["xyz"],rtol=0,atol=1e-13)
        assert len(np.unique(to_new))==old["nodes"]
    prior,pdata=checked_history(ROOT/"outputs/v354_validation/flow01/disturbed_half_1",OLD)
    current,cdata=checked_history(work/"flow/disturbed_half_1",candidate)
    assert prior["accepted"] and prior["cfgSHA256"]==current["cfgSHA256"]
    assert prior["elapsedPhysicalS"]==current["elapsedPhysicalS"]
    output.mkdir(exist_ok=False,parents=True)
    save(output/"flow_check.json",dict(mesh=g,originalAreaM2=model["originalAreaM2"],wallAreaM2=model["wallAreaM2"],
         volumeD3=model["volumeD3"],minimumVolumeD3=model["minCellVolumeD3"],farAndWallMarkers=model["markers"],
         patches=model["checkedPatches"],coarse=prior,refined=current))
    if cdata is not None:
        plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False})
        fig,axes=plt.subplots(1,2,figsize=(9.0,3.9),sharey=True,layout="constrained")
        for ax,values,title in zip(axes,(pdata,cdata),("原动量网 · 190,106 节点","加密动量网 · 547,066 节点")):
            for col,label,color,ls in zip(range(4),("压力","X 速度","Y 速度","Z 速度"),
                    ("#182f49","#bb502b","#178476","#6c57a2"),("-","-","-","--")):
                ax.plot(values[:,col],lw=1.4,ls=ls,color=color,label=label)
            ax.axhline(-8,ls=":",color="#444");ax.axvspan(len(values)-20,len(values)-1,color="#7e91a6",alpha=.15)
            ax.set(xlabel="内迭代次数",title=title);ax.spines[["top","right"]].set_visible(False)
        axes[0].set_ylabel(r"SU2 原生 $\log_{10}$ RMS");axes[0].legend(frameon=False,fontsize=8,ncol=2)
        fig.suptitle("0.05 m/s 风洞首步 · 目标物理时刻 0.300 s",fontsize=13)
        fig.savefig(output/"first_step.png",dpi=155);plt.close(fig)
    print("first step accepted:",current["accepted"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();work=args.work.resolve();work.relative_to(BASE);audit(work,args.output)
