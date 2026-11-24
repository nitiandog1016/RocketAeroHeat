"""保留真实运行与边界小表；体积储热分块公开，完整原帧及逐边量登记本机指纹。"""
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np
from build_time import ROOT,HERE,sha,read


def pack(base,work,out):
    small={};large={}
    def copy(path):
        small[path.relative_to(base).as_posix()]=path.read_bytes()
    for path in (base/"build01/CSpeciesSolver_process.json",base/"build01/CSpeciesSolver.log",
                 *list((base/"build02").glob("*.json")),work/"job.json",base/"audit.json",base/"post_process.json",base/"audit_process.json"):
        copy(path)
    for branch,step in (("full",1),("half",1),("half",2)):
        d=work/branch/f"step{step}";e=read(d/"energy.json")
        for folder in (d,d/"release"):
            for path in folder.iterdir():
                if path.is_dir():continue
                if path.name in ("case.cfg","debugger_process.json","debugger.log","native.json","observed.json","energy.json",
                                 "stages.csv","far_pairs.csv","wall_budget.npz","physical_step_initial.py") or path.name.startswith("history_"):
                    copy(path)
                if path.suffix in (".bin",".vtk") or path.name.startswith(("seed_","restart_")):
                    large[path.relative_to(base).as_posix()]=sha(path)
        old=np.loadtxt(d/f"seed_{step-1:05d}.csv",skiprows=1,delimiter=",",usecols=9)
        new=np.loadtxt(d/f"restart_{step:05d}.csv",skiprows=1,delimiter=",",usecols=9)
        n=np.fromfile(d/"nodes.bin","<f8").reshape(-1,18);order=np.argsort(n[:,0]);n=n[order]
        pre=np.fromfile(d/"before_time.bin","<f8")[order];post=np.fromfile(d/"after_time.bin","<f8")[order]
        start=np.arange(0,len(n),4096)
        rows=np.c_[start,np.minimum(start+4096,len(n)),np.add.reduceat(n[:,5]*(new-old)*e["scale"]/e["dt"],start),
                   np.add.reduceat((post-pre)*e["scale"],start)]
        tmp=io.StringIO();np.savetxt(tmp,rows,delimiter=",",fmt=["%d","%d","%.17g","%.17g"],
                                    header="firstID,endExclusive,storageFromFramesW,storageFromNativeW",comments="")
        small[f"storage/{branch}{step}.csv"]=tmp.getvalue().encode("ascii")
    for branch in ("full","half"):
        d=base/"ui01"/branch
        for name in ("ui.json","debugger_process.json"):copy(d/name)
    copy(base/"ui01/export_verified.json")
    small["manifest.json"]=json.dumps(dict(files={k:hashlib.sha256(v).hexdigest() for k,v in small.items()},largeLocal=large),indent=2).encode()
    with ZipFile(out,"x",compression=ZIP_DEFLATED,compresslevel=9) as z:
        for k,v in small.items():z.writestr(k,v)


def check(path,*,prefix="run01",nodes=261708):
    with ZipFile(path) as z:
        manifest=json.loads(z.read("manifest.json"))
        for k,h in manifest["files"].items():assert hashlib.sha256(z.read(k)).hexdigest()==h
        for branch,step in (("full",1),("half",1),("half",2)):
            root=f"{prefix}/{branch}/step{step}/";e=json.loads(z.read(root+"energy.json"))
            for sub in ("","release/"):
                rec=json.loads(z.read(root+sub+"debugger_process.json"))
                assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
                assert not(rec["survivingProcessIDs"] or rec["supervisorError"] or rec["stopReason"])
            with np.load(io.BytesIO(z.read(root+"wall_budget.npz"))) as v:
                assert len(v["wet"])==11906 and len(v["farID"])==2690
                heat=(-v["nodeInW"]).sum();far=(np.maximum(v["farMass"],0)@v["farTheta"])*e["scale"]
                np.testing.assert_allclose(v["nodeInW"],-(v["released"]+v["edgeCorrection"])*e["scale"],rtol=1e-9,atol=1e-10)
            s=np.genfromtxt(z.open(f"storage/{branch}{step}.csv"),names=True,delimiter=",")
            np.testing.assert_array_equal(s["firstID"][1:],s["endExclusive"][:-1]);assert s["firstID"][0]==0 and s["endExclusive"][-1]==nodes
            assert abs(s["storageFromFramesW"].sum()-e["storageW"])<1e-9 and abs(s["storageFromNativeW"].sum()-e["storageW"])<1e-9
            assert abs(heat-e["wallW"])<1e-10 and abs(far-e["farW"])<1e-10
            assert abs((heat-far-s["storageFromFramesW"].sum())/heat)<.005


if __name__=="__main__":
    base=ROOT/"outputs/v333_validation";archive=HERE/"evidence/receipts.zip"
    if "--check" not in sys.argv:pack(base,base/"run01",archive)
    check(archive);print("Verified three physical steps, separate release reactions and independently summed storage.")
