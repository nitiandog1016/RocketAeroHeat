"""从冻结的二进制 STL 与 MATLAB 返回的状态掩码复算五个飞行工况的面积。"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]


def check(area,p,q,summary):
    counts={"computedPressure":p==1,"prescribedPressure":p==2,"heat":q,"joint":(p==1)&q}
    assert len(area)==summary["faces"] and abs(area.sum()-summary["areaM2"])<1e-12
    for label,mask in counts.items():
        assert int(mask.sum())==summary[label+"Faces"]
        assert abs(area[mask].sum()-summary[label+"AreaM2"])<1e-12


def verify(source):
    case=ROOT/"cases/simple_launcher"
    definition=json.loads((case/"geometry.json").read_bytes())
    stl=case/definition["file"];actual=hashlib.sha256(stl.read_bytes()).hexdigest()
    record=json.loads((source/"cases.json").read_text(encoding="utf-8-sig"))
    assert actual==definition["fileSHA256"]==record["surfaceSHA256"]
    dtype=np.dtype([("normal","<f4",3),("vertices","<f4",(3,3)),("attribute","<u2")])
    tri=np.fromfile(stl,dtype=dtype,offset=84)["vertices"].astype(float)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert len(tri)==definition["faces"]==record["faces"] and abs(area.sum()-definition["areaM2"])<1e-12
    cases=record["cases"];assert [r["mach"] for r in cases]==[.30,.85,1.,1.20,2.]
    bundle=source/"status_masks.npz"
    packed={}
    for i,r in enumerate(cases,1):
        key=f"case_{i:02d}"
        if bundle.is_file():
            with np.load(bundle) as masks:
                p=masks[key+"_pressure"];q=masks[key+"_heat"]
        else:
            raw=loadmat(source/f"masks_{i:02d}.mat")
            p,q=raw["p"].ravel(),raw["q"].ravel()
        assert p.shape==q.shape==(len(area),) and np.isin(p,[0,1,2]).all() and np.isin(q,[0,1]).all()
        p=p.astype(np.uint8);q=q.astype(bool)
        assert r["altitudeM"]==5000 and r["wallGivenK"]==300 and abs(r["speedMps"]/r["soundSpeedMps"]-r["mach"])<1e-12
        check(area,p,q,r["coverage"])
        for g in definition["regions"]:
            ids=slice(g["firstFace"]-1,g["lastFace"])
            check(area[ids],p[ids],q[ids],r["regions"][g["name"]])
        assert r["fixedBoundaryError"]=="RAH:Launcher:Boundary"
        if not r["solverReturned"]:
            assert r["engineeringError"]=="RocketAeroHeat:SingleCase:Mach" and r["lastStage"]=="freestream"
            assert not p.any() and not q.any()
        if r["solverReturned"] and not np.all((p==1)&q):
            assert r["publicationError"]=="RocketAeroHeat:WorkflowIncompleteSurface"
        assert not r["resultPageAccepted"] and not r["wholeSurfacePowerKnown"]
        packed[key+"_pressure"]=p;packed[key+"_heat"]=q
    return record,packed


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path);parser.add_argument("--save-to",type=Path)
    args=parser.parse_args();record,masks=verify(args.source)
    if args.save_to:
        args.save_to.mkdir(exist_ok=False,parents=True)
        np.savez_compressed(args.save_to/"status_masks.npz",**masks)
        (args.save_to/"cases.json").write_bytes((json.dumps(record,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode())
    print(len(record["cases"]),"cases,",record["faces"],"STL faces each; all coverage checks passed")
