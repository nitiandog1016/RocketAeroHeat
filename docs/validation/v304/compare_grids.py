"""按相同 10,272 个面比较三次独立热敏珠计算，产出工作台的小型原面包。"""
import argparse
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np

from audit_third_mesh import FINE, HERE, MEDIUM, ROOT, WIRED, read, sha, store, terminal
from check_frames import inspect


def load_fields(folder):
    terminal(folder)
    s=read(folder/"surface/summary.json")
    assert s["sourceFrameSHA256"]==sha(folder/"restart.csv") and s["sourceReactionSHA256"]==sha(folder/"reaction0/restart.csv")
    assert s["wallSHA256"]==sha(folder/"surface/wall.csv") and s["meshSHA256"]==sha(folder/"mesh.su2")
    assert max(s["residualWindow"].values()) < -8 and max(s["plus20Window"].values()) < -8
    assert s["energy"]["sampledThetaExact"] and abs(s["energy"]["relativeDifference"]) < .005 and s["relativeChange20"] < .002
    w=np.genfromtxt(folder/"surface/wall.csv",names=True,delimiter=",")
    if folder==WIRED:
        flow=gzip.decompress((HERE.parent/"v297/evidence/nearwall_faces.csv.gz").read_bytes())
        assert hashlib.sha256(flow).hexdigest()=="2d92881abecb315b0ac2f10a0cf92f231ec86a07c3716db8324649e6bcc43c4e"
    else:
        f=read(folder/"surface/flow.json")
        assert f["sourceFrameSHA256"]==s["sourceFrameSHA256"] and f["tableSHA256"]==sha(folder/"surface/flow.csv")
        flow=(folder/"surface/flow.csv").read_bytes()
    reverse=np.genfromtxt(io.BytesIO(flow),names=True,delimiter=",")
    np.testing.assert_array_equal(w["faceID"],np.arange(10272))
    np.testing.assert_array_equal(reverse["face_id"],w["faceID"]+1)
    assert abs(-w["areaM2"]@w["heatFluxWm2"]-s["totalOutwardW"])<1e-13
    return s,w,reverse["reverse_both"].astype(bool)


def compare(out):
    checks=read(FINE.parent/"checkpoints.json")
    assert checks["fieldChecksPassed"] and max(checks["residualWindow"].values()) < -8
    for frame in checks["frames"]:
        assert frame["sha256"]==sha(FINE/frame["file"])
    with np.load(FINE/"nodes.npz") as nodes:
        plus=inspect(FINE/"plus20/restart.csv",nodes["xyz"],np.unique(nodes["wallMap"]))
    assert plus["outside"]==0 and plus["minPressurePa"]>0 and plus["maxWallSpeedMps"]<1e-10 and plus["maxWallThetaError"]==0
    folders=(WIRED,MEDIUM,FINE)
    loaded=[load_fields(f) for f in folders]
    base=loaded[0][1];area=base["areaM2"]
    assert loaded[0][0]==read(HERE.parent/"v297/evidence/fine_surface.json")
    assert loaded[1][0]==read(HERE.parent/"v300/evidence/summary.json")
    for s,w,_ in loaded:
        np.testing.assert_array_equal(w["regionID"],base["regionID"])
        np.testing.assert_array_equal(w["areaM2"],area)
        assert s["bodySHA256"]==sha(FINE/"body.json") and s["caseSHA256"]==sha(FINE/"case.json")
    out.mkdir(parents=True,exist_ok=False)
    s3,_,reverse=loaded[-1];body=read(FINE/"body.json")
    info=dict(areaM2=float(area.sum()),faces=len(area),newFaces=0,newAreaM2=0,jointExperimentVerifiedM2=0,plus20Frame=plus,
              grids=[dict(meshSHA256=s["meshSHA256"],frameSHA256=s["sourceFrameSHA256"],reactionSHA256=s["sourceReactionSHA256"],wallSHA256=s["wallSHA256"],
                          conductanceWK=s["totalOutwardW"]/12,beadAndJointsW=s["beadAndJointsW"],meanHWm2K=s["beadMeanHWm2K"],
                          relativeExperimentDifference=s["beadMeanHWm2K"]/s["experimentMeanHWm2K"]-1) for s,_,_ in loaded])
    changed=loaded[1][2]!=reverse
    info["reverseChangedFaces"]=int(changed.sum());info["reverseChangedAreaM2"]=float(area[changed].sum())
    with zipfile.ZipFile(out/"face_comparison.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for mode in ("givenWall","givenPower"):
            data=[base["faceID"]+1,base["regionID"],area]
            keys=["face_id","region_id","area_m2"];fields=[];totals=[];wall_k=[];peaks=[]
            for j,(s,w,_) in enumerate(loaded,1):
                k=1.0 if mode=="givenWall" else .004/s["totalOutwardW"]
                q=k*w["heatFluxWm2"];t=293.15+k*(w["nearWallCellTemperatureK"]-293.15)
                tw=293.15+12*k;fields.append(q);wall_k.append(tw);totals.append(float(-area@q));peaks.append(int(np.argmin(q))+1)
                if mode=="givenPower":assert abs(totals[-1]-.004)<1e-13
                keys.extend((f"p{j}_Pa",f"q{j}_into_W_m2",f"first{j}_K",f"Tw{j}_K"))
                data.extend((w["pressurePa"],q,t,np.full(len(q),tw)))
            regions=[]
            for i,name in enumerate(body["regionNames"],1):
                mask=base["regionID"]==i
                regions.append(dict(name=name,faces=int(mask.sum()),areaM2=float(area[mask].sum()),
                    outwardW=[float(-area[mask]@q[mask]) for q in fields],
                    peakFace=[int(np.flatnonzero(mask)[np.argmin(q[mask])])+1 for q in fields]))
            np.testing.assert_allclose(np.array([r["outwardW"] for r in regions]).sum(axis=0),totals,rtol=0,atol=1e-13)
            changes=[]
            for first,second in zip(fields,fields[1:]):
                d=second-first;j=int(np.argmax(abs(d)))
                changes.append(dict(maxAbsWm2=float(abs(d[j])),signedWm2=float(d[j]),faceID=j+1,rmsWm2=float(np.sqrt(np.average(d*d,weights=area)))))
            info[mode]=dict(wallK=wall_k,totalOutwardW=totals,peakFace=peaks,peakOutWm2=[float(-q.min()) for q in fields],regions=regions,changes=changes)
            buf=io.StringIO();np.savetxt(buf,np.column_stack(data),fmt=["%d","%d"]+["%.17g"]*(len(keys)-2),delimiter=",",header=",".join(keys),comments="")
            z.writestr(mode+".csv",buf.getvalue())
        z.writestr("sources.json",json.dumps(info,ensure_ascii=False).encode("utf8"))
    store(out/"comparison.json",info)
    archive=out/"leadsThird3.zip"
    with zipfile.ZipFile(archive,"x") as z:
        for p in (FINE/"body.json",FINE/"case.json",*(FINE/"surface"/n for n in ("summary.json","wall.csv","flow.csv","flow.json"))):
            z.writestr(zipfile.ZipInfo(p.name,date_time=(2026,1,1,0,0,0)),p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED)
    g=read(FINE/"geometry.json")
    contract={k:g[k] for k in ("nodes","cells","meshSHA256","bodySHA256","nodeMapSHA256")}
    contract.update(initialSHA256=sha(FINE/"initial_flow.csv"),configSHA256=sha(FINE/"case.cfg"),caseSHA256=s3["caseSHA256"],archiveSHA256=sha(archive))
    store(out/"leadsThird3.json",contract)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("output",type=Path)
    target=parser.parse_args().output.resolve();target.relative_to(ROOT/"outputs");compare(target)
