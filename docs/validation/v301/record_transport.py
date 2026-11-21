"""归档实际原生终点、温界失败的节点和与一阶基线的同面流场证据。"""
import argparse
import csv
import gzip
import io
import json
import zipfile
from pathlib import Path

import numpy as np

from check_frames import check
from transport_case import BASE, HERE, ROOT, read, sha, terminal
from thermistorBead import D, DT, RHO, U


def save(path,data):
    with path.open("x",encoding="utf8") as f:
        json.dump(data,f,indent=2,ensure_ascii=False,allow_nan=False)


def boundary_nodes(mesh):
    """从实际 SU2 边界面读取节点编号，判别低温点是否恰在远场或对称边界。"""
    result={}
    with mesh.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith("NMARK="):break
        else:raise ValueError("缺少体网边界记录")
        for _ in range(int(line.split("=")[1])):
            name=next(stream).split("=")[1].strip()
            count=int(next(stream).split("=")[1]);nodes=set()
            for _ in range(count):
                kind,*corners=map(int,next(stream).split())
                assert len(corners)=={5:3,9:4}[kind]
                nodes.update(corners)
            result[name]=nodes
    return result


def record(folder,out):
    receipt=terminal(folder)
    assert not receipt.get("supervisorError") and receipt["exitCode"]==0
    report=check(folder);frames=report["frames"];assert frames
    s=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1)
    base=np.loadtxt(BASE/"restart.csv",delimiter=",",skiprows=1)
    assert np.array_equal(s[:,:4],base[:,:4])
    current=frames[-1]
    state_ok=all(f["minAbsolutePressurePa"]>0 and f["maxWallScalarError"]==0 and f["wallSpeedMaxMps"]<1e-10*U for f in frames)
    eligible=report["checkpointsWithinBounds"] and report["residualPass"] and state_ok
    out.mkdir(parents=True,exist_ok=False)
    ends=boundary_nodes(folder/"mesh.su2");bad=np.flatnonzero((s[:,9]<-1e-12)|(s[:,9]>1+1e-12))
    xyz_m=s[:,1:4]*D
    counts={name:sum(int(i) in ids for i in bad) for name,ids in ends.items()}
    rows=io.StringIO(newline="");writer=csv.writer(rows)
    writer.writerow(("pointID","x_m","y_m","z_m","candidate_theta","base_terminal_theta","at_far","at_wall","at_symmetry"))
    symmetry=ends["sym0"]|ends["sym1"]
    for i in bad:
        writer.writerow((i,*xyz_m[i],s[i,9],base[i,9],int(i in ends["far"]),int(i in ends["wall"]),int(i in symmetry)))
    with (out/"out_of_bounds.csv.gz").open("xb") as f:
        with gzip.GzipFile(filename="",mode="wb",fileobj=f,mtime=0) as compressed:
            compressed.write(rows.getvalue().encode("ascii"))
    same_nodes=True
    for f in frames:
        theta=np.loadtxt(folder/f["file"],delimiter=",",skiprows=1,usecols=(9,))
        same_nodes &= bool(np.array_equal(np.flatnonzero((theta < -1e-12)|(theta > 1+1e-12)),bad))
    delta_t=(s[:,9]-base[:,9])*DT;peak=int(np.argmax(abs(delta_t)))
    with np.load(folder/"nodes.npz") as mapping:
        first=np.unique(mapping["firstMap"])
    report.update(acceptedForReactionAudit=bool(eligible),terminalBoundsPass=current["above"]+current["below"]==0,
        pressureAndWallChecksPass=state_ok,baseTerminalIteration=read(HERE.parent/"v300/evidence/frames.json")["iteration"],
        terminalToBaseTerminal=dict(maxPressureDeltaPa=float(RHO*U**2*np.max(abs(s[:,4]-base[:,4]))),
        maxVelocityVectorDeltaMps=float(U*np.linalg.norm(s[:,5:8]-base[:,5:8],axis=1).max()),
        maxFluidTemperatureDeltaK=float(abs(delta_t[peak])),deltaAtPeakK=float(delta_t[peak]),peakNodeID=peak,
        peakXYZM=xyz_m[peak].tolist(),maxAbsFirstLayerTemperatureDeltaK=float(np.max(abs(delta_t[first]))),
        sourceFrameSHA256=sha(BASE/"restart.csv")),
        sameViolatingNodeSetInAllSavedFrames=same_nodes,
        outsideNodesAtBoundary=counts,outsideNodesAwayFromAllMarkers=sum(not any(i in ids for ids in ends.values()) for i in bad),
        fractionOutside=len(bad)/len(s),outsideTableSHA256=sha(out/"out_of_bounds.csv.gz"),
        seconds=receipt["elapsedSeconds"],startedUTC=receipt["startedUTC"],endedUTC=receipt["endedUTC"],
        processID=receipt["processID"],resources=read(folder/"resources.json"))
    save(out/"frames.json",report)
    save(out/"inputs.json",read(folder/"inputs.json"))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for name in ("case.cfg","engine_request.json","engine_process.json","resources.json","history.csv","engine_stdout.log","supervisor.log"):
            z.write(folder/name,name)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("folder",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();folder=a.folder.resolve();out=a.output.resolve()
    for path in (folder,out):path.relative_to(ROOT/"outputs")
    report=record(folder,out)
    print(json.dumps({k:report[k] for k in ("lastHistoryIteration","residualPass","terminalBoundsPass","acceptedForReactionAudit","outsideNodesAtBoundary")}))
