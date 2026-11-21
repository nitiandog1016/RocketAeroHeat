"""将两种网距和两种热敏珠外形配对，分别积分七区、共同表面与几何替换区。"""
import argparse
import gzip
import hashlib
import zipfile
from pathlib import Path

import numpy as np

from mesh_audit import HERE, ROOT, WIRED, NEW, SMOOTH, REF_SMOOTH, read, save, sha, shared_faces
from compare_meshes import audited_field, table
from thermistorBead import AREA, DT


def integral(w, ids):
    return float(-w["areaM2"][ids] @ w["heatFluxWm2"][ids])


def run(out):
    inputs=read(HERE/"evidence/inputs.json")
    assert all(sha(NEW/k)==v for k,v in inputs["inputs"].items())
    prior=read(HERE.parent/"v297/evidence/fine_surface.json")
    assert prior==read(WIRED/"surface/summary.json") and prior["sourceFrameSHA256"]==sha(WIRED/"restart.csv")
    assert prior["wallSHA256"]==sha(WIRED/"surface/wall.csv")
    old=table((WIRED/"surface/wall.csv").read_bytes())
    raw=gzip.decompress((HERE.parent/"v297/evidence/nearwall_faces.csv.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest()=="2d92881abecb315b0ac2f10a0cf92f231ec86a07c3716db8324649e6bcc43c4e"
    old_flow=table(raw);old_r=old_flow["reverse_both"].astype(bool)
    summary,w,flow=audited_field(NEW);r=flow["reverse_both"].astype(bool)
    smooth0,_,_=audited_field(SMOOTH)
    smooth,s,sf=audited_field(REF_SMOOTH)
    assert smooth0==read(HERE.parent/"v298/evidence/surface.json") and smooth==read(HERE.parent/"v299/evidence/surface.json")
    assert np.array_equal(w["faceID"]+1,old_flow["face_id"]) and np.array_equal(old["faceID"],w["faceID"])
    assert np.array_equal(old["areaM2"],w["areaM2"]) and np.array_equal(w["regionID"],old["regionID"])
    body=read(NEW/"body.json");wi,si=shared_faces(body,read(REF_SMOOTH/"body.json")).T
    np.testing.assert_array_equal(w["areaM2"][wi],s["areaM2"][si])
    coeff=summary["experimentMeanHWm2K"]
    four={}
    for label,record in (("baseSmooth",smooth0),("baseWired",prior),("fineSmooth",smooth),("fineWired",summary)):
        h=record["beadAndJointsW"]/(AREA*DT)
        assert abs(h-record["beadMeanHWm2K"])<1e-10
        four[label]=dict(beadAndJointsW=record["beadAndJointsW"],totalOutwardW=record["totalOutwardW"],
                         meanHWm2K=h,relativeExperimentError=h/coeff-1,sourceFrameSHA256=record["sourceFrameSHA256"])
    x=np.asarray(body["V"])[np.asarray(body["F"])-1,0]
    area=w["areaM2"];regions=[]
    for i,name in enumerate(body["regionNames"],1):
        use=w["regionID"]==i
        peak=int(np.flatnonzero(use)[np.argmin(w["heatFluxWm2"][use])])
        rev=use&r
        regions.append(dict(name=name,faces=int(use.sum()),areaM2=float(area[use].sum()),
              baseW=integral(old,use),fineW=integral(w,use),deltaW=integral(w,use)-integral(old,use),
              reverseFaces=int(rev.sum()),reverseAreaM2=float(area[rev].sum()),
              reverseXRangeM=[float(x[rev].min()),float(x[rev].max())] if rev.any() else [],
              reverseChangedFaces=int(np.count_nonzero(use&(r!=old_r))),peakFaceID=peak+1,
              peakXRangeM=[float(x[peak].min()),float(x[peak].max())],peakIntoWallWm2=float(w["heatFluxWm2"][peak])))
    bead_delta=summary["beadAndJointsW"]-prior["beadAndJointsW"]
    geom_delta=summary["beadAndJointsW"]-smooth["totalOutwardW"]
    original_geom_delta=prior["beadAndJointsW"]-smooth0["totalOutwardW"]
    common_delta=integral(w,wi)-integral(s,si)
    replace_delta=integral(w,np.isin(w["regionID"],[3,5]))-integral(s,s["regionID"]!=2)
    assert abs(geom_delta-common_delta-replace_delta)<1e-13
    assert abs(sum(i["deltaW"] for i in regions[2:5])-bead_delta)<1e-13
    assert abs(sum(i["fineW"] for i in regions)-summary["totalOutwardW"])<1e-13
    assert abs(original_geom_delta-read(HERE.parent/"v298/evidence/comparison.json")["geometryDeltaW"])<1e-13
    result=dict(cases=four,regions=regions,commonFaces=len(wi),commonAreaM2=float(area[wi].sum()),commonDeltaW=common_delta,
                jointMinusCapW=replace_delta,baseGeometryDeltaW=original_geom_delta,fineGeometryDeltaW=geom_delta,
                geometryDifferenceChangeW=geom_delta-original_geom_delta,beadGridDeltaW=bead_delta,
                wholeGridDeltaW=summary["totalOutwardW"]-prior["totalOutwardW"],
                fineWholeGeometryDeltaW=summary["totalOutwardW"]-smooth["totalOutwardW"],
                commonReverseChangedFine=int(np.count_nonzero(r[wi]!=sf["reverse_both"][si])),
                reverseChangedAreaGridM2=float(area[r!=old_r].sum()),jointVerifiedAreaM2=0)
    out.mkdir(parents=True,exist_ok=False)
    save(out/"comparison.json",result)
    archive=out/"leadsMidbody3.zip"
    with zipfile.ZipFile(archive,"x") as z:
        for file in (NEW/"body.json",NEW/"case.json",*(NEW/"surface"/n for n in ("summary.json","wall.csv","flow.csv","flow.json"))):
            z.writestr(zipfile.ZipInfo(file.name,date_time=(2026,1,1,0,0,0)),file.read_bytes(),compress_type=zipfile.ZIP_DEFLATED)
    contract={k:inputs["geometry"][k] for k in ("nodes","cells","meshSHA256","bodySHA256","nodeMapSHA256")}
    contract.update(initialSHA256=inputs["inputs"]["initial_flow.csv"],configSHA256=inputs["inputs"]["case.cfg"],
                    caseSHA256=inputs["inputs"]["case.json"],archiveSHA256=sha(archive))
    save(out/"leadsMidbody3.json",contract)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("output",type=Path)
    output=parser.parse_args().output.resolve();output.relative_to(ROOT/"outputs")
    run(output)
