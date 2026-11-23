"""带线 B 珠 96 档专用：对同一标准面内的 STL 分片积分真实细壁节点场。"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

from thermistorAngularSpacing import overlap_weights
from thermistorBead import ROOT, read, sha
from thermistorWireSurface import residual

CASES = ROOT / "cases/thermistor_b"
PACK = CASES / "wallNodes96"


def extract(stage, frame, reaction, grid, spec):
    """按冻结原生帧和释放壁反力的节点顺序取压、热流及两类流体超温。"""
    field = stage / "restart.csv"
    rhs = stage / "reaction0/restart.csv"
    if sha(field) != frame or sha(rhs) != reaction:
        raise ValueError("原生终帧或壁反力指纹不符")
    state = np.loadtxt(field, delimiter=",", skiprows=1, usecols=(0,4,9))
    np.testing.assert_array_equal(state[:,0], np.arange(len(grid["xyz"])))
    raw = residual(rhs)
    wet = np.unique(grid["wallMap"])
    first = np.unique(grid["firstMap"])
    lump = np.bincount(grid["wallMap"].ravel(), weights=np.repeat(grid["area"]/3,3))
    thermal = spec.get("scalarReferenceDeltaK",spec["wallK"]-spec["airK"])
    scale = spec["rho"]*spec["speedMps"]*spec["cp"]*thermal*spec["diameterM"]**2
    values = dict(gauge=spec["rho"]*spec["speedMps"]**2*state[wet,1],
                  flux=scale*raw[wet]/lump[wet], wall=thermal*state[wet,2], first=thermal*state[first,2])
    assert all(np.all(np.isfinite(v)) for v in values.values())
    return values


def nodes(source, grid, metadata, mesh_hash, map_hash):
    frame=source["frame"]
    if not re.fullmatch("[a-f0-9]{64}",frame):
        raise ValueError("未登记的终帧格式")
    if frame in metadata:
        entry=metadata[frame];path=PACK/(frame+".npz")
        assert entry["mesh"]==mesh_hash and entry["nodeMap"]==map_hash
        assert not source["reaction"] or entry["reaction"]==source["reaction"]
        assert sha(path)==entry["sha256"]
        with np.load(path,allow_pickle=False) as bundle:
            values={k:bundle[k] for k in ("gauge","flux","wall","first")}
        return values,dict(frame=frame,reaction=entry["reaction"],nodeFieldsSHA256=entry["sha256"])
    # 新生成的原生帧直接从同一任务读取；不以最接近的登记算例代替。
    stage=Path(source["directory"]).resolve();stage.relative_to((ROOT/"outputs").resolve())
    prep=stage if (stage/"case.json").is_file() else stage.parent/"prepared"
    summary=stage/"surface/summary.json"
    if not summary.is_file():summary=stage.parent/"surface/summary.json"
    s=read(summary);reaction=s.get("sourceReactionSHA256",s.get("reactionSHA256"))
    assert s["sourceFrameSHA256"]==frame and (not source["reaction"] or reaction==source["reaction"])
    assert s["meshSHA256"]==mesh_hash and s["nodeMapSHA256"]==map_hash
    assert sha(prep/"case.json")==s["caseSHA256"]
    return extract(stage,frame,reaction,grid,read(prep/"case.json")),dict(frame=frame,reaction=reaction)


def project(request):
    """计算子面形函数积分，并按父面的每项已验收权重分别核对。"""
    prep=Path(request["prepared"]);manifest=read(prep/"prepared.json")
    assert sha(prep/"nodes.npz")==manifest["nodeMapSHA256"]==request["nodeMap"]
    assert sha(prep/"body.json")==manifest["bodySHA256"]==request["body"]
    assert manifest["meshSHA256"]==request["mesh"]
    body=read(prep/"body.json");spec=read(prep/"case.json")
    standard=np.asarray(body["V"])[np.asarray(body["F"])-1]
    parent=np.atleast_1d(request["parents"]).astype(int)-1
    area=np.atleast_1d(request["area"])
    bary=np.asarray(request["bary"]).reshape(-1,3,3);target=bary@standard[parent]
    with np.load(prep/"nodes.npz",allow_pickle=False) as grid:
        xyz=grid["xyz"]*spec["diameterM"];wall=grid["wallMap"];first=grid["firstMap"]
        src,near=xyz[wall],xyz[first];fine_area=grid["area"];old=grid["originalArea"]
        pi,ci,original=(grid[k] for k in ("overlapParent","overlapChild","shapeAreaM2"))
        shape=[];child=[];fine=[];worst=0.
        for j in np.unique(parent):
            ids=np.flatnonzero(parent==j)
            for e in np.flatnonzero(pi==j):
                k=ci[e];total=np.zeros(3)
                for i in ids:
                    w=overlap_weights(target[i],src[k],fine_area[k]);total+=w
                    if np.any(w):
                        shape.append(w);child.append(i);fine.append(k)
                worst=max(worst,float(np.max(abs(total-original[e]))/old[j]))
        assert worst<1e-10,("子面未覆盖父面的细壁形函数积分",worst)
        weights=np.asarray(shape);child=np.asarray(child);fine=np.asarray(fine)
        portions=np.bincount(child,weights=weights.sum(axis=1),minlength=len(target))
        expected=np.linalg.det(bary)*old[parent]
        np.testing.assert_allclose(portions,expected,atol=0,rtol=1e-9)
        np.testing.assert_allclose(portions,area,atol=0,rtol=1.0001e-5)
        def integrate(corners):
            return np.bincount(child, np.einsum("ij,ij->i",weights,corners[fine]),minlength=len(target))
        wi=np.searchsorted(np.unique(wall),wall);ni=np.searchsorted(np.unique(first),first)
        centers=np.column_stack([integrate(near[:,:,d])/portions for d in range(3)])
        metadata=read(PACK/"manifest.json")
        fields={};sources={}
        for name,source in request["sources"].items():
            val,sources[name]=nodes(source,grid,metadata,request["mesh"],request["nodeMap"])
            fields[name]=dict(gauge=(integrate(val["gauge"][wi])/portions).tolist(),
                              flux=(integrate(val["flux"][wi])/area).tolist(),
                              wall=(integrate(val["wall"][wi])/portions).tolist(),
                              first=(integrate(val["first"][ni])/portions).tolist())
    return dict(fields=fields,sources=sources,firstXYZ=centers.tolist(),shapeError=worst,
                overlaps=len(fine),childID=request["childID"],stlSHA256=request["stlSHA256"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("request",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();result=project(read(a.request))
    with a.output.open("x",encoding="utf8") as f:json.dump(result,f,allow_nan=False)
