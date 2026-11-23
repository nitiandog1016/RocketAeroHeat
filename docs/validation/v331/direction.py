"""Short Test 65 双支管银球：把两个已保存端点的近壁向量配对到相同原面。"""
import argparse
import io
import json
import os
from pathlib import Path
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

for env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[env] = "1"
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = list(map(str, (HERE.parent/"v329", HERE.parent/"v330", ROOT/"scripts")))
from map_surface import corner_weights
from compare import load_surface
from check_grid import COARSE, COARSE_KEY, FINE_KEY, sha
from auditSU2SSTState import vtk_tokens

DATA = ROOT/"cases/short_sphere"
PAIR = DATA/"supportedNearWall65.json"
FLAGS = ["axis", "noDirection", "targetEdge", "sourceEdge", "smallShear",
         "smallSpeed", "disagreement", "mixedShear", "mixedSpeed"]
EPS = 64*np.finfo(float).eps


def read(path):
    return json.loads(path.read_bytes())


def norm(v):
    return np.linalg.norm(v, axis=-1)


def tangent(v, normal):
    return v - np.sum(v*normal, axis=-1, keepdims=True)*normal


def basis(normal):
    direction = tangent(np.array([-1., 0., 0.]), normal)
    length = norm(direction)
    good = length > 1e-10
    return direction/np.where(good, length, 1.)[..., None], good


def edge_contact(corners, regions):
    # 跨区共享点影响该点相接的所有面；重复的三角片顶点也得到同一标签。
    ids = corners.ravel()
    groups = np.repeat(regions, corners.shape[1])
    lo = np.full(ids.max()+1, 6); hi = np.zeros_like(lo)
    np.minimum.at(lo, ids, groups); np.maximum.at(hi, ids, groups)
    return np.any(lo[corners] != hi[corners], axis=1)


def surface_graph(f):
    edges = np.sort(f[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    np.testing.assert_array_equal(edges[order[::2]], edges[order[1::2]])
    a, b = order[::2]//3, order[1::2]//3
    graph = coo_matrix((np.ones(len(a)*2), (np.r_[a,b], np.r_[b,a])), shape=(len(f),len(f))).tocsr()
    assert np.all(np.diff(graph.indptr)==3) and (graph.data==1).all()
    return graph


def geometry(pair,transfer_file=DATA/"supportedTransfer65.npz"):
    assert sha(DATA/"supported65.zip")==pair["bodySHA256"]
    with ZipFile(DATA/"supported65.zip") as archive:
        body = json.loads(archive.read("geometry.json"))
    assert sha(transfer_file)==pair["transferSHA256"]
    with np.load(transfer_file) as handle:
        m = dict(handle)
    f = np.array(body["F"])-1; v = np.array(body["V"]); p = v[f]
    normal = np.cross(p[:,1]-p[:,0], p[:,2]-p[:,0])
    np.testing.assert_allclose(norm(normal)/2, m["targetArea"], rtol=1e-11)
    normal /= norm(normal)[:,None]
    c = m["nativeCorners"]; xyz = m["nativeWallXYZD"]; quad = xyz[c]
    sn = np.cross(quad[:,2]-quad[:,0], quad[:,3]-quad[:,1])
    sn /= norm(sn)[:,None]
    source, target, weight = (m[k] for k in ("source","target","targetFraction"))
    np.testing.assert_array_equal(m["region"][target], m["nativeRegion"][source])
    assert len(np.unique(np.column_stack((target,source)), axis=0))==len(source) and min(weight)>0
    np.testing.assert_allclose(np.bincount(target, weight), 1., rtol=0, atol=2e-12)
    axis = np.any(norm(p[:,:,1:])<1e-12, axis=1) & (m["region"]==1)
    source_edge = edge_contact(c,m["nativeRegion"])[source]
    return dict(m=m, normal=normal, sn=sn, graph=surface_graph(f), weights=corner_weights(c,xyz),
                axis=axis, targetEdge=edge_contact(f,m["region"]),
                sourceEdge=np.bincount(target,source_edge,minlength=len(f))>0)


def native(name, pair, grid):
    coarse = name=="coarse"
    prep = COARSE if coarse else ROOT/"outputs/v330_validation/prepared01"
    stage = ROOT/("outputs/v329_validation" if coarse else "outputs/v330_validation")/"run01/steady03"
    wall_key = read(HERE.parent/"v330/evidence/wall_units.json")[name]
    assert sha(prep/"prepared.json")==(COARSE_KEY if coarse else FINE_KEY)
    assert wall_key["frame"]==pair[name+"Frame"]
    return native_vectors(prep,stage,grid,pair["condition"],wall_key)


def native_vectors(prep,stage,grid,condition,wall_key):
    spec = read(prep/"prepared.json")
    assert sha(prep/"mapping.npz")==spec["files"]["mapping.npz"] and spec["condition"]==condition
    record = read(stage/"record.json"); proc = read(stage/"debugger_process.json")
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not (proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
    assert "Exit Success (SU2_CFD)" in (stage/"debugger.log").read_text(errors="replace")
    assert max(record["window"])<-8 and sha(stage/"history.csv")==record["history"]
    frame = stage/"restart_000199.csv"; vtk = stage/"surface_000199.vtk"
    assert record["end"]==sha(frame)==wall_key["frame"] and sha(vtk)==wall_key["vtkSHA256"]
    with np.load(prep/"mapping.npz") as m:
        wet, first, xyz = m["wet"], m["first"], m["xyz"]
        np.testing.assert_array_equal(wet[grid["nativeCorners"]], m["corners"])
        np.testing.assert_array_equal(xyz[wet], grid["nativeWallXYZD"])
    points, fields = vtk_tokens(vtk)
    dist, idx = cKDTree(points).query(xyz[wet])
    assert max(dist)<1e-6 and len(np.unique(idx))==len(wet)==len(points)
    state = np.loadtxt(frame, delimiter=",", skiprows=1, usecols=range(8))
    np.testing.assert_array_equal(state[:,0], np.arange(len(xyz)))
    np.testing.assert_allclose(state[:,1:4],xyz,atol=2e-13,rtol=2e-14)
    np.testing.assert_allclose(fields["Pressure"][idx,0],state[wet,4],atol=1e-7,rtol=7e-7)
    cfg = spec["condition"]; speed = cfg["speedMps"]; dyn = .5*cfg["rho"]*speed**2
    vectors = [fields["Skin_Friction_Coefficient"][idx]*dyn, state[first,5:8]*speed]
    assert all(np.isfinite(x).all() for x in vectors)
    return vectors, wet, [dyn,speed], wall_key


def classify(components, tolerance):
    s = np.where(np.abs(components)>tolerance,np.sign(components),0).astype(int)
    return s, np.where(s[0]==s[1],s[0],0)


def mapped_fields(geom, vectors, scales):
    m = geom["m"]; t, s, w = (m[k] for k in ("target","source","targetFraction"))
    n = geom["normal"]; forward, good = basis(n)
    flags = dict(axis=geom["axis"], noDirection=~good, targetEdge=geom["targetEdge"], sourceEdge=geom["sourceEdge"])
    components, magnitudes, signs = [], [], []
    for k, values in enumerate(vectors):
        corners = values[m["nativeCorners"]]
        source = np.sum(geom["weights"][:,:,None]*tangent(corners,geom["sn"][:,None,:]), axis=1)
        # 对容差使用原输入模的平均值，抵消项不能把舍入带人为缩小。
        amplitude = np.sum(geom["weights"]*norm(corners),axis=1)
        precision = 1e-6 if k==0 else EPS
        bound = precision*amplitude + EPS*scales[k]
        projected = tangent(np.column_stack([np.bincount(t,w*source[s,j]) for j in range(3)]),n)
        along = np.sum(projected*forward,axis=1)
        sign = np.where(abs(along)>np.bincount(t,w*bound[s]), np.sign(along), 0).astype(int)
        part = np.sum(source[s]*forward[t],axis=1)
        both = (np.bincount(t,part>bound[s])>0) & (np.bincount(t,part< -bound[s])>0)
        flags[("mixedShear","mixedSpeed")[k]] = both & good
        flags[("smallShear","smallSpeed")[k]] = (sign==0) & good
        along[~good] = np.nan
        signs.append(sign); components.append(along); magnitudes.append(norm(projected))
    flags["disagreement"] = (signs[0]*signs[1]<0) & good
    bits = sum(flags[key].astype(int)<<j for j,key in enumerate(FLAGS))
    status = np.where(bits==0,signs[0],0)
    assert np.all(signs[0][bits==0]==signs[1][bits==0]) and np.all(status[bits==0]!=0)
    return np.column_stack((components[0],magnitudes[0],components[1],magnitudes[1],status,bits))


def station_point(ids, weights, degrees, vectors, scales):
    angle = np.deg2rad(degrees); n = np.array([np.cos(angle),0.,np.sin(angle)])
    forward, good = basis(n)
    signed, magnitude, bound = [], [], []
    for k, field in enumerate(vectors):
        v = tangent(weights@field[ids],n)
        signed.append(v@forward); magnitude.append(norm(v))
        bound.append((1e-6 if k==0 else EPS)*(weights@norm(field[ids]))+EPS*scales[k])
    sign, status = classify(signed,bound)
    reason = "" if status else ("disagreement" if sign.prod()<0 else "nearZero")
    return dict(angleDeg=float(degrees), status=int(status) if good else 0, reason=reason if good else "noDirection",
                shearPa=signed[0] if good else None, speedMps=signed[1] if good else None,
                shearMagnitudePa=magnitude[0], speedMagnitudeMps=magnitude[1])


def subtotal(mask, area, flux):
    power = -flux[mask]*area[mask]
    return dict(faces=int(mask.sum()), areaM2=area[mask].sum(), outW=power.sum(),
                releaseW=power[power>0].sum(), absorbedW=-power[power<0].sum())


def by_status(mask,status,area,flux):
    return {label:subtotal(mask & (status==k),area,flux) for label,k in (("forward",1),("reverse",-1),("unknown",0))}


def report_grid(name, geom, values, wet, vectors, scales, wall, stations, pair,sectors=64):
    m = geom["m"]; area, region = m["targetArea"],m["region"]; flux = wall["fluxInWm2"]
    status = values[:,4]; bits = values[:,5].astype(int); graph = geom["graph"]
    info = dict(regions={}, flags={}, stations=[], transitions=[], components=[])
    all_faces = np.ones(len(area),bool)
    info["total"] = by_status(all_faces,status,area,flux)
    assert abs(sum(x["outW"] for x in info["total"].values())-pair[name+"W"])<1e-10
    for i in range(1,6):
        info["regions"][str(i)] = by_status(region==i,status,area,flux)
    for j,key in enumerate(FLAGS):
        info["flags"][key] = subtotal((bits & (1<<j))!=0,area,flux)
    rev = np.flatnonzero(status==-1)
    count, labels = connected_components(graph[rev][:,rev],directed=False)
    for label in range(count):
        ids = rev[labels==label]; mask = np.zeros(len(area),bool); mask[ids]=True
        info["components"].append(dict(firstFaceID=int(ids.min()+1),regions=np.unique(region[ids]).tolist(),**subtotal(mask,area,flux)))
    info["components"].sort(key=lambda x:-x["areaM2"])
    def ring_stats(ids):
        adjacent = np.setdiff1d(np.unique(graph[ids].indices),ids)
        return dict(core=by_status(np.isin(np.arange(len(area)),ids),status,area,flux),
                    neighbours=by_status(np.isin(np.arange(len(area)),adjacent),status,area,flux))
    peak = pair[name+"Peak"]; ids = np.array(peak["nearFaces"])-1
    np.testing.assert_array_equal(ids,np.flatnonzero(-flux>=.99*max(-flux)))
    info["nearPeak"] = ring_stats(ids)
    j = peak["faceID"]-1
    info["peak"] = dict(faceID=j+1,status=int(status[j]),flags=int(bits[j]),**ring_stats([j]))
    xyz = m["nativeWallXYZD"]
    for row in stations:
        node = np.array([row["leftNode"],row["rightNode"]],int); w = row["rightWeight"]; mix = np.array([1-w,w])
        idx = np.searchsorted(wet,node); np.testing.assert_array_equal(wet[idx],node)
        p = mix@xyz[idx]
        assert abs(np.rad2deg(np.arctan2(p[2],p[0]))-row["angleDeg"])<1e-10 and abs(p[1])<1e-12
        point = station_point(idx,mix,row["angleDeg"],vectors,scales)
        face = int(row["originalFaceID"])-1
        point.update(faceID=face+1,faceStatus=int(status[face]),faceFlags=int(bits[face]),leftNode=int(node[0]),
                     rightNode=int(node[1]),rightWeight=w,computedNu=row["computedNu"],measuredNu=row["reportedNu"],
                     **ring_stats([face]))
        info["stations"].append(point)
    sphere = np.unique(m["nativeCorners"][m["nativeRegion"]==1])
    ring = sphere[abs(xyz[sphere,1])<1e-10]
    # 带有迎、背流轴点的测线仅在方向有定义的开半环上寻找符号变化。
    axial = abs(xyz[ring,2])<=1e-12
    line = ring[xyz[ring,2]>1e-12]
    phi = np.rad2deg(np.arctan2(xyz[line,2],xyz[line,0])); order=np.argsort(phi)
    assert len(ring)==sectors and 2*len(line)+axial.sum()==len(ring)
    trace = [station_point([j],np.ones(1),p,vectors,scales) for j,p in zip(line[order],phi[order])]
    for left,right in zip(trace,trace[1:]):
        if left["status"]!=right["status"]:
            info["transitions"].append(dict(left=left,right=right))
    return info


def main(out):
    pair = read(PAIR)
    assert sha(PAIR)=="d5c97ea022d968533a7e0b98c0b7cd599b14d6aa106860c575fedd8952df82d3"
    geom = geometry(pair); m = geom["m"]; area = m["targetArea"]
    result = dict(caseSHA256=pair["caseSHA256"],bodySHA256=pair["bodySHA256"], transferSHA256=pair["transferSHA256"],
                  areaM2=area.sum(),faces=len(area),flagNames=FLAGS,positiveHeat="intoWall")
    arrays = []; reverse = []
    for name,archive in (("coarse","supportedSurface65.zip"),("fine","supportedFine65.zip")):
        assert sha(DATA/archive)==pair[name+"ArchiveSHA256"]
        summary,wall,stations = load_surface(DATA/archive)
        np.testing.assert_array_equal(wall["faceID"],np.arange(1,len(area)+1))
        np.testing.assert_array_equal(wall["regionID"],m["region"]); np.testing.assert_array_equal(wall["areaM2"],area)
        vectors,wet,scales,keys = native(name,pair,m)
        assert keys["frame"]==summary["frameSHA256"]
        values = mapped_fields(geom,vectors,scales)
        result[name] = report_grid(name,geom,values,wet,vectors,scales,wall,stations,pair)
        result[name].update(frame=keys["frame"],vtk=keys["vtkSHA256"],dynamicPressurePa=scales[0],speedReferenceMps=scales[1])
        arrays.append(values); reverse.append(values[:,4]==-1)
    both = reverse[0]&reverse[1]; either = reverse[0]|reverse[1]
    result["overlap"] = dict(bothFaces=int(both.sum()),bothM2=area[both].sum(),eitherM2=area[either].sum(),
                             coarseOnlyM2=area[reverse[0]&~reverse[1]].sum(),fineOnlyM2=area[reverse[1]&~reverse[0]].sum())
    out.mkdir(parents=True,exist_ok=False)
    header = ["faceID","regionID"]+[name+k for name in ("coarse","fine") for k in
                    ("ShearPa","ShearMagnitudePa","SpeedMps","SpeedMagnitudeMps","Status","Flags")]
    text = io.StringIO()
    fmt = ["%d","%d"]+(["%.17g"]*4+["%d","%d"])*2
    np.savetxt(text,np.column_stack((np.arange(1,len(area)+1),m["region"],*arrays)),delimiter=",",fmt=fmt,header=",".join(header),comments="")
    file = out/"supportedDirection65.zip"
    with ZipFile(file,"x") as z:
        z.writestr(ZipInfo("direction.csv"),text.getvalue(),compress_type=ZIP_DEFLATED,compresslevel=9)
    result["archiveSHA256"]=sha(file)
    with (out/"supportedDirection65.json").open("x",encoding="ascii",newline="\n") as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write("\n")
    print(result["overlap"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("output",type=Path)
    main(parser.parse_args().output)
