"""从势流算例的四个原生场和边界导出分别核对两类压力、壁速和面通量。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]/"scripts"))
from openFoamConservation import array_entry
from openFoamFields import internal_field


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def listed(path):
    match = re.search(r"^\s*(\d+)\s*\n\(\s*\n(.*?)^\)",path.read_text(),re.M|re.S)
    if not match:
        raise ValueError(f"{path.name} 缺少完整列表")
    rows = match[2].strip().splitlines()
    assert len(rows) == int(match[1])
    return rows


def vtk_values(path):
    """仅读取原生转换工具输出的面数据；拒绝转为点数据的插值。"""
    text = path.read_text()
    assert "DATASET POLYDATA" in text and "POINT_DATA" not in text
    tokens = iter(text.split("CELL_DATA",1)[1].split())
    n = int(next(tokens))
    assert n == 128 and next(tokens) == "FIELD" and next(tokens) == "attributes"
    fields = {}
    for _ in range(int(next(tokens))):
        name, width = next(tokens), int(next(tokens))
        assert int(next(tokens)) == n and next(tokens) in ("float","double")
        fields[name] = np.array([float(next(tokens)) for _ in range(n*width)]).reshape(n,width)
    return fields


def face_geometry(points, face):
    vertices = points[face]
    origin = vertices.mean(axis=0)
    following = np.roll(vertices,-1,axis=0)
    areas = np.cross(vertices-origin,following-origin)/2
    weights = np.linalg.norm(areas,axis=1)
    centroid = np.average((origin+vertices+following)/3, axis=0, weights=weights)
    normal = areas.sum(axis=0)
    return centroid, normal/np.linalg.norm(normal)


def metrics(angle, weight, velocity, bernoulli, poisson):
    analytical_speed = 1.5*np.sin(angle)
    analytical_cp = 1-analytical_speed**2
    def rms(delta):
        return float(np.sqrt(np.average(delta**2,weights=weight)))
    result = {"tangentialRMS":rms(velocity-analytical_speed)}
    for name, cp in (("bernoulli",bernoulli),("nativePoisson",poisson)):
        delta = cp-analytical_cp
        worst = int(np.argmax(abs(delta)))
        result[name] = {"rmsCp":rms(delta),"maxAbsCp":float(abs(delta[worst])),
                        "worstAngleDeg":float(np.degrees(angle[worst]))}
        result[name]["passesWallGates"] = bool(rms(delta)<0.04 and abs(delta[worst])<0.10
                                              and result["tangentialRMS"]<0.03)
    return result


def common_nodes(angles, samples, power):
    """回到 v361 的共同极角；极点外推仅用数值面带，保留外推身份。"""
    nodes = np.linspace(0,np.pi,129)
    projected = np.interp(nodes, angles, samples)
    for index, ids in ((0,[0,1]),(-1,[-1,-2])):
        distances = abs(angles[ids]-nodes[index])**power
        v = samples[ids]
        projected[index] = (distances[1]*v[0]-distances[0]*v[1])/(distances[1]-distances[0])
    return projected


def read_case(job, csv_file):
    case = job/"case"
    fields = case/"0"
    assert not (fields/"phi").samefile(fields/"Phi")
    identities = {"U":("volVectorField","0 1 -1"),"p":("volScalarField","0 2 -2"),
                  "phi":("surfaceScalarField","0 3 -1"),"Phi":("volScalarField","0 2 -1")}
    for name,(kind,dim) in identities.items():
        data = (fields/name).read_text()
        assert re.search(r"\bobject\s+"+name+r";",data)
        assert re.search(r"\bclass\s+"+kind+r";",data)
        assert re.search(r"dimensions\s+\["+dim+r" 0 0 0 0\];",data)
    mesh = case/"constant/polyMesh"
    owner, neighbour = (np.array(listed(mesh/key),dtype=int) for key in ("owner","neighbour"))
    cells = int(owner.max()+1)
    assert cells == 128*72 and neighbour.max()<cells
    patches = {}
    for name,body in re.findall(r"(\w+)\s*\{([^{}]*)\}",(mesh/"boundary").read_text()):
        if "nFaces" in body:
            n = int(re.search(r"nFaces\s+(\d+);",body)[1])
            start = int(re.search(r"startFace\s+(\d+);",body)[1])
            patches[name] = slice(start,start+n)
    flux = np.full(len(owner),np.nan)
    flux[:len(neighbour)] = internal_field(fields/"phi",len(neighbour))
    boundary = (fields/"phi").read_text().split("boundaryField",1)[1]
    flux_totals = {}
    for name, ids in patches.items():
        n = ids.stop-ids.start
        if n:
            patch = re.search(r"\b"+name+r"\s*\{([^{}]*)\}",boundary,re.S)[1]
            values = array_entry(patch,"value",n)
            assert values is not None and len(values)==n and np.isnan(flux[ids]).all()
            flux[ids] = values
            flux_totals[name] = float(values.sum())
    assert np.isfinite(flux).all()
    imbalance = np.bincount(owner,weights=flux,minlength=cells)-np.bincount(neighbour,weights=flux[:len(neighbour)],minlength=cells)
    incoming = float(-np.minimum(flux[patches["far"]],0).sum())
    incident = np.bincount(owner,weights=abs(flux),minlength=cells)+np.bincount(neighbour,weights=abs(flux[:len(neighbour)]),minlength=cells)
    assert incoming>0 and (incident>0).all()

    points = np.array([row.strip("()").split() for row in listed(mesh/"points")],dtype=float)
    rows = listed(mesh/"faces")[patches["body"]]
    geom = []
    for row in rows:
        size,ids = row.rstrip(")").split("(")
        face = np.fromstring(ids,sep=" ",dtype=int)
        assert len(face)==int(size) and np.allclose(np.linalg.norm(points[face],axis=1),.5)
        geom.append(face_geometry(points,face))
    center, normal = map(np.array, zip(*geom))
    ids = owner[patches["body"]]
    native_u = internal_field(fields/"U",cells,True)[ids]
    wall_u = native_u - normal*np.sum(native_u*normal,axis=1)[:,None]
    native_p = internal_field(fields/"p",cells)[ids]
    export = vtk_values(case/"VTK/body/body_0.vtk")
    np.testing.assert_allclose(wall_u,export["U"],atol=6e-6,rtol=0)
    np.testing.assert_allclose(native_p,export["p"][:,0],atol=6e-6,rtol=0)
    theta = np.arctan2(center[:,1],center[:,0])
    order = theta.argsort()
    theta,wall_u,native_p = theta[order],wall_u[order],native_p[order]
    edges = np.linspace(0,np.pi,129)
    assert np.all((theta>edges[:-1]) & (theta<edges[1:]))
    tangent = -wall_u[:,0]*np.sin(theta) + wall_u[:,1]*np.cos(theta)
    cp_b = 1-np.sum(wall_u**2,axis=1)
    cp_p = 2*native_p
    bands = -np.diff(np.cos(edges))
    node_weight = (np.r_[bands,0]+np.r_[0,bands])/2
    node_fields = [common_nodes(theta,values,power) for values,power in ((tangent,1),(cp_b,2),(cp_p,2))]
    if csv_file:
        with csv_file.open("x",encoding="ascii",newline="") as out:
            writer = csv.writer(out,lineterminator="\n")
            writer.writerow(["angleDeg","ringWeight","utOverU","cpBernoulli","cpNativePoisson"])
            writer.writerows(zip(np.degrees(theta),bands,tangent,cp_b,cp_p))
    stages = {}
    for stage in ("blockMesh","checkMesh","potentialFoam","foamToVTK"):
        r = json.loads((job/stage/"report.json").read_bytes())
        assert r["passed"] and not r["liveGroupAfterCleanup"] and all(r["inputStable"].values())
        stages[stage] = {k:r[k] for k in ("elapsedSeconds","nativeExitCode","stopReason","peakObservedGroupRSSBytes","caseBytes")}
    log = (job/"potentialFoam/stdout.txt").read_text()
    residuals = {}
    for field in ("p","Phi"):
        matches = re.findall(r"Solving for "+field+r", Initial residual = ([\deE.+-]+), Final residual = ([\deE.+-]+), No Iterations (\d+)",log)
        assert len(matches)==5
        residuals[field] = [{"initial":float(a),"final":float(b),"iterations":int(n)} for a,b,n in matches]
    return dict(cells=cells,wallBands=len(theta),stages=stages,faceCenters=metrics(theta,bands,tangent,cp_b,cp_p),
                commonNodes=metrics(edges,node_weight,*node_fields),
                poleCp={k:[float(a[0]),float(a[-1])] for k,a in zip(("bernoulli","nativePoisson"),node_fields[1:])},
                fluxM3s=flux_totals,incomingFarM3s=incoming,maxCellImbalanceM3s=float(max(abs(imbalance))),
                maxCellImbalanceOverIncident=float(max(abs(imbalance)/incident)),farNetOverInflow=flux_totals["far"]/incoming,
                nativeContinuity=float(re.search(r"Continuity error = (\S+)",log)[1]),
                normalWallSpeedMax=float(np.max(np.abs(np.sum(wall_u*normal[order],axis=1)))),
                linearSolves=residuals,linearGatePassed=all(v[-1]["final"]<1e-8 for v in residuals.values()),
                outputSHA256={name:digest(fields/name) for name in identities},
                meshSHA256={p.name:digest(p) for p in mesh.iterdir() if p.is_file()})


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job",type=Path)
    parser.add_argument("--wall-csv",type=Path)
    args=parser.parse_args()
    print(json.dumps(read_case(args.job,args.wall_csv),indent=2,allow_nan=False))
