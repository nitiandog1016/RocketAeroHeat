"""读回固定工况的稳态场，核对通量，将子午带研究压力按面积投影至原面。"""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from run_outer import GRIDS, MODEL, ROOT, flight, read_meridian, sha
from read_sphere import face_geometry, listed
from openFoamConservation import array_entry
from openFoamFields import internal_field


def read_run(grid, job, recovery):
    """用原生面序独立复算各单元收支、壁速与压力；逐项读取监督终态。"""
    case = job/"case"
    mesh, fields = case/"constant/polyMesh", case/"0"
    assert not (fields/"phi").samefile(fields/"Phi")
    stages = {}
    for name in ("checkMesh", "potentialFoam"):
        receipt = json.loads((job/name/"report.json").read_bytes())
        assert receipt["passed"] and not receipt["liveGroupAfterCleanup"] and all(receipt["inputStable"].values())
        for path, value in receipt["inputSHA256"].items():
            assert sha(path) == value, path
        stages[name] = {k: receipt[k] for k in ("elapsedSeconds", "nativeExitCode", "peakObservedGroupRSSBytes", "caseBytes")}
    log = (job/"potentialFoam/stdout.txt").read_text()
    mesh_log = (job/"checkMesh/stdout.txt").read_text()
    assert "Mesh OK." in mesh_log and log.rstrip().endswith("End")
    solves = re.findall(r"Solving for Phi, Initial residual = (\S+), Final residual = (\S+), No Iterations (\d+)", log)
    assert len(solves) == 5
    owner, nei = [np.array(listed(mesh/key), int) for key in ("owner", "neighbour")]
    n = int(owner.max()+1)
    origin = ROOT / f"outputs/v350_validation/{grid}01/mesh.su2"
    assert sha(origin) == GRIDS[grid]
    if recovery:
        from run_recovery import recovery_grid
        origin, xy, cells, edges = recovery_grid(grid)
    else:
        xy, cells, edges = read_meridian(origin)
    condition, model = flight()
    xy *= model["maxDiameterM"]
    assert n == len(cells)
    patches = {}
    for name, block in re.findall(r"(\w+)\s*\{([^{}]*)\}", (mesh/"boundary").read_text()):
        if "nFaces" in block:
            start = int(re.search(r"startFace\s+(\d+);", block)[1])
            size = int(re.search(r"nFaces\s+(\d+);", block)[1])
            patches[name] = slice(start, start+size)
    flux_text = (fields/"phi").read_text()
    assert re.search(r"class\s+surfaceScalarField;", flux_text) and re.search(r"object\s+phi;", flux_text)
    assert re.search(r"dimensions\s+\[0 3 -1 0 0 0 0\];", flux_text)
    flux = np.full(len(owner), np.nan)
    flux[:len(nei)] = internal_field(fields/"phi", len(nei))
    for key, selection in patches.items():
        block = re.search(r"\b"+key+r"\s*\{([^}]+)\}", flux_text.split("boundaryField",1)[1])[1]
        vals = array_entry(block, "value", selection.stop-selection.start)
        assert vals is not None and np.isnan(flux[selection]).all()
        flux[selection] = vals
    assert np.isfinite(flux).all()
    net = np.bincount(owner, flux, n) - np.bincount(nei, flux[:len(nei)], n)
    incident = np.bincount(owner, abs(flux), n) + np.bincount(nei, abs(flux[:len(nei)]), n)
    far = flux[patches["far"]]
    incoming = -np.minimum(far, 0).sum()
    assert incoming > 0 and min(incident) > 0
    fluxes = dict(maxCellRelative=float(max(abs(net)/incident)), maxCellM3s=float(max(abs(net))),
                  farNetOverIn=float(far.sum()/incoming), wallAbsOverIn=float(abs(flux[patches["body"]]).sum()/incoming),
                  farInM3s=float(incoming))

    # 从原生网格的壁面顶点核对每一条源壁段，按滑移条件取速度的切向分量。
    pts = np.array([r.strip("()").split() for r in listed(mesh/"points")], float)
    normals = []
    wall = xy[edges["wall"]]
    for row, source in zip(listed(mesh/"faces")[patches["body"]], wall, strict=True):
        size, ids = row.rstrip(")").split("(")
        face = np.array(ids.split(), int)
        assert len(face) == int(size)
        verts = pts[face]
        xr = np.column_stack((verts[:,0], np.hypot(verts[:,1], verts[:,2])))
        delta = np.linalg.norm(xr[:,None]-source, axis=2)
        assert max(delta.min(axis=0)) < 1e-12 and max(delta.min(axis=1)) < 1e-12
        normals.append(face_geometry(pts, face)[1])
    normals = np.array(normals)
    u = internal_field(fields/"U", n, True)[owner[patches["body"]]]
    assert re.search(r"\bbody\s*\{\s*type\s+slip;", (fields/"U").read_text())
    uw = u - normals*np.einsum("ij,ij->i", u, normals)[:,None]
    cp = 1 - np.sum(uw**2, axis=1)/condition["speedMps"]**2
    linear = [dict(initial=float(i), final=float(f), iterations=int(k)) for i,f,k in solves]
    passed = linear[-1]["final"] < 1e-8 and fluxes["maxCellRelative"] < 1e-5 and abs(fluxes["farNetOverIn"]) < 1e-8 and fluxes["wallAbsOverIn"] < 1e-8
    angle = float(re.search(r"non-orthogonality Max: (\S+)", mesh_log)[1])
    severe = re.search(r"severely non-orthogonal .*?faces: (\d+)\.", mesh_log)
    assert angle <= 70 or severe
    receipt = dict(stages=stages, sourceMeshSHA256=sha(origin), cells=n, phiResiduals=linear, flux=fluxes,
                   nonOrthogonalMaxDeg=angle, nonOrthogonalAbove70=int(severe[1]) if severe else 0,
                   passesRegisteredGates=bool(passed), nativeSHA256={k:sha(fields/k) for k in ("U","phi","Phi")})
    return wall, cp, receipt


def overlap(model, wall):
    """在 STL 的平面三角形上积分分段常数。底面同心相似多边形的面积随半径平方变化。"""
    raw = np.dtype([("normal", "<f4", 3), ("v", "<f4", (3,3)), ("att", "<u2")])
    tri = np.fromfile(MODEL/model["file"], dtype=raw, offset=84)["v"].astype(float)
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]), axis=1)/2
    assert len(area) == model["faces"] and abs(area.sum()-model["areaM2"]) < 1e-12
    base = wall[:,1,0] == 0
    bounds = wall[~base,:,0]
    weights = np.zeros((len(area),len(wall)))
    for i, vertices in enumerate(tri):
        x = vertices[:,0]
        if max(x) == 0:
            radii = np.linalg.norm(vertices[:,1:], axis=1)
            assert np.count_nonzero(radii) == 2 and abs(radii.max()-wall[base,-1,1].max()) < 1e-8
            weights[i,base] = np.diff((wall[base,:,1]/wall[base,-1,1].max())**2, axis=1).ravel()
        else:
            lo, hi = min(x), max(x)
            assert len(np.unique(x)) == 2
            t = np.clip((bounds-lo)/(hi-lo), 0, 1)
            cdf = t**2 if np.count_nonzero(x == lo) == 1 else 1-(1-t)**2
            weights[i,~base] = np.diff(cdf,axis=1).ravel()
    np.testing.assert_allclose(weights.sum(axis=1), 1, atol=1e-12, rtol=0)
    assert weights.min() >= 0 and (weights.sum(axis=0) > 0).all()
    return area, weights


def summarize(folders, recovery):
    """比较同一批原面的两域差额，输出可还原原面编号的研究压力表。"""
    condition, model = flight()
    rho = condition["freePressurePa"] / (287.05287*condition["freeTemperatureK"])
    dynamic = 0.5*rho*condition["speedMps"]**2
    runs = [read_run(key, path, recovery) for key,path in zip(GRIDS, folders, strict=True)]
    np.testing.assert_array_equal(runs[0][0], runs[1][0])
    area, weights = overlap(model, runs[0][0])
    # 同一轴向三角形形状的周向面共享面积权重；反向索引显式保留 STL 原面序。
    groups, ids = np.unique(weights,axis=0,return_inverse=True)
    by_group = np.stack([groups @ cp for _,cp,_ in runs], axis=1)
    face_cp = by_group[ids]
    assert np.all(np.isfinite(face_cp))
    delta = face_cp[:,1] - face_cp[:,0]
    metrics = {}
    for r in [*model["regions"], dict(name="all", firstFace=1, lastFace=len(area))]:
        s = slice(r["firstFace"]-1, r["lastFace"])
        weights_area = area[s]
        cp = face_cp[s,1]
        metrics[r["name"]] = dict(faces=len(cp), areaM2=float(weights_area.sum()), cpMin=float(min(cp)), cpMax=float(max(cp)),
            cpMean=float(np.average(cp,weights=weights_area)), rmsDeltaCp=float(np.sqrt(np.average(delta[s]**2,weights=weights_area))),
            maxAbsDeltaCp=float(max(abs(delta[s]))))
    audit = {k:r[2] for k,r in zip(GRIDS,runs,strict=True)}
    for i,(key,(_,cp,receipt)) in enumerate(zip(GRIDS,runs,strict=True)):
        pressure = condition["freePressurePa"] + dynamic*face_cp[:,i]
        report = dict(minBandAbsolutePa=float(condition["freePressurePa"]+dynamic*min(cp)),
                      minFaceAbsolutePa=float(min(pressure)), nonPositiveFaces=int(np.sum(pressure<=0)),
                      maxBandSpeedMps=float(condition["speedMps"]*np.sqrt(1-min(cp))),
                      maxBandSpeedOverFreeSound=float(0.3*np.sqrt(1-min(cp))))
        receipt["pressureCheck"] = report
        if recovery and (report["minBandAbsolutePa"]<=0 or report["nonPositiveFaces"]):
            raise ValueError(f"{key} 原始压力非正 {report}")
    if not all(r["passesRegisteredGates"] for r in audit.values()):
        raise ValueError("原生场未通过事前数值门槛：" + json.dumps(audit))
    band_table = np.column_stack((runs[0][0].reshape(-1,4), runs[0][1], runs[1][1]))
    data = dict(surfaceSHA256=model["fileSHA256"],mach=0.3,altitudeM=5000,alphaDeg=0,betaDeg=0,wallGivenK=300,
                freePressurePa=condition["freePressurePa"],dynamicPressurePa=dynamic,regions=metrics,
                faceGroup=(ids+1).tolist(),groupCp14D=by_group[:,0].tolist(),groupCp22D=by_group[:,1].tolist(),
                qualifiedPressureFaces=0,heatFaces=0,wholeBodyPowerW=None)
    return data, audit, band_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fine", type=Path)
    parser.add_argument("outer", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--recovery", action="store_true",help="读取保留原壁和边界的三角网对照")
    args = parser.parse_args()
    result, audit, table = summarize([args.fine, args.outer],args.recovery)
    args.output.mkdir(parents=True,exist_ok=False)
    for name, data in (("surface", result), ("numerics", audit)):
        (args.output/f"{name}.json").write_text(json.dumps(data,allow_nan=False,separators=(",",":"))+"\n",encoding="ascii")
    np.savetxt(args.output/"bands.csv",table,delimiter=",",header="x0M,r0M,x1M,r1M,cp14D,cp22D",comments="")
    print(json.dumps(result["regions"],indent=2))
