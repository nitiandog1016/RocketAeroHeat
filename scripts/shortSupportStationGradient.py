"""Short Test 65 银球测线：从实际体网的连通温度层独立重建流体侧法向热流。"""
from collections import defaultdict
import hashlib
import io
import json
from zipfile import ZipFile

import numpy as np

from shortSphere65 import CASE, conductivity, read, sha
from thermistorLeadMesh import SIDES
from thermistorOuterDomain import load_mesh


class Unavailable(ValueError):
    """该测站的实际网格不能唯一支持所请求的法向模板。"""


def measurement_grid(prepared):
    """提取 Y=0 截面上的有序四边形，通过共有边寻找下一层，不假定体单元排列。"""
    spec = read(prepared/"prepared.json")
    assert spec["caseSHA256"] == sha(CASE) and spec["geometry"]["outerRadiusD"] == 14
    for name in ("mesh.su2", "mapping.npz"):
        assert sha(prepared/name) == spec["files"][name]
    cells, xyz, _ = load_mesh(prepared/"mesh.su2", ("wall", "far"))
    with np.load(prepared/"mapping.npz") as mapping:
        np.testing.assert_array_equal(xyz, mapping["xyz"])
        wet, first = mapping["wet"], mapping["first"]
        corners, groups = mapping["corners"], mapping["group"]
    sphere = set(np.unique(corners[groups == 1]))
    tubes = set(np.unique(corners[np.isin(groups, (2, 3))]))
    joins = xyz[sorted(sphere & tubes)]
    assert len(joins) and len(xyz) == spec["geometry"]["nodes"]
    hexes = np.asarray([c for c in cells if len(c)==8], dtype=np.int32)
    in_plane = (abs(xyz[:,1]) < 1e-12) & (np.linalg.norm(xyz,axis=1) < .55)
    hexes = hexes[in_plane[hexes].sum(axis=1) >= 4]
    adjacent = defaultdict(dict)
    for indices in SIDES[8]:
        faces = hexes[:,indices]
        for row in faces[in_plane[faces].all(axis=1)]:
            face_key = tuple(sorted(row))
            for left, right, next_right, next_left in zip(row,np.roll(row,-1),np.roll(row,-2),np.roll(row,-3)):
                adjacent[(left,right)][face_key] = (next_left,next_right)
                adjacent[(right,left)][face_key] = (next_right,next_left)
    return dict(spec=spec, xyz=xyz, adjacent=adjacent, first=dict(zip(wet,first)), sphere=sphere, joins=joins)


def derivative(d, t, i, j):
    """以壁点及两个非等距流体点的拉格朗日二次多项式求壁面一阶导数。"""
    a, b = d[i], d[j]
    return (b*(t[i]-t[0])/a - a*(t[j]-t[0])/b)/(b-a)


def profile(grid, theta, station, physical):
    count = grid["spec"]["geometry"].get("nearWallIntervals",5)
    angle = np.deg2rad(station["angleDeg"])
    normal = np.array([np.cos(angle),0.,np.sin(angle)])
    tangent = np.array([-np.sin(angle),0.,np.cos(angle)])
    edge = (int(station["leftNode"]), int(station["rightNode"]))
    if not set(edge) <= grid["sphere"] or np.any(abs(grid["xyz"][list(edge),1]) > 1e-12):
        raise Unavailable("测站的起始边不在外露球面测线上")
    first_edge = tuple(grid["first"][i] for i in edge)
    previous = None
    samples = []
    for layer in range(count+1):
        if layer == 1 and edge != first_edge:
            raise Unavailable("体网邻接与独立的首层映射不符")
        points = grid["xyz"][list(edge)]
        delta = points[1]-points[0]
        denominator = delta @ tangent
        if abs(denominator) < 1e-12:
            raise Unavailable("测量法线与网格边平行")
        weight = -(points[0] @ tangent)/denominator
        location = points[0] + weight*delta
        if not (-1e-12 <= weight <= 1+1e-12 and location @ normal > 0):
            raise Unavailable("测量法线超出本层连接边")
        if layer == 0:
            np.testing.assert_allclose(weight,station["rightWeight"],rtol=0,atol=1e-12)
            wall = location
            dihedral = np.rad2deg(np.arcsin(abs(delta @ normal)/np.linalg.norm(delta)))
        distance = (location-wall) @ normal * physical["diameterM"]
        temps = physical["airK"] + physical["scalarUnitK"]*theta[list(edge)]
        samples.append(dict(ids=list(map(int,edge)),xyzD=points.tolist(),nodeK=temps.tolist(),weight=weight,
                            distanceM=distance,temperatureK=temps[0]+weight*(temps[1]-temps[0])))
        if layer < count:
            options=[(f,pair) for f,pair in grid["adjacent"][edge].items() if f!=previous]
            if len(options)!=1:
                raise Unavailable("向流体一侧的相邻层不唯一")
            previous,edge=options[0]
    distances=np.array([s["distanceM"] for s in samples])
    temps=np.array([s["temperatureK"] for s in samples])
    if not (np.diff(distances)>0).all():
        raise Unavailable("法向距离没有逐层增加")
    np.testing.assert_allclose(temps[0],station["wallK"],rtol=0,atol=1e-10)
    q_scale=physical["kFilm"]
    nu_scale=-physical["diameterM"]/((temps[0]-physical["airK"])*conductivity(temps[0]))
    reaction=station["fluxInWm2"]
    np.testing.assert_allclose(reaction*nu_scale,station["computedNu"],rtol=1e-12)
    linear=q_scale*(temps[1]-temps[0])/distances[1]
    q=q_scale*derivative(distances,temps,1,2)
    wide=q_scale*derivative(distances,temps,2,4)
    return dict(status="available",samples=samples,layerCount=count,firstDistanceM=distances[1],
                wallK=temps[0],sectionNormalAngleDeg=dihedral,
                joinPlaneClearanceM=np.min(abs(grid["joins"][:,1]))*physical["diameterM"],
                gradientInWm2=q,linearInWm2=linear,wideInWm2=wide,gradientNu=q*nu_scale,
                linearNu=linear*nu_scale,wideNu=wide*nu_scale)


def reconstruct(grid, frame, archive):
    spec=grid["spec"]
    with ZipFile(archive) as z:
        summary=json.loads(z.read("summary.json"))
        raw=z.read("stations.csv")
    mesh=summary["meshSHA256"] if "energy" in summary else summary["thermal"]["meshSHA256"]
    assert spec["files"]["mesh.su2"]==mesh and summary["condition"]==spec["condition"]
    assert sha(frame)==summary["frameSHA256"] and summary["caseSHA256"]==sha(CASE)
    assert hashlib.sha256(raw).hexdigest()==summary["files"]["stations.csv"]
    stations=np.genfromtxt(io.BytesIO(raw),names=True,delimiter=",")
    case=read(CASE)
    np.testing.assert_array_equal(stations["angleDeg"],case["angleDeg"])
    np.testing.assert_array_equal(stations["reportedNu"],case["localNuSurface"])
    field=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=(0,1,2,3,9))
    np.testing.assert_array_equal(field[:,0],np.arange(len(grid["xyz"])))
    np.testing.assert_array_equal(field[:,1:4],grid["xyz"])
    assert np.isfinite(field[:,4]).all()
    rows=[]
    for st,wire in zip(stations,case["wireGradientFPerFt"]):
        row=dict(angleDeg=st["angleDeg"],faceID=int(st["originalFaceID"]),experimentNu=st["reportedNu"],
                 wireGradientMagnitudeKPerM=wire*(5/9)/.3048,reactionInWm2=st["fluxInWm2"],reactionNu=st["computedNu"])
        try:
            row.update(profile(grid,field[:,4],st,spec["condition"]))
        except Unavailable as error:
            row.update(status="unavailable",reason=str(error),samples=[])
        rows.append(row)
    return dict(meshSHA256=mesh,frameSHA256=summary["frameSHA256"],caseSHA256=sha(CASE),
                stationsSHA256=summary["files"]["stations.csv"],elapsedS=summary.get("elapsedSeconds"),
                kFilm=spec["condition"]["kFilm"],rows=rows)


def save_pair(prepared,work,out):
    """在 14D 生产入口从本次真实终帧附加诊断，独立于守恒主面表。"""
    grid=measurement_grid(prepared)
    sets=[]
    for branch,step in (("full",1),("half",2)):
        frame=work/branch/f"step{step}"/f"restart_{step:05d}.csv"
        sets.append(reconstruct(grid,frame,out/branch/"surface.zip"))
    file=out/"stationGradients.json"
    with file.open("x",encoding="utf8") as f:
        json.dump(dict(schema="short65-station-gradient-v1",sets=sets),f,ensure_ascii=False,allow_nan=False)
    return sha(file)
