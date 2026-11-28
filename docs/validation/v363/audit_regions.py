"""核对飞行档案、同形风洞面表与独立试件的来源，保留六区原面编号。"""
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FLIGHT = HERE.parent / "v358/evidence"
CASE = ROOT / "cases/simple_launcher"
sys.path.insert(0, str(FLIGHT.parent))
from verify_coverage import verify


# 针对 M0.30 的资料判断；气流状态由同工况证据决定，不能从区域名称推出。
SCOPES = {
    "sphericalNose": ("球形鼻端", "滑移球对照给出了理想球体的压力误差。鼻端各三角面覆盖有限面积，项目尚无该工况球冠的壁压或换热测量。"),
    "cone": ("锥形头罩", "项目没有该工况头罩的沿面流动记录。较高马赫数的工程壁压、热流均属各自工况。"),
    "cylinder": ("圆柱箭身", "等温平板试验在首条铜带促转捩。此工况圆柱段的转捩位置、附着状态和沿程压力均未确定。"),
    "contraction": ("收缩肩", "尚无该工况肩前至肩后的壁剪切或近壁流向记录，分离位置未知。"),
    "afterbody": ("尾段", "肩后流动会影响这段表面。该工况的附着、回流或再附位置均未确定。"),
    "base": ("平底", "同形风洞归档有四个半径的数值压力读数。M0.30 的底面压力、尾迹流态和外域尺寸误差均待核对。"),
}


def read(path):
    return json.loads(path.read_bytes())


def fingerprint(path):
    raw = path.read_bytes()
    if path.suffix in (".json",".md",".m"):
        raw = raw.replace(b"\r\n",b"\n")
    return hashlib.sha256(raw).hexdigest()


def audit():
    definition = read(CASE/"geometry.json")
    archive, masks = verify(FLIGHT)
    target, = (entry for entry in archive["cases"] if entry["mach"] == 0.30)
    assert not target["solverReturned"] and not target["wholeSurfacePowerKnown"]
    assert target["lastStage"] == "freestream"
    assert not masks["case_01_pressure"].any() and not masks["case_01_heat"].any()
    # 相邻原面状态相同才合并；MATLAB 拾取可由面号直接读到当时的计算/给定/缺项身份。
    face_states = []
    for index,entry in enumerate(archive["cases"],1):
        p = masks[f"case_{index:02}_pressure"]
        q = masks[f"case_{index:02}_heat"]
        last = np.r_[np.flatnonzero((p[1:]!=p[:-1]) | (q[1:]!=q[:-1]))+1,len(p)]
        lengths = np.diff(np.r_[0,last])
        np.testing.assert_array_equal(np.repeat(p[last-1],lengths),p)
        np.testing.assert_array_equal(np.repeat(q[last-1],lengths),q)
        face_states.append(dict(mach=entry["mach"],lastFace=last.tolist(),
                                pressureState=p[last-1].tolist(),heatComputed=q[last-1].tolist()))
    stl = CASE / definition["file"]
    dtype = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])
    tri = np.fromfile(stl, dtype=dtype, offset=84)["vertices"].astype(float)
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2

    # 同形的两套风洞对照逐行核对坐标几何的面积和面号，另查它们实际记录了哪些物理量。
    sources = [stl, CASE/"geometry.json", FLIGHT/"cases.json", FLIGHT/"status_masks.npz"]
    wind = []
    for z, condition, seal, seal_key in (
        (CASE/"wind_tunnel.zip", CASE/"tunnel.json", CASE/"receipt.json", "sha256"),
        (CASE/"wake/momentum_grid.zip", CASE/"wake/condition.json", HERE.parent/"v357/fields.json", "bundleSHA256"),
    ):
        bc = read(condition)
        assert fingerprint(z)==read(seal)[seal_key] and bc["surfaceSHA256"]==archive["surfaceSHA256"]
        with zipfile.ZipFile(z) as bundle:
            rows = np.genfromtxt(io.BytesIO(bundle.read("wall.csv")), delimiter=",", names=True)
        np.testing.assert_array_equal(rows["faceID"], np.arange(1,len(area)+1))
        np.testing.assert_allclose(rows["areaM2"],area,rtol=0,atol=1e-14)
        assert np.isfinite(rows["pressurePa"]).all() and np.isfinite(rows["heatFluxWm2"]).all()
        speed = float(np.linalg.norm(bc["velocityMps"]))
        assert speed!=target["speedMps"] and bc["wallTemperatureK"]!=target["wallGivenK"]
        wind.append(dict(condition=str(condition.relative_to(ROOT).as_posix()), speedMps=speed,
                         givenWallK=bc["wallTemperatureK"], wallColumns=list(rows.dtype.names),
                         originalFaces=len(rows), elapsedPhysicalTimeS=bc.get("elapsedPhysicalTimeS")))
        sources.extend((z,condition,seal))

    # 汇总实际 STL 面，不用圆柱、圆锥等解析面积替换离散分母。
    regions = []
    membership = np.zeros(len(area), dtype=int)
    by_key = {r["name"]:r for r in definition["regions"]}
    assert set(by_key)==set(SCOPES)
    domain = read(HERE/"wind_domain.json")
    assert domain["surfaceSHA256"]==archive["surfaceSHA256"] and domain["speedMps"]==wind[0]["speedMps"]
    sources.append(HERE/"wind_domain.json")
    for key,(label,note) in SCOPES.items():
        region = by_key[key]
        first,last = region["firstFace"],region["lastFace"]
        assert 1 <= first <= last <= len(area)
        ids = slice(first-1,last)
        membership[ids] += 1
        total = float(area[ids].sum())
        data = target["regions"][key]
        assert data["faces"]==last-first+1 and abs(data["areaM2"]-total)<1e-12
        assert abs(region["areaM2"]-total)<1e-12
        regions.append(dict(key=key,label=label,firstFace=first,lastFace=last,faces=data["faces"],
                            areaM2=total,percentOfBody=100*total/area.sum(),
                            xMinM=float(tri[ids,:,0].min()),xMaxM=float(tri[ids,:,0].max()),
                            attachmentEvidenceFaces=0,flowState="未判定",note=note,
                            windDomain=domain["regions"][key]))
    np.testing.assert_array_equal(membership,np.ones(len(area)))
    assert abs(sum(r["areaM2"] for r in regions)-archive["areaM2"])<1e-12
    # 参考文件绑定到它们自己的算例，明确适用范围；诊断用资料不授予飞行面元资格。
    sphere = HERE.parent/"v362/sphere.json"
    plate = HERE.parent/"v360/comparison.json"
    taps = CASE/"wind_tunnel_stations.json"
    assert read(sphere)["commonNodes"]["bernoulli"]["passesWallGates"]
    assert read(plate)["stationCount"]==22
    assert read(taps)["surfaceSHA256"]==archive["surfaceSHA256"]
    sources.extend((sphere,plate,taps,HERE.parent/"v351/RESULTS.md",FLIGHT.parent/"check_flight_cases.m"))
    return dict(surfaceSHA256=archive["surfaceSHA256"],faces=len(area),areaM2=float(area.sum()),
                mach=target["mach"],altitudeM=target["altitudeM"],speedMps=target["speedMps"],
                wallGivenK=target["wallGivenK"],matchingFlowStateEvidenceFaces=0,
                matchingPressureFaces=target["coverage"]["computedPressureFaces"],
                matchingHeatFaces=target["coverage"]["heatFaces"],wholeSurfacePowerW=None,
                regions=regions,faceStates=face_states,windTunnelSources=wind,
                sourceHashConvention="SHA-256; LF for text, exact bytes for binary files",
                sourceSHA256={p.relative_to(ROOT).as_posix():fingerprint(p) for p in sources})


if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,indent=2,allow_nan=False))
