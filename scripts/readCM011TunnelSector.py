"""将固定低速试件的三维扇区原生帧按镜像对称映射到全部原三角面。"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from auditSU2SSTState import vtk_tokens
from cm011Tunnel import CP, D, DT, P0, RHO, T0, U, meridian, write_json
from cm011TunnelSector import ANGLE, RAYS, REFERENCE
from readCM011Tunnel import evaluate
from runEmbeddedSU2 import sha


def checked_process(folder):
    rec = json.loads((folder/"engine_process.json").read_bytes())
    if (rec["status"] != "terminal" or rec.get("supervisorError") or rec["cancelled"] or
            not all(rec["inputsStable"].values()) or rec["runtimeFlavor"] != "official" or
            rec["executableVersion"] != "8.6.0"):
        raise ValueError("三维运行的来源或终态未确认")
    if rec["exitCode"] != 0 and not rec["timedOut"]:
        raise ValueError("求解器异常退出")
    for n in ("mesh.su2","case.cfg","initial_flow.csv","geometry.json","nodes.npz"):
        assert sha(folder/n)==rec["inputSHA256"][str(folder/n)], n
    assert rec["inputSHA256"][rec["command"][0]]=="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
    return rec


def state_at_frame(folder):
    rec = checked_process(folder)
    # 只接收本次固定工况、几何及完整初场的三维对照。
    fixed = {"mesh.su2":"f315415421c134e6609b88030190f1ecb4f9053025e7e6e8e63b77480fb2a97a",
             "case.cfg":"d07e10ea75b335b9abbfc1af380c87e684078d71915e64aae8323b7682a68da6",
             "initial_flow.csv":"94569fba64f67299f4b9130f55a80d44a1e086972767a4ca82aef7662bc53e69"}
    parent_index = None
    if sha(folder/"case.cfg")=="1b36120141e4e7eff9c3a39ae69f5a7f0253f85f140836585a3fac6a9701d6c6":
        # 接续读数的本地索引自零计，来源是另一进程的完整 13000 帧。
        parent_index = 13000
        fixed.update({"case.cfg":"1b36120141e4e7eff9c3a39ae69f5a7f0253f85f140836585a3fac6a9701d6c6",
                      "initial_flow.csv":"38acf6c69405b7c81201e4905bf90a1a4d45c31ac54c71a5b6638e4375706c75"})
        assert "1461b47ec49cdce1e6ec0ef169c561e9c0a96e901d594b47254c345b2eb9570a" in rec["inputSHA256"].values()
    assert all(sha(folder/n)==digest for n,digest in fixed.items())
    return frame_window(folder,rec,parent_index)


def frame_window(folder,rec,parent_index=None):
    assert sha(folder/"history.csv")==rec["outputSHA256"]["history.csv"]
    with (folder/"history.csv").open() as f:
        reader = csv.DictReader(f,skipinitialspace=True)
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        history = list(reader)
    last = int(history[-1]["Inner_Iter"])
    frames = sorted(int(p.stem.split("_")[-1]) for p in folder.glob("restart_*.csv")
                    if (folder/("surface_"+p.stem.split("_")[-1]+".vtk")).is_file())
    if not frames or frames[-1]>last:
        raise ValueError("缺少配对的原生保存帧")
    idx = frames[-1]
    normal = rec["exitCode"]==0 and not rec["timedOut"]
    names = [f"restart_{idx:06d}.csv",f"surface_{idx:06d}.vtk"]
    if normal:
        assert idx==last and all(sha(folder/n)==sha(folder/m) for n,m in zip(names,("restart.csv","surface.vtk")))
    for n in names: assert sha(folder/n)==rec["outputSHA256"][n]
    window = [row for row in history if idx-20 < int(row["Inner_Iter"]) <= idx]
    assert len(window)==20 and [int(r["Inner_Iter"]) for r in window]==list(range(idx-19,idx+1))
    res = [max(float(r[key]) for r in window) for key in ("rms[P]","rms[U]","rms[V]","rms[W]","rms[h]")]
    hf = np.array([float(r["HF"]) for r in window])
    spread = np.ptp(hf)/abs(hf[-1])
    return dict(normalTerminal=normal,timedOut=rec["timedOut"],savedIndex=idx,lastHistoryIndex=last,
                parentSavedIndex=parent_index,
                frameSHA256={n:sha(folder/n) for n in names},last20MaximumLogRMS=res,
                last20HeatRelativeRange=spread, historySectorWallMonitor=float(hf[-1]),
                passesResidualAndHeatWindow=bool(max(res)<-8 and spread<.002))


def wall_fields(folder,idx,m,rays):
    """独立核对三维原生坐标、壁压和焓，并从同一保存帧读取热流及首层温度。"""
    dims = np.load(folder/"nodes.npz")
    xyz,ids = dims["xyz"],dims["ids"]
    restart = np.loadtxt(folder/f"restart_{idx:06d}.csv",delimiter=",",skiprows=1)
    assert np.isfinite(restart).all() and np.array_equal(restart[:,0],np.arange(len(xyz)))
    np.testing.assert_allclose(restart[:,1:4],xyz,rtol=0,atol=1e-12)
    wall_ids = ids[:,:m]
    surf,fields = vtk_tokens(folder/f"surface_{idx:06d}.vtk")
    unique = np.unique(wall_ids)
    dist,order = cKDTree(surf).query(xyz[unique])
    assert len(unique)==len(surf)==len(np.unique(order)) and max(dist)<2e-6
    order = order[np.searchsorted(unique,wall_ids)]
    pressure = fields["Pressure"][:,0][order]*RHO*U**2
    heat = fields["Heat_Flux"][:,0][order]*RHO*CP*U*DT
    temp = fields["Temperature"][:,0][order]
    np.testing.assert_allclose(temp,restart[wall_ids,8]+298.15,rtol=0,atol=2e-6)
    np.testing.assert_allclose(pressure,restart[wall_ids,4]*RHO*U**2,rtol=0,atol=2e-8)
    normal_speeds = [np.max(abs(restart[ids[k],7]*np.cos(rays[k])-restart[ids[k],6]*np.sin(rays[k]))) for k in (0,-1)]
    no_slip = np.max(abs(restart[wall_ids,5:8]))
    near_temp = restart[ids[:,m:2*m],8]+298.15
    assert np.isfinite(heat).all() and (pressure+P0>0).all()
    # 在求解器的弦面上用分片线性热流作独立三角积分，再与原三角面映射分别比较。
    chord_area, chord_power = 0., 0.
    for k in range(len(rays)-1):
        coords = xyz[wall_ids[k:k+2]]
        for strip in (((0,0),(0,1),(1,1)),((0,0),(1,1),(1,0))):
            pts = [coords[ray,start:start+m-1] for ray,start in strip]
            areas = .5*np.linalg.norm(np.cross(pts[1]-pts[0],pts[2]-pts[0]),axis=1)*D**2
            flux = np.mean([heat[k+ray,start:start+m-1] for ray,start in strip],axis=0)
            chord_area += np.sum(areas)
            chord_power += areas@flux
    return pressure,heat,temp,near_temp,normal_speeds,no_slip,chord_area,chord_power,restart


def map_fields(folder, parent):
    refs = json.loads(REFERENCE.read_bytes())
    if any(sha(parent/n)!=h for n,h in refs["sourceSHA256"].items()):
        raise ValueError("二维对照和当前原面存档不匹配")
    baseline,wall2 = evaluate(parent)
    body = json.loads((parent/"body.json").read_bytes())
    geom = json.loads((folder/"geometry.json").read_bytes())
    if sha(parent/"body.json")!=geom["parentBodySHA256"] or sha(parent/"restart.csv")!=geom["parentRestartSHA256"]:
        raise ValueError("三维种子与原面存档不匹配")
    report = state_at_frame(folder)
    wall,_ = meridian(2)
    pressure,heat,temp,near_temp,normal_speeds,no_slip,chord_area,chord_power,restart = wall_fields(
        folder,report["savedIndex"],len(wall),RAYS)
    assert np.isclose(chord_area,geom["chordWallAreaM2"],rtol=1e-12,atol=0)

    # 每两个扇区互为镜像。顶点先对应母线行，再沿周向在相邻求解射线上插值。
    vertices = np.array(body["V"])/D
    faces = np.asarray(body["F"])-1
    row = np.asarray(body["nodeRow"])
    np.testing.assert_allclose(np.column_stack((vertices[:,0],np.hypot(vertices[:,1],vertices[:,2]))),
                               wall[row],atol=1e-12,rtol=0)
    phase = np.mod(np.arctan2(vertices[:,2],vertices[:,1]),2*ANGLE)
    fraction = np.minimum(phase,2*ANGLE-phase)/(ANGLE/2)
    left = np.minimum(np.floor(fraction).astype(int),1)
    weight = fraction-left
    assert np.all((weight>=0)&(weight<=1+1e-13))
    def to_faces(field):
        vertex = (1-weight)*field[left,row]+weight*field[left+1,row]
        return vertex[faces].mean(axis=1)
    gauge,q,tw,tnear = [to_faces(field) for field in (pressure,heat,T0+DT*(temp-1),T0+DT*(near_temp-1))]
    p2,q2 = wall2[:,1]-P0,wall2[:,2]
    v = np.asarray(body["V"])[faces]
    area = .5*np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)
    accepted = (report["normalTerminal"] and report["passesResidualAndHeatWindow"] and
                max(normal_speeds)<1e-8 and no_slip<1e-8 and np.all(heat<0) and np.all((near_temp>1)&(near_temp<2)))
    report.update(numericalWallAccepted=bool(accepted), physicalAccuracyValidated=False,
                  originalFaces=len(q), totalAreaM2=area.sum(), maxSymmetryNormalSpeedOverU=normal_speeds,
                  maxWallSpeedOverU=no_slip, fullTrianglePowerW=q@area, axisymmetricTrianglePowerW=q2@area,
                  fullChordTrianglePowerW=60*chord_power,
                  fullWallMonitorPowerW=60*report["historySectorWallMonitor"]*RHO*CP*U*DT*D**2,
                  wallFluidTemperatureRangeK=[tw.min(),tw.max()], firstFluidNodeTemperatureRangeK=[tnear.min(),tnear.max()],
                  volumeFluidTemperatureRangeK=[T0+DT*(x+297.15) for x in (restart[:,8].min(),restart[:,8].max())],
                  heatFluxRangeWm2=[q.min(),q.max()], gaugePressureRangePa=[gauge.min(),gauge.max()],
                  maxAbsFaceHeatDifferenceWm2=max(abs(q-q2)),maxAbsFaceGaugeDifferencePa=max(abs(gauge-p2)),
                  byRegion=[],geometry=geom,bodySHA256=sha(parent/"body.json"))
    for i,name in enumerate(body["regionNames"]):
        select = np.asarray(body["regionCode"])==i
        report["byRegion"].append(dict(name=name,faces=int(select.sum()),powerW=q[select]@area[select],
             referencePowerW=q2[select]@area[select],areaM2=area[select].sum(),
             meanGaugePressurePa=np.average(gauge[select],weights=area[select]),
             referenceMeanGaugePa=np.average(p2[select],weights=area[select]),
             maxAbsHeatDifferenceWm2=max(abs(q[select]-q2[select]))))
    assert abs(report["totalAreaM2"]-baseline["areaM2"])<1e-16
    return report,np.column_stack((wall2[:,0],gauge,q,tnear,tw))


def read_energy(path, report, copies=60, symmetry=("sym0","sym1"), schema="cm011-sector-energy-v1"):
    """只把同一保存态通过复核的三维原生能量账附入对照字段。"""
    energy = json.loads(path.read_bytes())
    idx = report["savedIndex"]
    if (energy["schema"]!=schema or not report["numericalWallAccepted"] or
            energy["sourceFrameSHA256"]!=report["frameSHA256"] or
            energy["sourceRestartSHA256"]!=report["frameSHA256"][f"restart_{idx:06d}.csv"] or
            energy["sampledInputExact"] is not True or not 0<=energy["updatedStateDifference"]<2e-7 or
            energy["physicallyValidated"] is not False or energy["savedIndex"]!=idx or energy["copies"]!=copies):
        raise ValueError("三维观测的来源或阶段与当前帧不符")
    increments = energy["allBoundaryIncrements"]
    href = energy["referenceNativeEnthalpy"]
    for row in increments.values():
        np.testing.assert_allclose(row["referenceAdvectionW"],row["massKgS"]*CP*DT*href,rtol=0,atol=1e-12)
        np.testing.assert_allclose(row["excessHeatW"],row["rawEnthalpyW"]-row["referenceAdvectionW"],rtol=0,atol=1e-12)
    wall = -increments["after_wall"]["excessHeatW"]
    far = increments["after_far"]["excessHeatW"]
    assert np.isfinite([wall,far,href]).all() and min(wall,far)>0
    gap = wall-far
    expected = [wall,far,gap,gap/wall,-gap]
    keys = ["nativeWallReleaseW","nativeFarExcessOutW","wallMinusFarW","wallMinusFarFraction","weightedResidualW"]
    np.testing.assert_allclose([energy[key] for key in keys],expected,rtol=0,atol=1e-12)
    for key in ("start","internal","before_far","before_wall",*symmetry,"end"):
        assert abs(increments[key]["excessHeatW"])<1e-12 and abs(increments[key]["massKgS"])<1e-12
    assert energy["belowHalfPercent"]==bool(abs(gap)<.005*wall)
    return energy


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("sector","parent","new_output"): parser.add_argument(arg,type=Path)
    parser.add_argument("--energy-observation",type=Path)
    args=parser.parse_args()
    stats, rows=map_fields(args.sector.resolve(),args.parent.resolve())
    if args.energy_observation:
        stats["nativeEnergyBalance"]=read_energy(args.energy_observation,stats)
    args.new_output.mkdir(parents=True,exist_ok=False)
    write_json(args.new_output/"summary.json",stats)
    np.savetxt(args.new_output/"wall.csv",rows,header="faceID,gaugePressurePa,heatFluxWm2,firstNodeK,wallNodeK",
               delimiter=",",comments="",fmt=["%d"]+["%.17g"]*4)
