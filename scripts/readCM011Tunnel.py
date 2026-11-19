"""回读同一次 CM-011 低速计算，核对终态、完整原面和有符号的净热功率。"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from auditSU2SSTState import vtk_tokens
from cm011Tunnel import CP, D, DT, GRIDS, P0, PE, RE, RHO, T0, U, meridian, surface, volume, write_json
from runEmbeddedSU2 import sha
from su2AxisymmetricAnalysis import integrate_scalar


def match_nodes(actual, expected):
    delta, ids = cKDTree(actual[:, :2]).query(expected)
    if len(actual) != len(expected) or max(delta) > 2e-5 or len(np.unique(ids)) != len(ids):
        raise ValueError("VTK 节点与本次计算的子午网格不一致")
    return ids


def shell_balance(xy, velocity, theta, j):
    # 以终场内侧数据重构闭合旋转面的质量、显热对流和扩散。
    # 这项独立面积分不等于引擎在远场面上实际使用的数值通量。
    points = xy[j]
    edge = np.diff(points, axis=0)
    length = np.linalg.norm(edge, axis=1)
    normal = np.column_stack((-edge[:, 1], edge[:, 0]))/length[:, None]
    lo, hi = j-1, min(j+1, len(xy)-1)
    radial = xy[hi]-xy[lo]
    system = np.stack((edge, (radial[:-1]+radial[1:])*.5), axis=1)
    delta = theta[hi]-theta[lo]
    gradient = np.linalg.solve(system, np.column_stack((np.diff(theta[j]), (delta[:-1]+delta[1:])*.5))[..., None])[..., 0]
    d_normal = np.sum(gradient*normal, axis=1)
    mass, heat, inflow = 0., 0., 0.
    for f in (.5-.5/np.sqrt(3), .5+.5/np.sqrt(3)):
        area = np.pi*length*((1-f)*points[:-1, 1]+f*points[1:, 1])
        speed = np.sum(((1-f)*velocity[j, :-1]+f*velocity[j, 1:])*normal, axis=1)
        excess = (1-f)*theta[j, :-1]+f*theta[j, 1:]-1
        mass += np.sum(area*speed)
        inflow += np.sum(area*np.maximum(-speed, 0))
        heat += np.sum(area*(excess*speed-d_normal/PE))
    return dict(radialRow=int(j), outwardHeatW=heat*RHO*CP*U*DT*D**2,
                relativeNetMassFlow=mass/inflow)


def evaluate(folder):
    read = lambda name: json.loads((folder/name).read_bytes())
    receipt = read("engine_process.json")
    if (receipt["status"] != "terminal" or receipt["exitCode"] != 0 or receipt["timedOut"] or
            receipt["cancelled"] or receipt.get("supervisorError") or not all(receipt["inputsStable"].values())):
        raise ValueError("引擎缺少可读的正常终态")
    for name in ("body.json", "case.cfg", "mesh.su2"):
        if sha(folder/name) != receipt["inputSHA256"][str(folder/name)]:
            raise ValueError("运行输入已更改："+name)
    for name in ("history.csv", "surface.vtk", "vol_solution.vtk", "restart.csv"):
        if sha(folder/name) != receipt["outputSHA256"][name]:
            raise ValueError("运行结果已更改："+name)
    body = read("body.json")
    level, layers, extent = GRIDS[body["grid"]]
    w, bands = meridian(level)
    assert np.array_equal(w, np.asarray(body["meridian"])) and body["diameterM"] == D
    reference = surface(w, bands, level)
    for key, field in reference.items():
        assert np.array_equal(field, body[key]), "试件原面与母线不匹配："+key
    config = dict(line.split("=", 1) for line in (folder/"case.cfg").read_text("ascii").splitlines() if "=" in line)
    physics = {"SOLVER": "INC_NAVIER_STOKES", "AXISYMMETRIC": "YES", "INC_ENERGY_EQUATION": "YES",
               "INC_DENSITY_MODEL": "CONSTANT", "INC_NONDIM": "DIMENSIONAL", "INC_DENSITY_INIT": "1",
               "INC_VELOCITY_INIT": "(-1, 0, 0)", "INC_TEMPERATURE_INIT": "1", "SPECIFIC_HEAT_CP": "1",
               "MARKER_ISOTHERMAL": "(wall, 2)", "MARKER_FAR": "(far)", "MARKER_SYM": "(axis)",
               "MU_CONSTANT": f"{1/RE:.17g}", "THERMAL_CONDUCTIVITY_CONSTANT": f"{1/PE:.17g}"}
    assert all(config[key].strip() == value for key, value in physics.items()), "保存算例的工况或物性不符"
    assert receipt["executableVersion"] == "8.6.0" and receipt["runtimeFlavor"] == "official"
    assert receipt["inputSHA256"][receipt["command"][0]] == "01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
    xy, _, _ = volume(w, layers, extent, level)
    v, fluid = vtk_tokens(folder/"vol_solution.vtk")
    order = match_nodes(v, xy)
    xy = xy.reshape(layers+1, len(w), 2)
    temp = fluid["Temperature"].ravel()[order].reshape(layers+1, len(w))
    vel = fluid["Velocity"][order, :2].reshape(xy.shape)
    s, wall = vtk_tokens(folder/"surface.vtk")
    ids = match_nodes(s, w)
    q = wall["Heat_Flux"].ravel()[ids]*RHO*CP*U*DT
    p = P0 + wall["Pressure"].ravel()[ids]*RHO*U**2
    observed = wall["Temperature"].ravel()[ids]
    assert np.allclose(observed, temp[0], atol=1e-6, rtol=0)
    edges = np.column_stack((np.arange(len(w)-1), np.arange(1, len(w))))
    curved_power = integrate_scalar(w*D, edges, q)
    v = np.asarray(body["V"])
    triangles = np.asarray(body["F"])-1
    source_rows = np.asarray(body["nodeRow"])[triangles]
    def to_faces(field):
        return np.mean(field[source_rows], axis=1)
    area = .5*np.linalg.norm(np.cross(v[triangles[:, 1]]-v[triangles[:, 0]], v[triangles[:, 2]]-v[triangles[:, 0]]), axis=1)
    qf = to_faces(q)
    with (folder/"history.csv").open(encoding="utf8") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        rows = list(reader)
    residuals = [max(float(r[k]) for r in rows[-20:]) for k in ("rms[P]", "rms[U]", "rms[V]", "rms[h]")]
    monitor = np.array([float(r["HF"]) for r in rows[-20:]])
    heat_change = np.ptp(monitor)/max(abs(monitor[-1]), np.finfo(float).tiny)
    stop_index = int(rows[-1]["Inner_Iter"])
    assert sha(folder/"surface.vtk") == sha(folder/f"surface_{stop_index:06d}.vtk")
    assert sha(folder/"restart.csv") == sha(folder/f"restart_{stop_index:06d}.csv")
    prefix = [sha(folder/"surface.vtk"), sha(folder/"vol_solution.vtk")]
    balances = [shell_balance(xy, vel, temp, j) for j in (layers//2, layers-1, layers)]
    for b in balances:
        b["heatDiscrepancyFraction"] = (b["outwardHeatW"]+curved_power)/abs(curved_power)
    accepted = (len(rows) >= 20 and max(residuals) < -8 and heat_change < .002 and
                np.all(np.isfinite(q)) and np.all(q < 0) and np.all(p > 0) and
                np.all(np.isfinite(temp)) and np.all(temp[1] < 2) and np.all(temp[1] > 1))
    groups = {}
    for i, name in enumerate(body["regionNames"]):
        select = np.asarray(body["regionCode"]) == i
        groups[name] = dict(faces=int(select.sum()), areaM2=float(area[select].sum()),
                            powerW=float(qf[select]@area[select]), fluxRangeWm2=[min(qf[select]), max(qf[select])])
    stats = dict(grid=body["grid"], acceptedNumericalWall=bool(accepted), physicalAccuracyValidated=False,
                 originalFaces=len(area), meridionalWallNodes=len(w), volumeNodes=len(order),
                 lastIteration=stop_index, nativeUpdatesInRun=len(rows), runtimeSeconds=receipt["elapsedSeconds"],
                 last20MaximumLogRMS=residuals, last20HeatMonitorRelativeRange=heat_change,
                 areaM2=area.sum(), curvedAreaM2=integrate_scalar(w*D, edges, np.ones(len(w))),
                 triangularPowerW=qf@area, curvedPowerW=curved_power, byRegion=groups,
                 pressureRangePa=[min(to_faces(p)), max(to_faces(p))], fluxRangeWm2=[min(qf), max(qf)],
                 fluidTemperatureRangeK=[T0+DT*(temp.min()-1), T0+DT*(temp.max()-1)],
                 wallNodeTemperatureRangeK=[T0+DT*(observed.min()-1), T0+DT*(observed.max()-1)],
                 firstFluidNodeTemperatureRangeK=[T0+DT*(temp[1].min()-1), T0+DT*(temp[1].max()-1)],
                 shellReconstruction=balances, sourceVTKSHA256=prefix, reynoldsD=RE, pecletD=PE)
    wall_rows = np.column_stack((np.arange(len(area)), to_faces(p), qf, T0+DT*(to_faces(temp[1])-1),
                                 T0+DT*(to_faces(observed)-1)))
    return stats, wall_rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("new_output", type=Path)
    args = parser.parse_args()
    stats, wall_rows = evaluate(args.folder.resolve())
    args.new_output.mkdir(parents=True, exist_ok=False)
    write_json(args.new_output/"summary.json", stats)
    if stats["acceptedNumericalWall"]:
        np.savetxt(args.new_output/"wall.csv", wall_rows, delimiter=",", comments="",
                   header="faceID,pressurePa,heatFluxWm2,nearWallCellTemperatureK,wallFluidTemperatureK",
                   fmt=["%d", "%.17g", "%.17g", "%.17g", "%.17g"])
