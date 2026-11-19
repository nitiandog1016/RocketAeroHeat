"""用独立的三角形面积、原生壁面表和源几何复核 CM-011 工作台并绘制母线分布。"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens


def read(path):
    return json.loads(path.read_bytes())


def reference_geometry():
    # 图纸外径和球半径在这里直接用英寸。母线分段按解析面积与体积积分。
    diam = 155.481
    scale = .0254/1000
    f, t, s = 180.890, 9.0118, (7.46+7.17)/2
    angle = math.radians(32+31/60+23/3600)
    cr = diam/2-s
    stop = math.asin(cr/(f-s))
    tangent = math.pi/2+angle
    pairs = [(t, 0., tangent, math.pi), (s, cr, stop, tangent), (f, 0., 0., stop)]
    cone_radii = np.array([t*math.cos(angle), cr+s*math.cos(angle)])
    delta_x = np.diff(cone_radii)[0]/math.tan(angle)
    area = math.pi*sum(cone_radii)*delta_x/math.cos(angle)
    volume = math.pi*delta_x*(sum(cone_radii**2)+np.prod(cone_radii))/3
    for radius, center, a, b in pairs:
        area += quad(lambda v: 2*math.pi*radius*(center+radius*math.sin(v)), a, b)[0]
        volume += quad(lambda v: math.pi*(center+radius*math.sin(v))**2*radius*math.sin(v), a, b)[0]
    return area*scale**2, volume*scale**3


def check_page(page, run):
    raw = (page/"case-data.js").read_text(encoding="utf8").strip()
    head = "window.ROCKET_CASE="
    assert raw.startswith(head) and raw.endswith(";")
    data = json.loads(raw[len(head):-1])
    body = read(run/"body.json")
    v, f = np.asarray(data["vertices"]), np.asarray(data["faces"])
    assert np.array_equal(v, body["V"]) and np.array_equal(f, np.asarray(body["F"])-1)
    edges = np.vstack((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]))
    _, group, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    assert np.all(counts == 2) and not np.any(np.bincount(group, weights=np.sign(edges[:, 1]-edges[:, 0])))
    assert len(v)-len(counts)+len(f) == 2
    normal = np.cross(v[f[:, 1]]-v[f[:, 0]], v[f[:, 2]]-v[f[:, 0]])
    area = np.linalg.norm(normal, axis=1)/2
    volume = np.sum(v[f[:, 0]]*normal)/6
    exact_a, exact_v = reference_geometry()
    assert volume > 0 and abs(area.sum()/exact_a-1) < .005 and abs(volume/exact_v-1) < .01
    assert abs(np.ptp(v[:, 1])-155.481*.0254/1000) < 1e-13
    fields = {item["key"]: item for item in data["fields"]}
    wall_xyz, src = vtk_tokens(run/"surface.vtk")
    w = np.asarray(body["meridian"])
    dist, ix = cKDTree(wall_xyz[:, :2]).query(w)
    assert max(dist)<1e-6 and len(np.unique(ix))==len(wall_xyz)
    face_rows = ix[np.asarray(body["nodeRow"])[f]]
    native_p = src["Pressure"].ravel()[face_rows].mean(axis=1)*1.177*.1**2
    native_q = src["Heat_Flux"].ravel()[face_rows].mean(axis=1)*1.177*1007*.1*.1
    checks = {"pressure": (101325+native_p)/1000, "pressureGauge": native_p, "heatFlux": native_q/1000,
              "wallFluidTemperature": 300+.1*(src["Temperature"].ravel()[face_rows].mean(axis=1)-1)}
    for key, expected in checks.items():
        fld = fields[key]
        assert all(i==1 for i in fld["status"]) and len(fld["methodNames"]) and np.allclose(fld["values"], expected, atol=2e-10, rtol=0)
        assert fld["computedStateLabel"]=="研究数值解" and "修正牛顿" not in fld["stateNote"]
    assert all(i==1 for i in fields["nearWallCellTemperature"]["status"])
    assert all(i==2 for i in fields["wallInput"]["status"]) and np.allclose(fields["wallInput"]["values"], 300.1)
    for key in ("edgeTemperature", "recoveryTemperature", "wallSolved"):
        assert all(i==0 for i in fields[key]["status"]) and all(x is None for x in fields[key]["values"])
    result = data["wallCoverage"]
    pwr = float(area@native_q)
    assert result["fullSurfaceComputed"] and result["jointComputedFaceCount"] == len(f)
    assert result["jointMissingAreaM2"] == 0 and abs(result["wholeSurfaceHeatPowerW"]/pwr-1)<1e-11
    assert data["evidenceIdentity"]["origin"]=="research"
    assert all(key in data["caseNote"] for key in ("CM-011", "100.00%", "300 K"))
    return dict(faces=len(f), areaM2=area.sum(), analyticAreaM2=exact_a, volumeM3=volume, analyticVolumeM3=exact_v,
                independentTrianglePowerW=pwr, wattUnit="positive into wall"), (w, src, ix)


def curves(run):
    body = read(run/"body.json")
    w = np.asarray(body["meridian"])
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(w, axis=0), axis=1))]
    frames = sorted(run.glob("surface_*.vtk"))[-3:]
    history = []
    for frame in frames:
        xyz, data = vtk_tokens(frame)
        order = cKDTree(xyz[:, :2]).query(w)[1]
        flux = 1.177*1007*.1*.1*data["Heat_Flux"].ravel()[order]
        p = 1.177*.1**2*data["Pressure"].ravel()[order]
        # 两点高斯公式对每个直线母线段的 r*q 乘积精确积分。
        power = 0.
        for z in (.5-.5/math.sqrt(3), .5+.5/math.sqrt(3)):
            power += np.pi*(155.481*.0254/1000)**2*np.sum(
                np.diff(arc)*(w[:-1, 1]*(1-z)+w[1:, 1]*z)*(flux[:-1]*(1-z)+flux[1:]*z))
        history.append([int(frame.stem.split("_")[-1]), power])
    assert len(history)==3 and np.ptp(np.asarray(history)[:, 1])/abs(power)<.002
    return (arc[-1]-arc)/arc[-1], p, flux, history


def main(base, out):
    out.mkdir(parents=True, exist_ok=False)
    checks = {}
    for view, run in [("entry01", "entry01/flow"), ("fine_view01", "fine02")]:
        checks[view], _ = check_page(base/view/"workbench_final", base/run)
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 5.1), sharex=True, constrained_layout=True)
    steps = {}
    for run, name, style in [("coarse03","Coarse, 12D","--"), ("outer02","Fine, 20D",":"),
                             ("fine02","Fine, 12D","-")]:
        x, p, q, steps[run] = curves(base/run)
        axes[0].plot(x, p, style, label=name, lw=1.5)
        axes[1].plot(x, q, style, lw=1.5)
    axes[0].set(ylabel="Pressure above ambient [Pa]", title="CM-011 laboratory model   |   0.1 m/s air, 300.1 K wall")
    axes[0].legend(ncol=3, fontsize=8)
    axes[1].set(xlabel="Meridional distance from shield pole / total arc length", ylabel="Heat flux into wall [W/m²]", xlim=(0, 1))
    for ax in axes: ax.grid(alpha=.23)
    fig.savefig(out/"wall_distribution.png", dpi=160)
    plt.close(fig)
    manifest = {}
    for n in ("coarse01","coarse02","coarse03","fine01","fine02","outer01","outer02","entry01/flow"):
        raw = base/n/"engine_process.json"; j=read(raw)
        manifest[n] = dict(receiptSHA256=hashlib.sha256(raw.read_bytes()).hexdigest(),
                           exitCode=j["exitCode"], timedOut=j["timedOut"], seconds=j["elapsedSeconds"],
                           originalInputHashesStable=all(j["inputsStable"].values()),
                           bytes=sum(p.stat().st_size for p in raw.parent.iterdir() if p.is_file()))
    with (out/"export_check.json").open("x", encoding="utf8") as stream:
        json.dump(dict(pageChecks=checks, nativeAxisymmetricPowerHistory=steps, runs=manifest,
                       physicallyValidatedAreaM2=0, customRocketCoverageGainM2=0, browserInteractionTested=False),
                  stream, indent=2, allow_nan=False)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("new_output", type=Path)
    args=parser.parse_args()
    main(args.base, args.new_output)
