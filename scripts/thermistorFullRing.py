"""生成 B 型热敏珠 6 m/s 的全周非定常首步输入。使用名义珠体原面及恒物性空气。"""
import argparse
from pathlib import Path
import shutil

import numpy as np
from cm011TunnelSector import mesh, write_volume
from thermistorBead import ROOT, D, RHO, CP, MU, K, wall, write_json, sha

SPEED = 6.
ANGLE_CELLS = 48
MERIDIAN, LAYERS = 64, 48
PARENT = ROOT / "outputs/v283_validation/strong01"
SOURCE = {
    "mesh.su2": "f59482ab24bcdc3d886fac3e442ff474a13c2d34c0596c822ebfa2e7cb81e723",
    "nodes.npz": "8d3911634283e2df8d4fe61d3ee50a89cba4503b2eb16721cee365c97de01629",
    "body.json": "5e66bb2f70f0774b8ebb87a77b76061fc84b89b9f379d54e90879e81628bf62f",
    "restart.csv": "5f7eef6507a851ac0e6bc5c832dd2207a92bab55f0b83ea2ddfd4ace644957ba",
    "case.cfg": "b61456e3fc30a5161d4ff74717fb100d1b8418ac2822e027d2975e447d33c6f7",
}


def config(reynolds, peclet):
    """沿用已核对的层流与强温差标量，按新风速设置输运系数，在完整域走一个实际时间步。"""
    replacement = dict(MU_CONSTANT=1/reynolds, THERMAL_CONDUCTIVITY_CONSTANT=1/peclet,
                       DIFFUSIVITY_CONSTANT=1/peclet, SOLUTION_FILENAME="seed",
                       SCREEN_WRT_FREQ_INNER=1, OUTPUT_WRT_FREQ="(1, 1)",
                       OUTPUT_FILES="(RESTART_ASCII, SURFACE_PARAVIEW_ASCII)",
                       SCREEN_OUTPUT="(TIME_ITER, INNER_ITER, RMS_PRESSURE, RMS_VELOCITY-X, "
                                     "RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY, RMS_SPECIES_0)")
    lines = []
    for line in (PARENT/"case.cfg").read_text("ascii").splitlines():
        key = line.partition("=")[0].strip()
        if key in ("ITER", "MARKER_SYM", "WRT_RESTART_OVERWRITE", "WRT_SURFACE_OVERWRITE"):
            continue
        lines.append(f"{key}= {replacement.pop(key)}" if key in replacement else line)
    assert not replacement
    return "\n".join(lines)+"""
TIME_DOMAIN= YES
TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER
TIME_STEP= 0.05
MAX_TIME= 1.0
RESTART_ITER= 1
TIME_ITER= 2
INNER_ITER= 180
HISTORY_WRT_FREQ_INNER= 1
WRT_RESTART_COMPACT= NO
"""


def prepare(out):
    """焊接周向接缝，保存无对称边界的原生体网、非轴对称初猜和 SI 工况。"""
    if any((out/name).exists() for name in ("mesh.su2", "case.cfg", "seed_00000.csv", "geometry.json")):
        raise FileExistsError("请使用新的计算目录")
    for name, digest in SOURCE.items():
        if sha(PARENT/name) != digest:
            raise ValueError(f"热敏珠初猜来源不符 {name}")
    rays = np.linspace(0, 2*np.pi, ANGLE_CELLS+1)
    xy, xyz, source, phi, ids, cells, markers, geom = mesh(rays, 2, LAYERS, wall(MERIDIAN), 20., D)
    assert set(markers) == {"wall", "far"} and np.array_equal(ids[0], ids[-1])
    assert len(cells) == MERIDIAN*LAYERS*ANGLE_CELLS and len(xyz) == (LAYERS+1)*(2+(MERIDIAN-1)*ANGLE_CELLS)
    write_volume(out, xyz, cells, markers)
    np.savez_compressed(out/"nodes.npz", xyz=xyz, ids=ids)
    shutil.copyfile(PARENT/"body.json", out/"body.json")

    # 初始流动和温差形状来自 3 m/s 的无量纲参照帧，横向速度随本射线旋转。
    # 扰动只施加在尾迹环带，沿圆周变号；壁、远场以及极轴节点保持原初猜。
    with np.load(PARENT/"nodes.npz") as nodes:
        mother = nodes["ids"][0]
        np.testing.assert_allclose(nodes["xyz"][mother,:2], xy, atol=1e-13, rtol=0)
    base = np.loadtxt(PARENT/"restart.csv", delimiter=",", skiprows=1)
    assert np.array_equal(base[:,0], np.arange(len(base)))
    np.testing.assert_allclose(base[mother,1:3], xy, atol=1e-13, rtol=0)
    template = base[mother[source]]
    x, radius = xy[source].T
    ellipse = ((x+1.7)/.8)**2 + ((radius-.55)/.45)**2
    envelope = np.maximum(1-ellipse,0.)**2
    du = 1e-3*envelope*np.cos(phi)
    assert np.max(du)>.0008 and abs(np.sum(du))<1e-12
    assert not np.any(du[ids[:, :MERIDIAN+1]]) and not np.any(du[ids[:, -(MERIDIAN+1):]])
    values = np.column_stack((np.arange(len(xyz)),xyz,template[:,4],template[:,5]+du,
                              template[:,6]*np.cos(phi),template[:,6]*np.sin(phi),template[:,8:10]))
    with (PARENT/"restart.csv").open(encoding="ascii") as stream:
        header = stream.readline().strip()
    np.savetxt(out/"seed_00000.csv",values,delimiter=",",fmt=["%d"]+["%.17g"]*9,header=header,comments="")
    reynolds = RHO*SPEED*D/MU
    peclet = RHO*CP*SPEED*D/K
    (out/"case.cfg").write_text(config(reynolds,peclet), encoding="ascii")
    geom.update(cellCount=len(cells), angularCells=ANGLE_CELLS, seamInternalFaces=MERIDIAN*LAYERS,
                physicalSurfaces=list(markers), speedMps=SPEED, airK=293.15, nominalWallK=303.15,
                airPa=101325., rho=RHO, mu=MU, cp=CP, k=K, transverseDiameterM=D,
                reynoldsD=reynolds, pecletD=peclet, dtSeconds=.05*D/SPEED,
                peakAxialPerturbationOverU=float(np.max(np.abs(du))),
                initialGuessFrameSHA256=SOURCE["restart.csv"], originalBodySHA256=SOURCE["body.json"])
    write_json(out/"geometry.json",geom)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    directory=parser.parse_args().output.resolve()
    directory.relative_to(ROOT/"outputs")
    directory.mkdir(parents=True,exist_ok=True)
    prepare(directory)
