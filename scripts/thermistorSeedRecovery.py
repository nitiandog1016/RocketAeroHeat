"""按已验收的物理启动过程重建 B 型热敏珠的 48 层旧态。只读取算例输入和本次生成物。"""
import os
from pathlib import Path
import shutil
import sys

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
import numpy as np

import thermistorBead as bead
import thermistorBeadSector as sector
import thermistorTemperature as strong
from runSU2NativeDebugger import observe
from thermistorNearWall import FIELDS
from thermistorStartupSurface import columns, write

CASE = bead.ROOT / "cases/thermistor_b/startup6"
COARSE_SEED = "ea94eb8a3b668a47e1d3e66e59065e704ab05c57ae752fcd6ab75fa112f42445"
BASE_FRAME = "a9286dea2d009bfe7d3a78040c6b84b5cc7e8de4b4460adea6d5a86c121b028d"
RELAX_FRAME = "37c9ce962b213684c307bffed39fd524887ce910b739a6a27a226fa6246ad6b1"


def identical(path, expected):
    actual = bead.sha(path)
    if actual != expected:
        raise ValueError(f"{path} 的 SHA-256 为 {actual}，登记值为 {expected}")
    return path


def write_seed(file, state):
    with file.open("x", encoding="ascii") as f:
        np.savetxt(f, state, fmt=["%d"] + ["%.17g"]*9, delimiter=",", comments="",
                   header=",".join('"' + key + '"' for key in FIELDS))


def cold_sector(work):
    """从冷态逐级取得 7.5° 扇区的强温差参考帧，与 v283 的实际保存帧逐字节核对。"""
    axis, weak, hot = (work / s for s in ("axis3", "sector3", "strong3"))
    bead.prepare(axis, "outer")
    bead.solve(axis)
    sector.prepare(axis, weak, 4)
    bead.solve(weak)
    strong.prepare(weak, hot)
    bead.solve(hot)
    identical(hot / "restart.csv", "5f7eef6507a851ac0e6bc5c832dd2207a92bab55f0b83ea2ddfd4ace644957ba")
    return hot


def native(folder, threads, seconds, mb):
    identical(bead.ENGINE, bead.ENGINE_SHA)
    locks = [bead.ENGINE, Path(__file__).resolve(), folder/"case.cfg", folder/"mesh.su2", folder/"initial_flow.csv"]
    rec = observe([str(bead.ENGINE),"case.cfg","--threads",str(threads)], folder,
                  {str(p):bead.sha(p) for p in locks}, seconds, mb*1024**2)
    assert rec["debuggerExitCode"] == 0 and not rec["stopReason"] and not rec["survivingProcessIDs"]
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")


def matched(work, hot):
    """使用同速扇区输运系数求解独立基流，用真实父态构造全周初始速度与温差。"""
    out = work / "base6"
    out.mkdir()
    for src, dest in (("mesh.su2","mesh.su2"), ("restart.csv","initial_flow.csv")):
        os.link(hot/src,out/dest)
    shutil.copyfile(identical(CASE/"base6.cfg","0b54fe154ad50209f8550ab5b6f469aa36c0a8d7078e62a497f4c9b4d9731045"),out/"case.cfg")
    native(out,1,600,180)
    return columns(identical(out/"restart.csv",BASE_FRAME), FIELDS)


def recover(work, coarse):
    hot = cold_sector(work)
    base = matched(work, hot)
    with np.load(hot/"nodes.npz") as z:
        mother = z["ids"]
    with np.load(coarse/"nodes.npz") as z:
        xyz, ids = z["xyz"], z["ids"]
    # 五条扇区射线交替反射。扰动保留原两次浮点加减的顺序，与既有初场逐点对应。
    nodes, first = np.unique(ids[:-1].ravel(),return_index=True)
    ray, station = np.divmod(first,ids.shape[1])
    folded = ray % 8
    k = np.minimum(folded,8-folded)
    n = base[mother[k,station]]
    phi = ray*np.pi/24
    angle = k*np.pi/24
    radial = n[:,6]*np.cos(angle) + n[:,7]*np.sin(angle)
    tangent = np.where(folded>4,-1,1)*(-n[:,6]*np.sin(angle)+n[:,7]*np.cos(angle))
    x, radius = xyz[ids[0,station],:2].T
    ellipse = ((x+1.7)/.8)**2 + ((radius-.55)/.45)**2
    pulse = 1e-3*np.maximum(1-ellipse,0.)**2*np.cos(np.linspace(0,2*np.pi,49)[ray])
    axial = base[mother[0,station],5]
    impulse = (n[:,5]+((axial+pulse)-axial)) - n[:,5]
    relaxed = work / "relax6"
    relaxed.mkdir()
    os.link(coarse/"mesh.su2",relaxed/"mesh.su2")
    state = np.column_stack((nodes,xyz,n[:,4:6],radial*np.cos(phi)-tangent*np.sin(phi),
                             radial*np.sin(phi)+tangent*np.cos(phi),n[:,8:10]))
    write_seed(relaxed/"initial_flow.csv",state)
    identical(relaxed/"initial_flow.csv","939e864fdad2caac2a38c61a72c66eb2a3a7c9725428cdf44426a583333ea559")
    np.save(relaxed/"impulse.npy",impulse)
    identical(relaxed/"impulse.npy","d4094998e05206ca8f7473a468e83874098dcce64088ce1606a547c9308e75fb")
    shutil.copyfile(identical(CASE/"relax.cfg","f758faa53e7e7549aee627f9a958d14e543c4da85e833915a821ad0573980d51"),relaxed/"case.cfg")
    native(relaxed,2,1200,768)
    state = columns(identical(relaxed/"restart.csv",RELAX_FRAME),FIELDS)
    state[:,5] += impulse
    write_seed(work/"seed_00000.csv",state)
    identical(work/"seed_00000.csv",COARSE_SEED)
    write(work/"seed_recovery.json",dict(source="internal-cold", baseSHA256=BASE_FRAME,
          relaxationSHA256=RELAX_FRAME,seedSHA256=COARSE_SEED))
