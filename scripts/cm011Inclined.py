"""生成并运行 +10° 的 CM-011 低速等温试件；整个迎流与背流侧由同一半域求解。"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

import numpy as np
import psutil

from cm011Tunnel import D, ROOT, U, meridian, surface, write_json
from cm011TunnelSector import mesh, write_volume
from runEmbeddedSU2 import sha


ALPHA = 10.
RAYS = np.deg2rad(np.linspace(-90.,90.,13))
FLOW = np.array([-np.cos(np.deg2rad(ALPHA)),0.,np.sin(np.deg2rad(ALPHA))])
TEMPLATE = ROOT/"docs/validation/v271/evidence/process/coarse03/case.cfg"
BODY_SHA = "07588689cfec2d88933a516520cdb905eff1976403b0990bb86ac0fc12946c19"
OFFICIAL_SHA = "01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"


def config(iterations):
    text = TEMPLATE.read_text("ascii")
    replacements = {"AXISYMMETRIC":"NO","MARKER_SYM":"(mirror)","INC_VELOCITY_INIT":
                    "("+", ".join(f"{v:.17g}" for v in FLOW)+")","CFL_NUMBER":"18","ITER":str(iterations),
                    "CONV_RESIDUAL_MINVAL":"-8.4","OUTPUT_WRT_FREQ":f"(500, 500, {iterations})"}
    for key,value in replacements.items():
        text,count = re.subn("^"+key+r"= .*?$",key+"= "+value,text,flags=re.M)
        assert count==1,key
    assert text.count("RMS_VELOCITY-Y, RMS_ENTHALPY")==2
    return text.replace("RMS_VELOCITY-Y, RMS_ENTHALPY","RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY")


def half_domain_rays(sectors):
    """两种半周离散共用同一对称面和原有的 13 个方位。"""
    if sectors not in (12, 24):
        raise ValueError("CM-011 半周网格只支持 12 或 24 个角向单元")
    return np.deg2rad(np.linspace(-90., 90., sectors+1))


def prepare(folder, imported_body=None, registration=None, sectors=12):
    rays = half_domain_rays(sectors)
    if imported_body is None and sectors != 12:
        raise ValueError("加密对照需要导入试件的实际壁面与配准记录")
    folder.mkdir(parents=True,exist_ok=False)
    wall,bands = meridian(2)
    if imported_body is None:
        body = surface(wall,bands,2)
        body.update(meridian=wall.tolist(),diameterM=D,grid="fine",layers=42,extentD=12.)
    else:
        body = imported_body
        wall = np.asarray(body["meridian"])
    write_json(folder/"body.json",body)
    if imported_body is None:
        assert sha(folder/"body.json")==BODY_SHA
    coarse = wall[::2] if imported_body is not None else meridian(1)[0]
    _,xyz,_,_,ids,cells,markers,stats = mesh(rays,level=1,layers=24,wall=coarse)
    # 两端位于同一 XZ 平面，有同一外法向。共享轴线只受一个法向约束。
    markers["mirror"] = markers.pop("sym0")+markers.pop("sym1")
    for face in markers["mirror"]:
        points = xyz[list(face)]
        assert np.max(abs(points[:,1]))<1e-12
        assert np.cross(points[1]-points[0],points[-1]-points[0])[1]<0
    write_volume(folder,xyz,cells,markers)
    np.testing.assert_allclose(wall[::2],coarse,atol=1e-13,rtol=0)
    stats.update(boundaryFaces={k:len(v) for k,v in markers.items()},alphaDeg=ALPHA,betaDeg=0,
                 velocityMps=(U*FLOW).tolist(),grid="inclined",reflection="Y=0",copies=2,
                 originalFaces=len(body["F"]),originalBodySHA256=sha(folder/"body.json"),
                 rayRadians=rays.tolist(),meridianNodes=len(coarse),radialLayers=24)
    if registration is not None:
        stats["importedSurface"]=registration
    write_json(folder/"geometry.json",stats)
    np.savez_compressed(folder/"nodes.npz",xyz=xyz,ids=ids)
    fields = np.tile([0.,*FLOW,1.-298.15],(len(xyz),1))
    wall_ids = np.unique(ids[:,:len(coarse)])
    fields[wall_ids,1:4] = 0
    fields[wall_ids,4] = 2.-298.15
    np.savetxt(folder/"initial_flow.csv",np.column_stack((np.arange(len(xyz)),xyz,fields)),
               fmt=["%d"]+8*["%.17g"],delimiter=",",comments="",
               header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy"')
    (folder/"case.cfg").write_text(config(12000),encoding="ascii")
    locks = [Path(__file__),Path(__file__).with_name("cm011TunnelSector.py"),Path(__file__).with_name("cm011Tunnel.py"),
             TEMPLATE,folder/"geometry.json",folder/"body.json",folder/"nodes.npz"]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               maximumSeconds=720,threads=1,jobToken=uuid.uuid4().hex,
               sourceFiles=[p.resolve().relative_to(ROOT).as_posix() for p in locks]))


def solve(folder):
    """观测本次求解进程树和输出；超预算时交由原生监督入口终止并写回执。"""
    request = json.loads((folder/"engine_request.json").read_bytes())
    if (folder/"engine_process.json").exists():
        raise ValueError("此目录已有运行记录")
    seen,peak,peak_disk,reason = {},0,0,None
    start = time.monotonic()
    proc = subprocess.Popen([sys.executable,str(ROOT/"scripts/runEmbeddedSU2.py"),"--request",str(folder/"engine_request.json")],
                            creationflags=subprocess.CREATE_NO_WINDOW)
    root = psutil.Process(proc.pid)
    while proc.poll() is None:
        try: members = [root]+root.children(recursive=True)
        except psutil.NoSuchProcess: members = []
        rss = 0
        for p in members:
            try:
                rss += p.memory_info().rss
                seen[p.pid] = p
            except psutil.NoSuchProcess: pass
        peak=max(peak,rss)
        peak_disk=max(peak_disk,sum(p.stat().st_size for p in folder.iterdir() if p.is_file()))
        if not reason and (peak>2*1024**3 or peak_disk>160*1024**2):
            reason="memory" if peak>2*1024**3 else "disk"
            write_json(folder/"engine_cancel.json",dict(jobToken=request["jobToken"]))
        time.sleep(.25)
    _,alive = psutil.wait_procs(list(seen.values()),timeout=3)
    write_json(folder/"resources.json",dict(elapsedSeconds=time.monotonic()-start,exitCode=proc.returncode,
               peakRSSBytes=peak,peakDiskBytes=peak_disk,stopReason=reason,livePids=[p.pid for p in alive]))
    if proc.returncode or reason or alive: raise RuntimeError("半域求解停止，见原生进程回执")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("prepare","solve"))
    parser.add_argument("directory",type=Path)
    args=parser.parse_args()
    folder=args.directory.resolve()
    folder.relative_to(ROOT/"outputs")
    if args.action=="prepare": prepare(folder)
    else: solve(folder)
