"""在热敏珠的同工况三维扇区中推进全部动量和能量分量。"""
import argparse
from pathlib import Path
import shutil
import uuid

import numpy as np
from cm011TunnelSector import mesh, write_volume
from thermistorBead import D, GRIDS, ROOT, body, read, sha, solve, wall, write_json

RAYS = np.deg2rad([0.,15.,30.])
COPIES = 12


def sector_rays(angular_cells):
    """同一 30° 扇区的两种已登记离散，端点始终为对称面。"""
    if angular_cells not in (2, 4):
        raise ValueError("热敏珠扇区支持 2 或 4 个周向单元")
    return np.deg2rad(np.linspace(0., 30., angular_cells+1))


def prepare_workbench(out, angular_cells):
    """在边界确认前生成真实的体网及壁面预览，输入仅含固定珠体几何和分辨率。"""
    m,n,extent,level = GRIDS["outer"]
    rays = sector_rays(angular_cells)
    _,xyz,_,_,ids,cells,markers,geom = mesh(rays,level,n,wall(m),extent,D)
    out.mkdir(parents=True,exist_ok=False)
    write_volume(out,xyz,cells,markers)
    np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=ids)
    fixed = read(ROOT/"cases/thermistor_b/grids.json")["grids"]["outer" if angular_cells==2 else "outerAngular4"]
    for name in ("mesh.su2","nodes.npz"):
        assert sha(out/name)==fixed[name], name
    # 只把原生壁面边界三角化用于预览；体网中的四边形面保留原拓扑。
    triangles=[]
    for face in markers["wall"]:
        triangles.append([face[0],face[1],face[2]])
        if len(face)==4:
            triangles.append([face[0],face[2],face[3]])
    nodes,local=np.unique(triangles,return_inverse=True)
    geom.update(angularCells=angular_cells,meshSHA256=fixed["mesh.su2"],nodesSHA256=fixed["nodes.npz"],
                preview=dict(V=(D*xyz[nodes]).tolist(),F=(local.reshape(-1,3)+1).tolist()))
    write_json(out/"prepared.json",geom)


def use_workbench_grid(directory, expected, out, angular_cells, rays):
    """复用用户确认的确切体网字节，在本次独立计算目录核对输入并恢复初猜映射。"""
    assert sha(directory/"prepared.json")==expected
    geom=read(directory/"prepared.json")
    fixed=read(ROOT/"cases/thermistor_b/grids.json")["grids"]["outer" if angular_cells==2 else "outerAngular4"]
    assert geom.pop("angularCells")==angular_cells
    assert geom.pop("meshSHA256")==fixed["mesh.su2"] and geom.pop("nodesSHA256")==fixed["nodes.npz"]
    geom.pop("preview")
    for name in ("mesh.su2","nodes.npz"):
        shutil.copyfile(directory/name,out/name)
        assert sha(out/name)==fixed[name], name
    with np.load(out/"nodes.npz") as data:
        ids,xyz=data["ids"],data["xyz"]
    source=np.empty(len(xyz),int);phi=np.empty(len(xyz))
    for k in reversed(range(len(rays))):
        source[ids[k]]=np.arange(ids.shape[1])
        phi[ids[k]]=rays[k]
    return xyz[ids[0],:2],xyz,source,phi,geom


def prepare(parent, out, angular_cells=2, prepared=None, expected=None):
    spec = read(parent/"case.json")
    rec = read(parent/"engine_process.json")
    assert rec["status"]=="terminal" and rec["exitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["timedOut"] or rec["cancelled"] or rec.get("supervisorError"))
    for name in ("case.json","body.json","mesh.su2","case.cfg"):
        assert sha(parent/name)==rec["inputSHA256"][str(parent/name)]
    assert sha(parent/"restart.csv")==rec["outputSHA256"]["restart.csv"] and read(parent/"body.json")==body()
    m,n,extent,level = GRIDS[spec["grid"]]
    rays = sector_rays(angular_cells)
    if angular_cells == 4 and spec["grid"] != "outer":
        raise ValueError("7.5° 对照使用 20D 外域")
    out.mkdir(parents=True,exist_ok=False)
    if prepared is None:
        xy,xyz,source,phi,ids,cells,markers,geom = mesh(rays,level,n,wall(m),extent,D)
        write_volume(out,xyz,cells,markers)
        np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=ids)
    else:
        assert spec["grid"]=="outer"
        xy,xyz,source,phi,geom=use_workbench_grid(prepared,expected,out,angular_cells,rays)
    native = np.loadtxt(parent/"restart.csv",delimiter=",",skiprows=1)
    assert np.array_equal(native[:,0],np.arange(len(xy)))
    np.testing.assert_allclose(native[:,1:3],xy,rtol=0,atol=1e-13)
    p,u,v,h = native[source,3:7].T
    z = h+298.15
    assert np.isfinite(native).all() and z.min()>=1-1e-10 and z.max()<=2+1e-10
    state = np.column_stack((np.arange(len(xyz)),xyz,p,u,v*np.cos(phi),v*np.sin(phi),h))
    np.savetxt(out/"initial_flow.csv",state,fmt=["%d"]+["%.17g"]*8,delimiter=",",comments="",
               header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy"')
    cfg=(parent/"case.cfg").read_text("ascii")
    for a,b in (("AXISYMMETRIC= YES","AXISYMMETRIC= NO"),("MARKER_SYM= (axis)","MARKER_SYM= (sym0, sym1)"),
                ("RMS_VELOCITY-Y, RMS_ENTHALPY","RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY")):
        assert a in cfg
        cfg=cfg.replace(a,b)
    cfg+="RESTART_SOL= YES\nREAD_BINARY_RESTART= NO\nSOLUTION_FILENAME= initial_flow\nRESTART_FILENAME= restart\n"
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    write_json(out/"body.json",body())
    write_json(out/"case.json",spec)
    geom.update(parentFrameSHA256=sha(parent/"restart.csv"),parentBodySHA256=sha(parent/"body.json"),
                raysRadians=rays.tolist(),rotationCopies=COPIES)
    write_json(out/"geometry.json",geom)
    inputs=[Path(__file__).resolve(),ROOT/"scripts/thermistorBead.py",ROOT/"scripts/cm011Tunnel.py",
            ROOT/"scripts/cm011TunnelSector.py",parent/"restart.csv",out/"body.json",out/"case.json",
            out/"initial_flow.csv",out/"nodes.npz",out/"geometry.json"]
    if prepared is not None:
        inputs.extend(prepared/name for name in ("mesh.su2","nodes.npz","prepared.json"))
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
               threads=1,maximumSeconds=360,jobToken=str(uuid.uuid4()),
               sourceFiles=[p.relative_to(ROOT).as_posix() for p in inputs]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent",type=Path)
    parser.add_argument("new_folder",type=Path)
    parser.add_argument("--angular-cells",type=int,choices=(2,4),default=2)
    args=parser.parse_args()
    out=args.new_folder.resolve()
    out.relative_to(ROOT/"outputs")
    prepare(args.parent.resolve(),out,args.angular_cells)
    solve(out)
