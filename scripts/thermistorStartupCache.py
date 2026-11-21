"""为同一个 6 m/s 热敏珠启动工况准备可重建的本机输入缓存。"""
from contextlib import contextmanager
import errno
import msvcrt
import os
from pathlib import Path
import shutil
import sys
import tempfile
import zipfile

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
import numpy as np

from cm011TunnelSector import mesh, write_volume
from thermistorBead import D, ROOT, wall, read, sha
from thermistorNearWall import FIELDS, geometry, read_grid, refine
from thermistorSeedRecovery import CASE, COARSE_SEED, identical, recover, write_seed
from thermistorStartupSurface import columns, write
from runSU2NativeDebugger import observe

CACHE = ROOT/"outputs/thermistor_inputs/B6-56"
COARSE = CACHE/"mesh48"
GRID = CACHE/"grid56"
BUILD = CACHE/"observer"
PREVIOUS = ROOT/"docs/validation/v293"
GRID_SPEC = read(PREVIOUS/"evidence/grid.json")
OBSERVER_SPEC = read(PREVIOUS/"evidence/readout_build.json")
COARSE_MESH = "673151f2681de2a150f7de689bea3096e3b6fa6d9b91ddb80dc635e563fe951c"
COARSE_NODES = "2617470c05476d59a00ad33c509d6c9a5a2db7ce1d22b0ce86ae61d9ae0392b6"
BODY_SHA = "5e66bb2f70f0774b8ebb87a77b76061fc84b89b9f379d54e90879e81628bf62f"
GEOM_SHA = "1587ac97a3ed59ab75ce2b35266039c46c58ef837480f66fcb6cdeb2a2cfb7fe"


def supervised(stage, inputs):
    work=Path(tempfile.mkdtemp(prefix=stage+"-",dir=CACHE))
    script=Path(__file__).resolve()
    locks=[script,ROOT/"scripts/thermistorNearWall.py",ROOT/"scripts/thermistorSeedRecovery.py",
           ROOT/"scripts/cm011TunnelSector.py",PREVIOUS/"evidence/grid.json",*inputs]
    rec=observe([sys.executable,str(script),"--worker",stage,str(work)],work,
                {str(p):sha(p) for p in locks},180,200*1024**2)
    if rec["debuggerExitCode"] or rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]:
        raise RuntimeError("热敏珠输入生成未完成，现场 "+str(work))
    return work


def publish(target, work, files):
    """保持原生工作目录可追溯，只把已核对的文件发布到稳定缓存路径。"""
    temp=Path(tempfile.mkdtemp(prefix="publish-",dir=CACHE))
    for name in files:
        link_or_copy(work/name,temp/name)
    write(temp/"installed.json",dict(sourceDirectory=str(work)))
    temp.rename(target)


def link_or_copy(source, dest):
    try:
        os.link(source,dest)
    except OSError as exc:
        if exc.errno!=errno.EXDEV and getattr(exc,"winerror",None)!=17:
            raise
        # 跨卷时先复制并核验。即使磁盘写满，中断的种子也不会占用正式缓存文件名。
        with tempfile.TemporaryDirectory(dir=dest.parent, prefix="transfer-") as folder:
            pending=Path(folder)/dest.name
            shutil.copyfile(source,pending)
            identical(pending,sha(source))
            pending.rename(dest)


def available(folder, files):
    """旧安装是可选的迁移来源。保留损坏文件，从登记配方重建缺失输入。"""
    for name,digest in files:
        file=folder/name
        if not file.is_file():
            return False
        if sha(file)!=digest:
            print("旧输入指纹不匹配，将从计算配方恢复："+str(file),flush=True)
            return False
    return True


@contextmanager
def installation_lock(name="install"):
    CACHE.mkdir(parents=True,exist_ok=True)
    with (CACHE/(name+".lock")).open("a+b") as f:
        f.seek(0)
        try:
            msvcrt.locking(f.fileno(),msvcrt.LK_NBLCK,1)
        except OSError as exc:
            raise RuntimeError("另一任务正在准备热敏珠输入，请在准备结束后使用缓存。") from exc
        try:
            yield
        finally:
            f.seek(0)
            msvcrt.locking(f.fileno(),msvcrt.LK_UNLCK,1)


def create_coarse(out):
    """由名义尺寸建立全周网格和节点映射，从随工程发布的同一几何恢复原面编号。"""
    _,xyz,_,_,ids,cells,markers,_ = mesh(np.linspace(0,2*np.pi,49),2,48,wall(64),20.,D)
    write_volume(out,xyz,cells,markers)
    np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=ids)
    identical(out/"mesh.su2",COARSE_MESH)
    identical(out/"nodes.npz",COARSE_NODES)
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/startup_6ms.zip") as archive:
        (out/"body.json").write_bytes(archive.read("body.json"))
    identical(out/"body.json",BODY_SHA)


def coarse_grid():
    files=(("mesh.su2",COARSE_MESH),("nodes.npz",COARSE_NODES),("body.json",BODY_SHA))
    if not COARSE.exists():
        work=ROOT/"outputs/v289_validation/window01"
        if not available(work,files):
            work=supervised("mesh48",[ROOT/"cases/thermistor_b/startup_6ms.zip",ROOT/"scripts/thermistorBead.py"])
        publish(COARSE,work,[name for name,_ in files])
    for name,digest in files:
        identical(COARSE/name,digest)


def bound_source(name):
    files = {"geometry.json": (CASE/"geometry.json",GEOM_SHA), "body.json": (COARSE/"body.json",BODY_SHA),
             "nodes.npz": (COARSE/"nodes.npz",COARSE_NODES)}
    return identical(*files[name])


def fine_mesh(out):
    """在原网最内八个区间增加中点，重算对偶体积；网格步骤只生成几何。"""
    xyz,cells,markers = read_grid(COARSE/"mesh.su2")
    with np.load(COARSE/"nodes.npz") as z:
        assert np.array_equal(z["xyz"],xyz)
        ids=z["ids"].reshape(49,49,65)
    new_xyz,fine,parent,_,new_cells = refine(xyz,ids,cells)
    write_volume(out,new_xyz,new_cells,markers)
    split=np.bincount(parent,minlength=len(cells))==2
    dual,_,_,geom=geometry(new_xyz,new_cells,markers,split[parent])
    assert geom==GRID_SPEC["newGeometry"]
    _,indices=np.unique(ids[:,0].ravel(),return_index=True)
    first=ids[:,1].ravel()[indices]
    near=fine[:,1].ravel()[indices]
    np.savez_compressed(out/"nodes.npz",xyz=new_xyz,ids=fine.reshape(49,-1),originalFirst=first,
                        refinedFirst=near,oldGridNodes=len(xyz),dualVolume=dual)
    # 复用同字节网格已有的科学审核记录。初场身份在计算步骤独立核对。
    for name,key in (("mesh.su2","gridSHA256"),("nodes.npz","nodeMapSHA256")):
        identical(out/name,GRID_SPEC[key])
    shutil.copyfile(PREVIOUS/"evidence/grid.json",out/"grid.json")


def mesh_info(grid=GRID):
    assert grid.resolve()==GRID and read(grid/"grid.json")==GRID_SPEC
    for name,key in (("mesh.su2","gridSHA256"),("nodes.npz","nodeMapSHA256")):
        identical(grid/name,GRID_SPEC[key])
    return GRID_SPEC


def grid_info(grid=GRID):
    info=mesh_info(grid)
    identical(grid/"seed_00000.csv",info["seedSHA256"])
    return info


def ensure_mesh():
    with installation_lock():
        coarse_grid()
        if not GRID.exists():
            work=ROOT/"outputs/v293_validation/grid01"
            if not available(work,(("mesh.su2",GRID_SPEC["gridSHA256"]),("nodes.npz",GRID_SPEC["nodeMapSHA256"]),
                                   ("grid.json","f796d2d3e3aee6682e2196ac625ab439e1098516bc0dac0832260952faecbf0b"))):
                work=supervised("grid56",[COARSE/"mesh.su2",COARSE/"nodes.npz"])
            publish(GRID,work,("mesh.su2","nodes.npz","grid.json"))
        return mesh_info()


def prolong(seed, dest):
    """按缓存中的真实节点编号恢复新增中点的两个父节点，保持旧节点的六分量。"""
    with np.load(GRID/"nodes.npz") as z:
        ids=z["ids"].reshape(49,57,65)
        xyz=z["xyz"]
    old=columns(identical(seed,COARSE_SEED),FIELDS)
    middle=ids[:,1:16:2].ravel()-len(old)
    pairs=np.empty((len(xyz)-len(old),2),int)
    assert np.array_equal(np.unique(middle),np.arange(len(pairs)))
    pairs[middle,0]=ids[:,0:16:2].ravel()
    pairs[middle,1]=ids[:,2:18:2].ravel()
    assert np.array_equal(xyz[len(old):],xyz[pairs].mean(axis=1))
    values=np.vstack((old[:,4:],old[pairs,4:].mean(axis=1)))
    write_seed(dest,np.column_stack((np.arange(len(xyz)),xyz,values)))
    identical(dest,GRID_SPEC["seedSHA256"])


def ensure_grid():
    ensure_mesh()
    with installation_lock():
        target=GRID/"seed_00000.csv"
        if not target.exists():
            old_fine=ROOT/"outputs/v293_validation/grid01"/target.name
            if available(old_fine.parent,[(target.name,GRID_SPEC["seedSHA256"])]):
                link_or_copy(old_fine,target)
                return grid_info()
            seed=CACHE/"initial"/target.name
            legacy=ROOT/"outputs/v289_validation/window01"/target.name
            if not seed.exists():
                work=Path(tempfile.mkdtemp(prefix="initial-",dir=CACHE))
                if available(legacy.parent,[(target.name,COARSE_SEED)]):
                    shutil.copyfile(identical(legacy,COARSE_SEED),work/target.name)
                    write(work/"seed_recovery.json",dict(source="verified-local-migration",seedSHA256=COARSE_SEED))
                else:
                    recover(work,COARSE)
                publish(seed.parent,work,(target.name,"seed_recovery.json"))
            work=supervised("prolong",[seed,GRID/"nodes.npz"])
            os.link(work/target.name,target)
        return grid_info()


def install_seed(file):
    """把另一安装中由程序计算的同一 48 层初态装入缓存，不接受其他物理场。"""
    source=identical(file.resolve(),COARSE_SEED)
    with installation_lock():
        target=CACHE/"initial"
        if target.exists():
            identical(target/"seed_00000.csv",COARSE_SEED)
            return
        work=Path(tempfile.mkdtemp(prefix="publish-initial-",dir=CACHE))
        shutil.copyfile(source,work/"seed_00000.csv")
        identical(work/"seed_00000.csv",COARSE_SEED)
        write(work/"installed.json",dict(source="verified-internal-checkpoint",seedSHA256=COARSE_SEED))
        work.rename(target)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("热敏珠输入审核需要开启 Python 断言")
    import argparse
    cli=argparse.ArgumentParser(description=__doc__)
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--install-seed",type=Path,help="可选：复用另一安装中已生成的准确旧态")
    mode.add_argument("--worker",choices=("mesh48","grid56","prolong"))
    cli.add_argument("work",type=Path,nargs="?")
    args=cli.parse_args()
    if args.worker:
        work=args.work.resolve()
        assert work.parent==CACHE and work.name.startswith(args.worker+"-")
        if args.worker=="mesh48":create_coarse(work)
        elif args.worker=="grid56":fine_mesh(work)
        else:prolong(CACHE/"initial/seed_00000.csv",work/"seed_00000.csv")
    elif args.install_seed is None:
        ensure_grid()
    else:
        install_seed(args.install_seed)
