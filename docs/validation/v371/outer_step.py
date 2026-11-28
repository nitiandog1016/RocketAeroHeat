"""细分固定低速算例的两个外侧网格间隔，并在共同初场上推进一个物理步。"""
import argparse
from pathlib import Path
import shutil
import sys
import numpy as np
from scipy.spatial import cKDTree

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
BASE=ROOT/"outputs/v371_validation"
sys.path[:0]=[str(DOC.parent/"v369"),str(DOC.parent/"v367")]
from short_step import EXE,EXE_SHA,GRID_HASH,SEED_HASH,OLD,ORIGIN,RES,observe,numbers,read,save,sha,terminal,residual_window,vtk_tokens
from refine_tail import integral,write_volume_mesh,CONS

COARSE=ROOT/"outputs/v369_validation/pair01/half1"
CFG_HASH="ae6b1dfab53e1676eeb72c899e219c7e859717fe1eb35fb027250633604a56af"
COARSE_FRAME="7c783fbdee774c07729c86f6cc39daa778f84955a3111ef5de0fd10bea88f970"
LAYERS=np.sort(np.r_[np.arange(89),86.5,87.5])
SHAPE=(len(LAYERS),109)


def connectivity(rows):
    ids=np.arange(rows*109).reshape(rows,109)
    cells=np.stack((ids[1:,:-1],ids[1:,1:],ids[:-1,1:],ids[:-1,:-1]),axis=-1).reshape(-1,4)
    markers=dict(wall=np.column_stack((ids[0,1:],ids[0,:-1])),far=np.column_stack((ids[-1,:-1],ids[-1,1:])),
                 axis=np.vstack((np.column_stack((ids[:-1,0],ids[1:,0])),np.column_stack((ids[1:,-1],ids[:-1,-1])))))
    return cells,markers


def mesh_points(path,rows):
    """逐行核对单元方向、网格点和三个边界的线段与方向，取得原始求解坐标。"""
    cells,markers=connectivity(rows)
    with path.open(encoding="ascii") as f:
        assert f.readline().strip()=="NDIME= 2" and f.readline().strip()==f"NELEM= {len(cells)}"
        raw=np.array([f.readline().split() for _ in cells],int)
        assert (raw[:,0]==9).all()
        np.testing.assert_array_equal(raw[:,1:],cells)
        assert f.readline().strip()==f"NPOIN= {rows*109}"
        points=np.array([f.readline().split() for _ in range(rows*109)],float)
        assert f.readline().strip()=="NMARK= 3"
        for name,segments in markers.items():
            assert f.readline().strip()==f"MARKER_TAG= {name}" and f.readline().strip()==f"MARKER_ELEMS= {len(segments)}"
            raw=np.array([f.readline().split() for _ in segments],int)
            assert (raw[:,0]==3).all()
            np.testing.assert_array_equal(raw[:,1:],segments)
        assert not f.read().strip()
    return points.reshape(rows,109,2)


def insert_rows(values):
    lower=LAYERS.astype(int)
    upper=np.ceil(LAYERS).astype(int)
    w=(LAYERS-lower)[:,None,None]
    return (1-w)*values[lower]+w*values[upper]


def prepare(out):
    terminal(OLD);terminal(COARSE)
    assert sha(OLD/"mesh.su2")==GRID_HASH and sha(ORIGIN)==SEED_HASH and sha(COARSE/"restart_00001.csv")==COARSE_FRAME
    xy=mesh_points(OLD/"mesh.su2",89)
    coarse=numbers(ORIGIN)
    np.testing.assert_array_equal(coarse["PointID"],np.arange(9701))
    np.testing.assert_allclose(np.stack([coarse["x"],coarse["y"]],axis=-1),xy.reshape(-1,2),atol=2e-14,rtol=0)
    u=np.column_stack([coarse[k] for k in CONS]).reshape(89,109,4)
    points,new_u=insert_rows(xy),insert_rows(u)
    keep=LAYERS==np.floor(LAYERS)
    np.testing.assert_array_equal(points[keep],xy)
    np.testing.assert_array_equal(new_u[keep],u)
    totals,scales=integral(xy,u)
    fine_totals,_=integral(points,new_u)
    error=abs(fine_totals-totals)/scales
    assert max(error)<1e-12
    cells,markers=connectivity(SHAPE[0])
    corners=points.reshape(-1,2)[cells]
    sides=np.roll(corners,-1,axis=1)-corners
    following=np.roll(sides,-1,axis=1)
    assert np.min(sides[:,:,0]*following[:,:,1]-sides[:,:,1]*following[:,:,0])>0
    data=new_u.reshape(-1,4)
    assert min(data[:,0])>0 and min(data[:,3]-.5*(data[:,1]**2+data[:,2]**2)/data[:,0])>0
    out.mkdir(parents=True,exist_ok=False)
    mesh,seed=out/"mesh.su2",out/"initial.csv"
    write_volume_mesh(mesh,points.reshape(-1,2),cells,markers)
    np.testing.assert_array_equal(mesh_points(mesh,SHAPE[0]),points)
    np.savetxt(seed,np.column_stack((np.arange(len(data)),points.reshape(-1,2),data)),delimiter=",",
               comments="",fmt=["%d"]+["%.17g"]*6,header=",".join(f'"{k}"' for k in ("PointID","x","y",*CONS)))
    written=numbers(seed)
    np.testing.assert_array_equal(np.column_stack([written[k] for k in CONS]),data)
    save(out/"prepared.json",dict(originalMeshSHA256=GRID_HASH,originalSeedSHA256=SEED_HASH,
         fineMeshSHA256=sha(mesh),fineSeedSHA256=sha(seed),originalNodes=len(u.reshape(-1,4)),fineNodes=len(data),
         cells=len(cells),originalLayerCoordinates=LAYERS.tolist(),integrals=["volume","density","axialMomentum","radialMomentum","energy"],
         absoluteScales=scales.tolist(),oldIntegrals=totals.tolist(),newIntegrals=fine_totals.tolist(),relativeChanges=error.tolist()))


def run(prepared,out):
    manifest=read(prepared/"prepared.json")
    assert sha(EXE)==EXE_SHA and sha(COARSE/"case.cfg")==CFG_HASH
    mesh,seed=prepared/"mesh.su2",prepared/"initial.csv"
    assert sha(mesh)==manifest["fineMeshSHA256"] and sha(seed)==manifest["fineSeedSHA256"]
    out.mkdir(parents=True,exist_ok=False)
    for src,name in ((COARSE/"case.cfg","case.cfg"),(mesh,"mesh.su2"),(seed,"seed_00000.csv")):
        shutil.copyfile(src,out/name)
        assert sha(src)==sha(out/name)
    paths=[EXE,mesh,seed,COARSE/"case.cfg",OLD/"condition.json",prepared/"prepared.json",DOC/"PLAN.md",Path(__file__).resolve(),
           *(out/k for k in ("case.cfg","mesh.su2","seed_00000.csv"))]
    observe([str(EXE),"case.cfg","--threads","2"],out,{str(p):sha(p) for p in paths},150,64*1024**2)
    proc=terminal(out)
    h=numbers(out/"history_00001.csv");frame=out/"restart_00001.csv";wall=out/"surface_00001.vtk"
    np.testing.assert_array_equal(h["Inner_Iter"],np.arange(len(h["Inner_Iter"])))
    assert len(h["Inner_Iter"])<=2500 and (h["Time_Iter"]==1).all() and (h["Time_Step"]==6.25e-6).all() and (h["Cur_Time"]==6.25e-6).all()
    assert {k for k in h if k.startswith("rms[")}==set(RES)
    rec=residual_window(np.column_stack([h[k] for k in RES]))
    a,b=numbers(frame),numbers(seed)
    np.testing.assert_array_equal(a["PointID"],b["PointID"])
    for key in ("x","y"):
        np.testing.assert_allclose(a[key],b[key],atol=2e-14,rtol=0)
    rho=a["Density"]
    assert len(rho)==manifest["fineNodes"] and np.isfinite(np.column_stack(list(a.values()))).all() and min(rho)>0
    p=.4*(a["Energy"]-.5*(a["Momentum_x"]**2+a["Momentum_y"]**2)/rho)
    assert min(p)>0 and np.any(a["Energy"]!=b["Energy"])
    np.testing.assert_allclose(p,a["Pressure"],rtol=1e-10)
    np.testing.assert_allclose(p/rho,a["Temperature"],rtol=1e-10)
    assert not np.any(a["Momentum_x"][:109]) and not np.any(a["Momentum_y"][:109])
    xyz,fields=vtk_tokens(wall)
    d,ix=cKDTree(xyz[:,:2]).query(np.column_stack((a["x"][:109],a["y"][:109])))
    assert len(xyz)==len(np.unique(ix))==109 and max(d)<1e-6
    for native,vol in (("Pressure","Pressure"),("Heat_Flux","Heat_Flux"),("Temperature","Temperature"),
                       ("Skin_Friction_Coefficient","Skin_Friction_Coefficient_x")):
        np.testing.assert_allclose(fields[native][ix,0],a[vol][:109],rtol=6e-6,atol=1e-9)
    rec.update(elapsedUs=6.25,meshSHA256=sha(mesh),seedSHA256=sha(seed),configSHA256=CFG_HASH,frameSHA256=sha(frame),
               historySHA256=sha(out/"history_00001.csv"),wallSHA256=sha(wall),seconds=proc["elapsedSeconds"],peakRSSBytes=proc["peakFamilyRSSBytes"])
    save(out/"step.json",rec)
    print(rec["passed"],rec["iterations"],rec["last20Max"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared",type=Path);p.add_argument("--run",type=Path)
    a=p.parse_args();prepared=a.prepared.resolve();prepared.relative_to(BASE)
    if a.run:
        dest=a.run.resolve();dest.relative_to(BASE);run(prepared,dest)
    else:
        prepare(prepared)
