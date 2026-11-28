"""在既有低速全域中二分三个子午间隔，复用原物理步和面积映射。"""
import argparse
from pathlib import Path
import shutil
import numpy as np
from match_state import PARENT,BASE,DOC,ROOT,EXE,EXE_SHA,RES,numbers,read,save,sha,terminal,observe
from physical_step import STEPS, step_config
from energy_budget import patch,GAUSS
from cm011Tunnel import write_volume_mesh
from read_wall import mapping
from simpleLauncherTemperature import transfer

KNOTS=np.sort(np.r_[np.arange(109),11.5,12.5,13.5])
CONS=("Density","Momentum_x","Momentum_y","Energy")


def prolong(values):
    lo=KNOTS.astype(int);hi=np.ceil(KNOTS).astype(int)
    w=(KNOTS-lo)[None,:,None]
    return (1-w)*values[:,lo]+w*values[:,hi]


def integral(xy,values):
    sums=np.zeros(5);scales=np.zeros(5)
    for s in GAUSS:
        for t in GAUSS:
            x,xs,xt=patch(xy,s,t)
            weight=.5*np.pi*x[:,:,1]*(xs[:,:,1]*xt[:,:,0]-xs[:,:,0]*xt[:,:,1])
            assert (weight>0).all()
            u=patch(values,s,t)[0]
            sums[0]+=weight.sum();scales[0]+=weight.sum()
            sums[1:]+=(weight[:,:,None]*u).sum(axis=(0,1))
            scales[1:]+=(weight[:,:,None]*abs(u)).sum(axis=(0,1))
    return sums,scales


def prepare(dest):
    dest=dest.resolve();dest.relative_to(BASE)
    original=PARENT/"mesh.su2";data=PARENT/"seed_00000.csv"
    assert sha(original)=="f5656636b80e5a47a7d0cdb04b52f7a97e24fefb274f298ae1c31f3a42e4da68"
    assert sha(data)=="998600a5f54616298783e53f78370ecd51b203f64e639835f598d4d8a7871346"
    match=read(BASE/"matched01/sample.json");assert max(match["logDifferences"].values())<.001
    lines=original.read_text().splitlines();row=lines.index("NPOIN= 9701")+1
    xy=np.array([line.split()[:2] for line in lines[row:row+9701]],float).reshape(89,109,2)
    a=numbers(data)
    np.testing.assert_allclose(np.column_stack((a["x"],a["y"])),xy.reshape(-1,2),rtol=0,atol=2e-14)
    values=np.column_stack([a[k] for k in CONS]).reshape(89,109,4)
    pts,u=prolong(xy),prolong(values)
    np.testing.assert_array_equal(pts[:,KNOTS==np.floor(KNOTS)],xy)
    np.testing.assert_array_equal(u[:,KNOTS==np.floor(KNOTS)],values)
    coarse,scale=integral(xy,values);fine,_=integral(pts,u);err=abs(fine-coarse)/scale
    assert max(err)<1e-12
    assert u[:,:,0].min()>0 and np.min(u[:,:,3]-.5*(u[:,:,1]**2+u[:,:,2]**2)/u[:,:,0])>0
    ids=np.arange(pts.shape[0]*pts.shape[1]).reshape(pts.shape[:2])
    cells=np.stack((ids[1:,:-1],ids[1:,1:],ids[:-1,1:],ids[:-1,:-1]),axis=-1).reshape(-1,4)
    edge=lambda x,y: np.column_stack((x,y))
    markers=dict(wall=edge(ids[0,1:],ids[0,:-1]),far=edge(ids[-1,:-1],ids[-1,1:]),
        axis=np.vstack((edge(ids[:-1,0],ids[1:,0]),edge(ids[1:,-1],ids[:-1,-1]))))
    corners=pts.reshape(-1,2)[cells];sides=np.roll(corners,-1,axis=1)-corners
    assert np.all(np.cross(sides,np.roll(sides,-1,axis=1))>0)
    dest.mkdir(parents=True,exist_ok=False)
    write_volume_mesh(dest/"mesh.su2",pts.reshape(-1,2),cells,markers)
    np.savetxt(dest/"initial.csv",np.column_stack((ids.ravel(),pts.reshape(-1,2),u.reshape(-1,4))),
        delimiter=",",fmt=["%d"]+["%.17g"]*6,comments="",header=",".join(f'"{k}"' for k in ("PointID","x","y",*CONS)))
    save(dest/"prepared.json",dict(seedSHA256=sha(data),parentMeshSHA256=sha(original),meshSHA256=sha(dest/"mesh.su2"),
        initialSHA256=sha(dest/"initial.csv"),originalNodes=9701,nodes=ids.size,cells=len(cells),oldKnots=109,knots=KNOTS.tolist(),
        integralQuantities=["volume","density","momentum-x","momentum-r","energy"],relativeIntegralDifferences=err.tolist()))


def run(prepared,out,kind):
    prepared=prepared.resolve();out=out.resolve()
    for p in (prepared,out):p.relative_to(BASE)
    init=read(prepared/"prepared.json");assert sha(EXE)==EXE_SHA and sha(prepared/"mesh.su2")==init["meshSHA256"]
    dt,i=STEPS[kind];source=prepared/"initial.csv"
    assert sha(source)==init["initialSHA256"]
    if i==2:
        parent=out/"half1";r=read(parent/"step.json");terminal(parent)
        source=parent/"restart_00001.csv"
        assert r["passed"] and r["seedSHA256"]==init["initialSHA256"] and sha(source)==r["frameSHA256"]
    folder=out/kind;folder.mkdir(parents=True,exist_ok=False)
    seed=folder/f"seed_{i-1:05d}.csv";mesh=folder/"mesh.su2";cfg=folder/"case.cfg"
    shutil.copyfile(source,seed);shutil.copyfile(prepared/"mesh.su2",mesh)
    cfg.write_text("".join(f"{k}= {v}\n" for k,v in step_config("newton",kind).items()),encoding="ascii")
    paths=[EXE,source,seed,mesh,cfg,prepared/"prepared.json",DOC/"REFINEMENT.md",Path(__file__).resolve(),DOC.parent/"v366/physical_step.py"]
    observe([str(EXE),cfg.name,"--threads","2"],folder,{str(p):sha(p) for p in paths},150,64*1024**2)
    rc=terminal(folder);history=folder/f"history_{i:05d}.csv";frame=folder/f"restart_{i:05d}.csv"
    h=numbers(history);f=numbers(frame);seed_data=numbers(seed)
    np.testing.assert_array_equal(h["Inner_Iter"],np.arange(len(h["Inner_Iter"])))
    assert np.all(h["Time_Iter"]==i) and np.all(h["Time_Step"]==dt)
    assert len(f["PointID"])==init["nodes"] and np.isfinite(np.column_stack(list(f.values()))).all()
    for key in ("PointID","x","y"):
        np.testing.assert_allclose(f[key],seed_data[key],rtol=0,atol=2e-14)
    assert min(f["Density"])>0 and min(f["Pressure"])>0 and min(f["Temperature"])>0
    assert np.max(abs(f["Velocity_x"][:len(KNOTS)]))==np.max(abs(f["Velocity_y"][:len(KNOTS)]))==0
    last={k:float(max(h[k][-20:])) for k in RES}
    native=f["Residual_Energy"].reshape(89,-1)
    # 转置限制把新增节点的右端分给父段两端，保留带符号总和。
    restricted=np.zeros((89,109))
    for j,knot in enumerate(KNOTS):
        left=int(knot);w=knot-left
        restricted[:,left]+=(1-w)*native[:,j]
        if w:restricted[:,left+1]+=w*native[:,j]
    np.testing.assert_allclose(restricted.sum(),native.sum(),rtol=1e-12,atol=1e-14)
    model,tri,area,_,wall=mapping()
    matrix=transfer(tri,wall,KNOTS,1.)
    wallq=matrix@f["Heat_Flux"][:len(KNOTS)]
    assert len(wallq)==9216 and np.isfinite(wallq).all()
    energy_square=restricted**2
    rec=dict(kind=kind,dtSeconds=dt,elapsedFromSeedSeconds=i*dt,nativeTimeSeconds=float(h["Cur_Time"][-1]),inner=len(h["Inner_Iter"]),
        last20Max=last,passed=bool(len(h["Inner_Iter"])>=20 and max(last.values())< -8),frameSHA256=sha(frame),
        seedSHA256=sha(seed),meshSHA256=sha(mesh),historySHA256=sha(history),seconds=rc["elapsedSeconds"],peakRSSBytes=rc["peakFamilyRSSBytes"],
        restrictedEnergyLogRms=float(.5*np.log10(energy_square.mean())),restrictedTailShare=float(energy_square[:,12:14].sum()/energy_square.sum()),
        originalFaces=len(area),candidateCylinderIntoWallW=float(area[2016:6912]@wallq[2016:6912]))
    save(folder/"step.json",rec);print(rec,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("prepared",type=Path)
    p.add_argument("mode",choices=("prepare",*STEPS));p.add_argument("--out",type=Path)
    a=p.parse_args()
    if a.mode=="prepare":prepare(a.prepared)
    else:run(a.prepared,a.out,a.mode)
