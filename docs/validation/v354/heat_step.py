"""同一低速时间层的三维显热、流体储热、原面壁压与入壁热流。"""
import argparse
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator,bicgstab,splu
from scipy.spatial import cKDTree

from flux_transfer import DOC,ROOT,build_model,project,read
from runEmbeddedSU2 import sha
from simpleLauncherAngularHeat import ANGLE,SECTORS,wall_results
from simpleLauncherHeatMesh import sampler
from simpleLauncherTemperature import transfer


def advance(m,operators,old,dt):
    matrices,right,left,boundary,projection=operators
    n=len(m["xy"]);axis=m["axis"];wet=m["wet"];far=m["far"]
    off=np.ones(n,bool);off[axis]=False
    cut=SECTORS*off.sum()
    def pack(a):
        return np.r_[a[:,off].ravel(),a[:,axis].mean(axis=0)]
    def unpack(v):
        a=np.empty((SECTORS,n));a[:,off]=v[:cut].reshape(SECTORS,-1);a[:,axis]=v[cut:]
        return a
    fixed=np.zeros((SECTORS,n),bool);fixed[:,far]=boundary[:,far]<0
    fixed[:,wet]=True
    # 共轴外界节点也只允许一个温度，按其完整控制体的净外流符号处理。
    shared_far=np.intersect1d(axis,far)
    fixed[:,shared_far]=boundary[:,shared_far].sum(axis=0)<0
    free=pack(~fixed).astype(bool);common_free=~fixed.all(axis=0)
    prescribed=np.zeros_like(old);prescribed[:,wet]=1.
    storage=m["mass"]/dt
    def spatial(theta):
        result=np.stack([a@t for a,t in zip(matrices,theta)])
        result+=right*(np.roll(theta,-1,axis=0)-theta)
        result+=left*(np.roll(theta,1,axis=0)-theta)
        return result
    def lift(a):
        t=np.zeros(len(free));t[free]=a
        return unpack(t)
    def matvec(x):
        t=lift(x)
        return pack(spatial(t)+storage*t)[free]
    size=int(free.sum());op=LinearOperator((size,size),matvec=matvec,dtype=float)
    rhs=pack(storage*(old-prescribed)-spatial(prescribed))[free]
    # 轴对称分块仅作线性预条件。完整残差含随方位变化的速度和周向平流项。
    factors=[]
    for k in range(SECTORS//2+1):
        mask=common_free.copy()
        if k:mask[axis]=False
        a=m["matrix"]+diags(storage+4*np.sin(k*ANGLE/2)**2*m["conductance"])
        factors.append((mask,splu(a[mask][:,mask].tocsc())))
    def inverse(v):
        f=np.fft.rfft(lift(v),axis=0);sol=np.zeros_like(f)
        for k,(mask,factor) in enumerate(factors):
            b=f[k,mask];part=factor.solve(np.column_stack((b.real,b.imag)))
            sol[k,mask]=part[:,0]+1j*part[:,1]
        return pack(np.fft.irfft(sol,n=SECTORS,axis=0))[free]
    prec=LinearOperator((size,size),matvec=inverse,dtype=float)
    iterations=[]
    x,status=bicgstab(op,rhs,x0=pack(old)[free],M=prec,tol=1e-14,atol=0.,maxiter=100,
                      callback=lambda _:iterations.append(1))
    theta=lift(x)+prescribed
    transport=spatial(theta);capacity=(theta-old)*storage;row=transport+capacity
    residual=np.abs(pack(row)[free]).max()/np.abs(rhs).max()
    # 从同一时间层分别累加壁行、开放远界的对流/入流边界反力以及全流体储热。
    inlet=fixed.copy();inlet[:,wet]=False
    wall=row[:,wet]
    remote=np.sum(boundary*theta)-row[inlet].sum()
    stored=capacity.sum()
    error=abs(wall.sum()-remote-stored)/abs(wall.sum())
    flux,integral=wall_results(m,wall)
    power=integral["outwardW"]/wall.sum()
    measure=dict(**projection,iterations=len(iterations),linearStatus=int(status),originalResidual=float(residual),
        minTheta=float(theta.min()),maxTheta=float(theta.max()),freeDOF=size,storageW=float(stored*power),
        farW=float(remote*power),relativeEnergyImbalance=float(error),**integral)
    good=bool(status==0 and np.isfinite(theta).all() and residual<1e-8 and error<.005
              and theta.min()>=-1e-12 and theta.max()<=1+1e-12)
    return theta,flux,measure,good


def surface_pressure(m,native):
    corners=m["corners"];diameter=m["diameter"]
    ids=m["flowIds"][:,:m["flowWallCount"]]
    tree=cKDTree(native[ids.ravel(),1:4]*diameter)
    distance,node=tree.query(corners.reshape(-1,3))
    assert distance.max()<1e-7
    pressure=native[ids.ravel()[node],4].reshape(-1,3).mean(axis=1)
    bottom=(corners[:,:,0]==0).all(axis=1)
    mid=corners[bottom].mean(axis=1)
    j=np.floor(np.mod(np.arctan2(mid[:,2],mid[:,1]),2*np.pi)/ANGLE).astype(int)
    base_end=np.flatnonzero(np.diff(m["wall"][:,0])>0)[0];base=ids[:,:base_end+1]
    r=m["wall"][:base_end+1,1]/m["wall"][base_end,1];a,b=r[:-1],r[1:]
    face_nodes=native[base,4]
    ray=np.sum((b-a)*((2*a+b)*face_nodes[:,:-1]+(a+2*b)*face_nodes[:,1:]),axis=1)/3
    pressure[bottom]=.5*(ray[j]+ray[(j+1)%SECTORS])
    taps=np.stack([np.interp([0.,.286,.572,.858],r,p) for p in face_nodes])
    return 101325+1.177*.05**2*pressure,taps*1.177*.05**2*1000


def near_wall(m,theta,prepared):
    """从真实热三角形插至第一层流体节点，再以同一角向扇区映回原面。"""
    ids=m["flowIds"];w=m["flowWallCount"]
    with np.load(prepared/"stations.npz") as data:
        xyz=data["xyz"][ids[0,w:2*w]]
    interp=sampler(m["xy"],m["tri"],np.column_stack((xyz[:,0],np.hypot(xyz[:,1],xyz[:,2]))))
    tri=m["corners"];mid=tri.mean(axis=1)
    sector=np.floor(np.mod(np.arctan2(mid[:,2],mid[:,1]),2*np.pi)/ANGLE).astype(int)
    mapping=transfer(tri,m["wall"],np.arange(w),m["diameter"])
    near=np.empty(len(tri))
    for j in range(SECTORS):
        mask=sector==j;near[mask]=300+.001*(mapping[mask]@interp(theta[j]))
    assert np.isfinite(near).all()
    return near


def run(flow_job,momentum,prepared,heat_root,name,initial_field=None,heat_folder=None):
    step=read(momentum/name/"step.json")
    frame=momentum/name/f"restart_{step['step']:05d}.csv"
    assert step["accepted"] and sha(frame)==step["frameSHA256"]
    assert step["meshSHA256"]==read(prepared/"prepared.json")["files"]["mesh.su2"]
    m=build_model(flow_job,prepared,heat_folder)
    initial=initial_field or ROOT/"outputs/v353_validation/coarse01/fields.npz"
    if step["step"]==1:
        proof=read(initial.parent/"results.json")
        assert proof["accepted"] and all(proof["gates"].values()) and sha(initial)==proof["fieldsSHA256"]
        assert proof["sourceFlowSHA256"]==m["flowSHA256"] and proof["sourceHeatMeshSHA256"]==m["heatMeshSHA256"]
        with np.load(initial) as f:old=f["theta"][0]
    else:
        parent=heat_root/f"{step['branch']}_{step['stepSize']}_{step['step']-1}"
        previous=read(parent/"heat.json");initial=parent/"thermal.npz"
        assert previous["accepted"] and sha(initial)==previous["fieldSHA256"]
        assert previous["sourceFrameSHA256"]==step["seedSHA256"]
        if heat_folder is not None:assert previous["heatMeshSHA256"]==m["heatMeshSHA256"]
        with np.load(initial) as f:old=f["theta"]
    dst=heat_root/name
    native=np.loadtxt(frame,delimiter=",",skiprows=1)
    operators=project(m,native)
    theta,q,summary,ok=advance(m,operators,old,step["dtDU"])
    p,taps=surface_pressure(m,native)
    tri=m["corners"]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    near=near_wall(m,theta,prepared)
    np.savez_compressed(dst/"thermal.npz",theta=theta)
    np.savetxt(dst/"wall.csv",np.column_stack((np.arange(1,len(q)+1),area,p,q,near)),delimiter=",",
        fmt=["%d"]+4*["%.17g"],header="faceID,areaM2,pressurePa,heatFluxWm2,nearWallK",comments="")
    summary.update(accepted=ok,branch=step["branch"],step=step["step"],stepSize=step["stepSize"],
        elapsedPhysicalS=step["elapsedPhysicalS"],sourceFrameSHA256=sha(frame),initialThetaSHA256=sha(initial),
        heatMeshSHA256=m["heatMeshSHA256"],
        fieldSHA256=sha(dst/"thermal.npz"),wallSHA256=sha(dst/"wall.csv"),baseTapGaugeMilliPa=taps.tolist())
    (dst/"heat.json").write_bytes((json.dumps(summary,indent=2,allow_nan=False)+"\n").encode())
    print(name,ok,"residual",summary["originalResidual"],"energy",summary["relativeEnergyImbalance"],
          "bounds",summary["minTheta"],summary["maxTheta"],flush=True)
    return ok


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("flow_job","momentum","prepared","heat_root"):p.add_argument(arg,type=Path)
    p.add_argument("name",choices=[f"{b}_{z}_{i}" for b in ("control","disturbed")
                                   for z,steps in (("full",1),("half",2)) for i in range(1,steps+1)])
    p.add_argument("--initial-field",type=Path)
    p.add_argument("--heat-folder",type=Path,help="同一粗速度终帧上单独建立的热网与稳态初场目录")
    a=p.parse_args()
    sys.exit(0 if run(a.flow_job.resolve(),a.momentum.resolve(),a.prepared.resolve(),a.heat_root.resolve(),a.name,a.initial_field,a.heat_folder) else 2)
