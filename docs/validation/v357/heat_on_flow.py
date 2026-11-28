"""在原 0.05 m/s 热网和原初温上，使用加密动量网原生速度及壁压。"""
import argparse
import sys
from pathlib import Path
from refine_flow import BASE,DOC,OLD,ROOT,SEED,UP,read,save,sha
import numpy as np
from flux_transfer import SECTORS,boundary_sampler,build_model
from heat_step import run_on_mesh,surface_pressure
from simpleLauncherHeatMesh import sampler
from simpleLauncherTemperature import transfer
from wall_sampling import original_face_map
from runSU2NativeDebugger import observe


def attach(prepared):
    g=read(prepared/"refinement.json")
    assert sha(prepared/"meridian.npz")==g["meridianSHA256"]
    for name in ("control.csv","stations.npz"):
        assert sha(prepared/name)==g["newFiles"][name]
    assert sha(SEED/"steady/temperature.npz")==g["heatMesh"] and sha(SEED/"initial/fields.npz")==g["savedInitial"]
    m=build_model(UP,OLD,SEED/"steady")
    coarse=m["seedVelocity"][m["flowIds"][0]][:,:2]
    with np.load(prepared/"meridian.npz") as saved:
        coords=saved["xy"];knots=saved["station"]
    width=len(knots);xy=coords.reshape(-1,2)
    grid=np.arange(len(xy)).reshape(coords.shape[:2])
    quads=np.stack((grid[:-1,:-1],grid[:-1,1:],grid[1:,1:],grid[1:,:-1]),axis=-1).reshape(-1,4)
    tris=np.vstack((quads[:,[0,1,2]],quads[:,[0,2,3]]))
    axis=xy[:,1]==0;ids=np.empty((SECTORS,len(xy)),int)
    ids[:,axis]=np.arange(axis.sum());ids[:,~axis]=axis.sum()+np.arange(SECTORS*(~axis).sum()).reshape(SECTORS,-1)
    with np.load(prepared/"stations.npz") as data:
        np.testing.assert_array_equal(ids[:,:width],data["wall"])
        station_xyz=data["xyz"][data["wall"]]*m["diameter"]
        for j,line in enumerate(ids):
            phi=j*2*np.pi/SECTORS
            actual=data["xyz"][line]
            np.testing.assert_allclose(actual,np.column_stack((xy[:,0],xy[:,1]*np.cos(phi),xy[:,1]*np.sin(phi))),atol=1e-13,rtol=0)
    velocity=np.loadtxt(prepared/"control.csv",delimiter=",",skiprows=1,usecols=(5,6,7))
    control=velocity[ids[0],:2]
    # 热算子的基态是旧稳态守恒通量。显式加上两种网格对同一控制种子插值的差，
    # 使输入等于新网原生速度与旧热算子速度在相同物理位置之差。
    far=m["xy"][m["far"]]
    far_mid=(far[:-1]+far[1:])/2
    biases={}
    for key,interp in (("atCentroid",sampler(xy,tris,m["centroid"])),("atFar",boundary_sampler(xy,tris,far_mid))):
        offset=interp(control)-m[key](coarse)
        biases[key]=float(abs(offset).max())
        m[key]=lambda field,fn=interp,bias=offset:fn(field)+bias
    m["atNode"]=boundary_sampler(xy,tris,m["xy"])
    m["seedVelocity"]=velocity;m["flowIds"]=ids;m["flowWallCount"]=width
    m["nativeWall"]=coords[0];m["knotsFlow"]=knots
    m["wallObservation"]=original_face_map(m["corners"],m["diameter"],m["wall"],coords[0],
                                          np.arange(SECTORS*width).reshape(SECTORS,width))
    centroid_error=np.max(abs(m["wallObservation"]@station_xyz.reshape(-1,3)-m["corners"].mean(axis=1)))
    assert centroid_error<1e-8
    m["nearInterpolator"]=sampler(m["xy"],m["tri"],coords[1])
    m["nearFaceMap"]=transfer(m["corners"],m["wall"],knots,m["diameter"])
    return m,dict(controlInterpolationMaxDifferenceOverU=biases,affineCoordinateMaxDifferenceM=float(centroid_error))


def wall_pressure(m,native):
    wall=m["nativeWall"];base_end=np.flatnonzero(np.diff(wall[:,0])>0)[0]
    values=native[m["flowIds"][:,:len(wall)],4]
    p=m["wallObservation"]@values.ravel()
    radii=wall[:base_end+1,1]/wall[base_end,1]
    taps=np.stack([np.interp([0.,.286,.572,.858],radii,v[:base_end+1]) for v in values])
    return 101325+1.177*.05**2*p,taps*1.177*.05**2*1000


def near_temperature(m,theta,prepared):
    mid=m["corners"].mean(axis=1)
    sector=np.floor(np.mod(np.arctan2(mid[:,2],mid[:,1]),2*np.pi)/(2*np.pi/SECTORS)).astype(int)
    out=np.empty(len(mid))
    for j in range(SECTORS):
        mask=sector==j
        out[mask]=300+.001*(m["nearFaceMap"][mask]@m["nearInterpolator"](theta[j]))
    return out


def preflight(work):
    prepared=work/"prepared";m,info=attach(prepared)
    reference=build_model(UP,OLD,SEED/"steady")
    prev=np.loadtxt(OLD/"control.csv",delimiter=",",skiprows=1)
    current=np.loadtxt(prepared/"control.csv",delimiter=",",skiprows=1)
    a,pa=surface_pressure(reference,prev);b,pb=wall_pressure(m,current)
    info["initialOriginalFacePressureMaxDifferencePa"]=float(abs(a-b).max())
    info["initialBaseTapMaxDifferenceMilliPa"]=float(abs(pa-pb).max())
    assert abs(a-b).max()<3e-11 and abs(pa-pb).max()<1e-12
    save(work/"sampling.json",info)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path);parser.add_argument("stage",choices=("sampling","1","2"))
    parser.add_argument("--worker",action="store_true",help=argparse.SUPPRESS)
    args=parser.parse_args();work=args.work.resolve();work.relative_to(BASE)
    if args.stage=="sampling":
        preflight(work)
    else:
        n=int(args.stage);name=f"disturbed_half_{n}";prepared=work/"prepared"
        for j in (1,2):assert read(work/"flow"/f"disturbed_half_{j}"/"step.json")["accepted"]
        if not args.worker:
            step=work/"flow"/name
            paths=[DOC/"PLAN.md",Path(__file__).resolve(),DOC/"wall_sampling.py",DOC.parent/"v354/heat_step.py",
                   DOC.parent/"v354/flux_transfer.py",ROOT/"scripts/simpleLauncherAngularHeat.py",prepared/"refinement.json",
                   prepared/"meridian.npz",prepared/"stations.npz",prepared/"control.csv",work/"sampling.json",
                   step/"step.json",step/f"restart_{n:05d}.csv",UP/"coarse/restart.csv",OLD/"control.csv",OLD/"stations.npz",
                   SEED/"steady/temperature.npz",SEED/"initial/fields.npz",SEED/"initial/results.json"]
            if n==2:paths.extend([work/"heat/disturbed_half_1/thermal.npz",work/"heat/disturbed_half_1/heat.json"])
            target=work/"heat"/name;target.mkdir(parents=True,exist_ok=False)
            rec=observe([sys.executable,str(Path(__file__).resolve()),str(work),str(n),"--worker"],target,
                        {str(p):sha(p) for p in paths},240,300*1024**2)
            assert rec["status"]=="terminal" and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
            assert rec["debuggerExitCode"]==0 and read(target/"heat.json")["accepted"]
            sys.exit(0)
        m,_=attach(prepared)
        sys.exit(0 if run_on_mesh(m,work/"flow",prepared,work/"heat",name,SEED/"initial/fields.npz",5,
                                 wall_pressure,near_temperature,True) else 2)
