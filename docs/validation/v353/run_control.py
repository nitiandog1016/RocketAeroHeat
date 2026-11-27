"""对等温壁及两个方位的局部热边界求解，按源帧、全表面和守恒门槛验收。"""
import argparse
import json
import os
from pathlib import Path
import sys
from time import perf_counter
from zipfile import ZipFile

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from simpleLauncherAngularHeat import ANGLE,PHI,SECTORS,domain,steady_fields,heat_rows,wall_results
from simpleLauncherTunnel import CASE,CP,RHO,SPEED,T_AIR,T_WALL
from runEmbeddedSU2 import sha


def reference_annuli(model, reference_row, area):
    """由 v350 的二维壁反力和完整面积恢复径向线性通量，在真实半径上独立作环带积分。"""
    wet_xy=model["xy"][model["wet"]]
    bottom=wet_xy[:,0]==0
    assert np.array_equal(np.flatnonzero(bottom),np.arange(bottom.sum()))
    r=wet_xy[bottom,1]/wet_xy[bottom,1].max()
    dual=np.asarray(model["faceMap"].T@area).ravel()
    scale=2*np.pi*RHO*SPEED*CP*(T_WALL-T_AIR)*model["diameter"]**2
    q=reference_row[bottom]*scale/dual[bottom]
    face_base=(model["corners"][:,:,0]==0).all(axis=1)
    annuli=[]
    for a,b in zip([0,.25,.5,.75],[.25,.5,.75,1.]):
        knots=np.r_[a,r[(r>a)&(r<b)],b]
        left,right=knots[:-1],knots[1:]
        middle=(left+right)/2
        integrand=lambda s: 2*s*np.interp(s,r,q)
        integral=np.sum((right-left)*(integrand(left)+4*integrand(middle)+integrand(right))/6)
        annuli.append(area[face_base].sum()*integral)
    return np.array(annuli)


def run(job,grid,out):
    start=perf_counter()
    out.relative_to(ROOT/"outputs")
    if not out.is_dir():
        out.mkdir(parents=True)
    assert not (out/"fields.npz").exists() and not (out/"results.json").exists()
    seal=json.loads((CASE/"receipt.json").read_bytes())
    assert sha(CASE/"wind_tunnel.zip")==seal["sha256"]
    with ZipFile(CASE/"wind_tunnel.zip") as z: released=json.loads(z.read("summary.json"))
    assert released["accepted"] and released["condition"]==json.loads((CASE/"tunnel.json").read_bytes())
    reference=json.loads((DOC.parent/"v352/baseline.json").read_bytes())[grid]
    assert released["runs"][grid]["restartSHA256"]==reference["flowSHA256"]
    heat=job/("heat_"+grid)
    assert sha(heat/"wall.csv")==released["runs"][grid]["wallSHA256"]
    primary=np.loadtxt(heat/"wall.csv",delimiter=",",skiprows=1)
    model=domain(job,grid,reference)
    np.testing.assert_array_equal(primary[:,0],np.arange(1,len(primary)+1))
    tri=model["corners"]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    np.testing.assert_allclose(primary[:,1],area,rtol=0,atol=1e-15)
    with np.load(heat/"temperature.npz") as old:
        annuli=reference_annuli(model,old["reaction"],area)
        reference_theta=old["theta"]
    regions=model["definition"]["regions"]
    cylinder=next(r for r in regions if r["name"]=="cylinder")
    xs=tri[cylinder["firstFace"]-1:cylinder["lastFace"],:,0]
    x=model["xy"][model["wet"],0]*model["diameter"]
    envelope=np.zeros(len(x))
    active=(x>xs.min())&(x<xs.max())
    envelope[active]=np.sin(np.pi*(x[active]-xs.min())/np.ptp(xs))**2
    uniform=np.ones((SECTORS,len(x)))
    pattern=uniform+.1*np.cos(PHI[:,None]-.37)*envelope
    walls=np.stack((uniform,pattern,np.roll(pattern,SECTORS//4,axis=0)))
    theta=steady_fields(model,walls)
    wet=model["wet"];axis=model["axis"]
    np.testing.assert_allclose(theta[:,:,wet],walls,rtol=0,atol=1e-12)
    q=[]; reactions=[]; summaries=[]
    for field in theta:
        row,budget=heat_rows(model,field)
        flux,power=wall_results(model,row)
        q.append(flux);reactions.append(row);summaries.append(dict(**budget,**power))
    ref_power=-area@primary[:,3]
    ref_regions={}
    for region in regions:
        faces=slice(region["firstFace"]-1,region["lastFace"])
        ref_regions[region["name"]]=float(-area[faces]@primary[faces,3])
    base=summaries[0]
    region_error={k:abs(base["regionW"][k]/v-1) for k,v in ref_regions.items()}
    ring_error=abs(np.array(base["baseAnnulusW"])/annuli-1)
    rms=np.sqrt((area@(q[0]-primary[:,3])**2)/(area@primary[:,3]**2))
    rotate_temp=float(np.max(abs(theta[2]-np.roll(theta[1],SECTORS//4,axis=0))))
    rotate_wall=float(np.max(abs(reactions[2]-np.roll(reactions[1],SECTORS//4,axis=0)))/np.max(abs(reactions[1])))
    # 独立能量恒等式：取笛卡尔坐标 y=r cosφ，已知连续角向耗散积分为 π∫r dA/Pe。
    y=np.cos(PHI[:,None])*model["xy"][:,1]
    dissipation=ANGLE*np.sum(model["conductance"]*(y-np.roll(y,-1,axis=0))**2)
    expected=np.pi*model["volumePerRadian"]/model["pe"]*np.sinc(1/SECTORS)**2
    free=np.ones(theta.shape[-1],bool);free[np.union1d(wet,model["inlet"])]=False
    interior_range=float(np.max(np.ptp(theta[1][:,free],axis=0)))
    pass_bounds=all(s["minTheta"]>=-1e-12 and s["maxTheta"]<=w.max()+1e-12
                    for s,w in zip(summaries,walls))
    gates=dict(temperatureBounds=pass_bounds,allFaceValuesFinite=bool(np.isfinite(q).all()),
        axisSingleValued=all(s["maxAxisSpread"]<1e-12 for s in summaries),
        residual=all(s["relativeFreeRowResidual"]<1e-8 for s in summaries),
        heatBalance=all(s["energyImbalance"]<.005 for s in summaries),
        powerDifference=bool(abs(base["outwardW"]/ref_power-1)<.01),
        regionalPowers=bool(max(region_error.values())<.02),
        annularPowers=bool(max(ring_error)<.02),surfaceFluxRMS=bool(rms<.02),
        interiorResponds=bool(interior_range>1e-6 and summaries[1]["maxFreeAngularExchangePerRadian"]>1e-12),
        rotation=bool(rotate_temp<1e-9 and rotate_wall<1e-8),angularDissipation=bool(abs(dissipation/expected-1)<1e-10))
    np.savez_compressed(out/"fields.npz",theta=theta,wallReactions=np.array(reactions),wet=wet,axis=axis,
                        xy=model["xy"],triangles=model["tri"],sectorCentreRad=PHI,wallTheta=walls)
    np.savetxt(out/"wall.csv",np.column_stack((primary[:,:3],np.array(q).T)),delimiter=",",comments="",
        fmt=["%d"]+["%.17g"]*5,header="faceID,areaM2,pressurePa,uniformHeatWm2,patternHeatWm2,rotatedHeatWm2")
    summary=dict(grid=grid,accepted=all(gates.values()),gates=gates,
        thermalDOF=SECTORS*(len(model["xy"])-len(axis))+len(axis),originalFaces=len(area),areaM2=area.sum(),
        sectorCount=SECTORS,sourceFlowSHA256=reference["flowSHA256"],sourceHeatMeshSHA256=reference["heatMeshSHA256"],
        sourceFaceSHA256=sha(heat/"wall.csv"),wallSHA256=sha(out/"wall.csv"),fieldsSHA256=sha(out/"fields.npz"),
        flowProjectionError=model["projection"],uniform=base,pattern=summaries[1],rotated=summaries[2],
        regionRelativeDifferences=region_error,referenceBaseAnnulusW=annuli.tolist(),
        baseAnnulusRelativeDifferences=ring_error.tolist(),relativeWallFluxRms=rms,referenceOutwardW=ref_power,
        uniformMaxThetaDifference=float(np.max(abs(theta[0]-reference_theta))),patternInteriorThetaRange=interior_range,
        quarterTurnMaxThetaError=rotate_temp,quarterTurnRelativeReactionError=rotate_wall,
        angularDissipationRelativeError=float(abs(dissipation/expected-1)),seconds=perf_counter()-start)
    (out/"results.json").write_bytes((json.dumps(summary,indent=2,allow_nan=False)+"\n").encode())
    print(json.dumps({"grid":grid,"accepted":summary["accepted"],"gates":gates,"seconds":summary["seconds"]}))
    return summary["accepted"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job",type=Path);parser.add_argument("grid",choices=("coarse","fine","outer"))
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    sys.exit(0 if run(args.job.resolve(),args.grid,args.out.resolve()) else 1)
