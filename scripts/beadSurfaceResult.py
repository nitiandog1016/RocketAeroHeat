"""将热敏珠的原生压力与强温差边界反力映回全珠原三角面。"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from thermistorBead import A,B,CP,D,DT,K,P0,PE,RHO,ROOT,T0,U,body,read,sha,write_json
from thermistorBeadSector import COPIES,RAYS,sector_rays
from thermistorReaction import terminal


def projection(vertices,faces,ids,rays=RAYS):
    """在原三角面内做等面积子三角积分，得到到各原生壁节点的非负插值权重。"""
    tri=vertices[faces]
    wall_ids=np.unique(ids)
    assert ids.shape[0]==len(rays)
    sectors=(rays[-1],rays[1]-rays[0])
    rows=[];cols=[];weights=[]
    subdivisions=6
    for a in range(subdivisions):
        for b in range(subdivisions-a):
            offsets=[1./3] if a+b+1==subdivisions else [1./3,2./3]
            for shift in offsets:
                u,v=(a+shift)/subdivisions,(b+shift)/subdivisions
                x,y,z=((1-u-v)*tri[:,0]+u*tri[:,1]+v*tri[:,2]).T
                mer=np.arctan2(np.hypot(y,z)/B,-x/A)*(ids.shape[1]-1)/np.pi
                phase=np.mod(np.arctan2(z,y),2*sectors[0])
                ray=np.minimum(phase,2*sectors[0]-phase)/sectors[1]
                j=np.minimum(mer.astype(int),ids.shape[1]-2)
                k=np.minimum(ray.astype(int),len(rays)-2)
                for dj,wm in ((0,1-mer+j),(1,mer-j)):
                    for dk,wa in ((0,1-ray+k),(1,ray-k)):
                        rows.append(np.arange(len(tri)))
                        cols.append(np.searchsorted(wall_ids,ids[k+dk,j+dj]))
                        weights.append(wm*wa/subdivisions**2)
    mat=coo_matrix((np.concatenate(weights),(np.concatenate(rows),np.concatenate(cols))),
                   shape=(len(tri),len(wall_ids))).tocsr()
    assert mat.data.min()>=-2e-15
    np.testing.assert_allclose(np.asarray(mat.sum(axis=1)).ravel(),1,atol=3e-14,rtol=0)
    return mat,wall_ids


def thermal_frame(folder):
    rec=terminal(folder)
    with (folder/"history.csv").open() as f:
        reader=csv.DictReader(f,skipinitialspace=True)
        reader.fieldnames=[k.strip() for k in reader.fieldnames]
        history=list(reader)
    last=int(history[-1]["Inner_Iter"])
    assert len(history)>=20 and [int(r["Inner_Iter"]) for r in history[-20:]]==list(range(last-19,last+1))
    keys=["rms["+s+"]" for s in ("P","U","V","W","h","rho*Y_0")]
    window=np.asarray([[float(row[key]) for key in keys] for row in history[-20:]])
    assert np.isfinite(window).all()
    residuals=dict(zip(keys,window.max(axis=0)))
    assert max(residuals.values()) < -8
    assert sha(folder/f"restart_{last:06d}.csv")==sha(folder/"restart.csv")
    return rec,residuals


def state_and_reaction(parent,probe,wall_ids,xyz):
    rec,residuals=thermal_frame(parent)
    terminal(probe)
    assert sha(probe/"initial_flow.csv")==sha(parent/"restart.csv") and sha(probe/"mesh.su2")==sha(parent/"mesh.su2")
    cfg=(parent/"case.cfg").read_text("ascii")
    assert "MARKER_WALL_SPECIES= (wall, VALUE, 1.0)" in cfg and "SPECIES_CLIPPING= NO" in cfg
    assert f"DIFFUSIVITY_CONSTANT= {1/PE:.17g}" in cfg
    native=np.loadtxt(parent/"restart.csv",delimiter=",",skiprows=1)
    assert native.shape==(len(xyz),10) and np.isfinite(native).all()
    assert np.array_equal(native[:,0],np.arange(len(xyz)))
    np.testing.assert_allclose(native[:,1:4],xyz,rtol=0,atol=1e-13)
    theta=native[:,-1]
    assert theta.min()>=-1e-12 and theta.max()<=1+1e-12 and np.all(theta[wall_ids]==1)
    assert np.max(abs(native[wall_ids,5:8]))<1e-10 and min(P0+native[:,4]*RHO*U**2)>0
    with (probe/"restart.csv").open() as f:
        names=next(csv.reader(f))
    raw=np.loadtxt(probe/"restart.csv",delimiter=",",skiprows=1,usecols=(0,names.index("Residual_Species_0")))
    assert np.array_equal(raw[:,0],native[:,0]) and np.isfinite(raw).all()
    reaction=raw[:,1]*RHO*U*CP*DT*D**2*COPIES
    assert np.all(reaction[wall_ids]<0)
    return native,reaction,residuals,rec


def collect(parent,probes,out):
    geometry=read(parent/"geometry.json");spec=read(parent/"case.json");shape=read(parent/"body.json")
    rec=terminal(parent)
    for name in ("geometry.json","body.json","nodes.npz","case.json"):
        assert sha(parent/name)==rec["inputSHA256"][str(parent/name)]
    # 固定三套名义试件体网、量纲、物性和标记。旧的弱壁解不能用于这里的结果。
    contract=read(ROOT/"cases/thermistor_b/grids.json")
    angular_cells=len(geometry["raysRadians"])-1
    rays=sector_rays(angular_cells)
    fixed=contract["grids"][spec["grid"] if angular_cells==2 else "outerAngular4"]
    assert all(sha(parent/name)==digest for name,digest in fixed.items())
    assert sha(parent/"case.cfg")==contract["strongConfigSHA256"]
    continuation=probes/"plus20"
    assert sha(continuation/"mesh.su2")==sha(parent/"mesh.su2")
    main_config=(parent/"case.cfg").read_text("ascii")
    assert (continuation/"case.cfg").read_text("ascii")==main_config.replace("ITER= 36000\n","ITER= 20\n")
    for source,probe in ((parent,probes/"reaction0"),(continuation,probes/"reaction20")):
        config=(source/"case.cfg").read_text("ascii").splitlines()
        changes={"MARKER_WALL_SPECIES": "(wall, FLUX, 0.0)", "ITER": "1",
                 "VOLUME_OUTPUT": "(COORDINATES, SOLUTION, RESIDUAL)"}
        expected=[line.split("=",1)[0]+"= "+changes[line.split("=",1)[0]]
                  if line.split("=",1)[0] in changes else line for line in config]
        assert (probe/"case.cfg").read_text("ascii").splitlines()==expected+["WRT_RESTART_COMPACT= NO"]
    assert shape==body() and geometry["raysRadians"]==rays.tolist() and geometry["rotationCopies"]==COPIES
    with np.load(parent/"nodes.npz") as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    m=spec["meridianIntervals"]+1
    vertices,faces=np.array(shape["V"]),np.array(shape["F"])-1
    transfer,wall_ids=projection(vertices,faces,ids[:,:m],rays)
    tri=vertices[faces];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=np.linalg.norm(cross,axis=1)/2
    assert np.all(np.sum(cross*tri.mean(axis=1),axis=1)>0)
    s,r,res,rec=state_and_reaction(parent,probes/"reaction0",wall_ids,xyz)
    normal_speed=max(np.max(np.abs(s[ids[k],7]*np.cos(rays[k])-s[ids[k],6]*np.sin(rays[k]))) for k in (0,-1))
    assert normal_speed<1e-8
    _,r20,res20,_=state_and_reaction(probes/"plus20",probes/"reaction20",wall_ids,xyz)
    assert sha(probes/"plus20/initial_flow.csv")==sha(parent/"restart.csv")
    # 面积集总的保守传递：每个原生边界节点的反力，分给在该点有权重的显示面。
    dual=transfer.T@area
    assert np.all(dual>0)
    flux=transfer@(r[wall_ids]/dual)
    power=float(area@flux);assert abs(power-r[wall_ids].sum())<1e-13
    delta=abs(r20[wall_ids].sum()/power-1)
    assert delta<.002 and abs((r.sum()-power)/power)<.005
    near_id=np.zeros(len(xyz),int);near_id[ids[:,:m]]=ids[:,m:2*m]
    near=T0+DT*(transfer@s[near_id[wall_ids],-1])
    gauge=RHO*U**2*(transfer@s[wall_ids,4])
    assert np.all(flux<0) and np.all((near>T0)&(near<T0+DT))
    out.mkdir(parents=True,exist_ok=False)
    path=out/"wall.csv"
    np.savetxt(path,np.column_stack((np.arange(len(faces)),area,P0+gauge,gauge,flux,near)),
               fmt=["%d"]+["%.17g"]*5,delimiter=",",comments="",
               header="faceID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,nearWallCellTemperatureK")
    report=dict(grid=spec["grid"],originalFaces=len(area),areaM2=area.sum(),fullSurfaceComputed=True,
                verifiedJointAreaM2=0,meanCoefficientWm2K=-power/(DT*area.sum()),meanNuD=-power*D/(DT*area.sum()*K),
                nominalExperimentAreaCoefficientWm2K=-power/(DT*spec["areaSourceM2"]),powerW=power,
                probeAllResidualW=r.sum(),boundaryRelativeDifference=abs(r.sum()/power-1),
                heatChange20=delta,last20MaxLogRMS=res,continued20MaxLogRMS=res20,
                temperatureRangeK=[T0+DT*s[:,-1].min(),T0+DT*s[:,-1].max()],
                weakEnthalpyTemperatureRangeK=[T0+DT*(s[:,8].min()+297.15),T0+DT*(s[:,8].max()+297.15)],
                pressurePaRange=[(P0+gauge).min(),(P0+gauge).max()],fluxRangeWm2=[flux.min(),flux.max()],
                firstNodeRangeK=[near.min(),near.max()],wallTableSHA256=sha(path),bodySHA256=sha(parent/"body.json"),
                sourceFrameSHA256=sha(parent/"restart.csv"),sourceReactionSHA256=sha(probes/"reaction0/restart.csv"),
                configSHA256=sha(parent/"case.cfg"),elapsedSeconds=rec["elapsedSeconds"])
    if angular_cells==4:
        report.update(angularCells=angular_cells,cellAngleDeg=30./angular_cells,
                      meshSHA256=fixed["mesh.su2"],maxSymmetryNormalSpeedOverU=normal_speed)
    write_json(out/"summary.json",report)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("parent","probes","output"):parser.add_argument(k,type=Path)
    args=parser.parse_args()
    out=args.output.resolve();out.relative_to(ROOT/"outputs")
    collect(args.parent.resolve(),args.probes.resolve(),out)
