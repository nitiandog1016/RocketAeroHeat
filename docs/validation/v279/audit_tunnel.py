"""从实际终场重建圆柱专属热通量、校验重启列及同阶段功率，再映射到 2160 个原面。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens
from runEmbeddedSU2 import sha


def read(path):
    return json.loads(path.read_bytes())


def table(path):
    with path.open() as f:
        header=[s.strip() for s in next(csv.reader(f,skipinitialspace=True))]
    return header,np.loadtxt(path,delimiter=",",skiprows=1,ndmin=2)


def project(x,wall_values):
    """在原三角面的仿射 (周向,跨向) 坐标中，对原生格点双线性场分片精确积分。"""
    result=[]
    degree=np.arange(361)
    # 原生壁面每 5° 一站，拼接完整圆周。
    value=np.array([np.interp(degree,np.arange(0,361,5),np.r_[row,row[0]]) for row in wall_values])
    original=np.array([-.076,-.033,.033,.076])
    for left,right in zip(original[:-1],original[1:]):
        bounds=np.r_[left,x[(x>left)&(x<right)],right]
        first=np.zeros(360);second=np.zeros(360)
        for lo,hi in zip(bounds[:-1],bounds[1:]):
            # 对给定跨向位置，先沿周向积分两个互补三角形，再用两点 Gauss 积分三次式。
            for sign in (-1,1):
                point=(lo+hi)/2+sign*(hi-lo)/(2*np.sqrt(3))
                k=np.searchsorted(x,point)-1
                fraction=(point-x[k])/(x[k+1]-x[k])
                line=(1-fraction)*value[k]+fraction*value[k+1]
                base=line[:-1];slope=np.diff(line);t=(point-left)/(right-left)
                weight=(hi-lo)/(right-left)
                first+=weight*(base*(1-t)+slope*(1-t*t)/2)
                second+=weight*(base*t+slope*t*t/2)
        result.extend((first,second))
    return np.concatenate(result)


def audit(base,out):
    run=base/"flow04"
    rec=read(run/"terminal.json")
    assert rec["exitCode"]==0 and rec["stopReason"] is None and rec["inputsStable"]
    assert all(sha(run/name)==h for name,h in rec["outputs"].items())
    c=read(run/"condition.json")
    histnames,hist=table(run/"history_00001.csv")
    assert np.array_equal(hist[:,2],np.arange(120)) and np.all(hist[:,0]==1) and np.isfinite(hist).all()
    keys,volume=table(run/"restart_00001.csv")
    oldkeys,old=table(base/"flow02/restart_00001.csv")
    seed_names,seed=table(run/"seed_00000.csv")
    assert oldkeys==seed_names[:13] and np.array_equal(volume[:,[keys.index(k) for k in oldkeys]],old)
    assert np.array_equal(volume[:,0],seed[:,0]) and np.isfinite(volume).all()
    np.testing.assert_allclose(volume[:,1:4],seed[:,1:4],rtol=0,atol=1e-15)
    assert sha(run/"history_00001.csv")==sha(base/"flow02/history_00001.csv")
    assert sha(run/"cylinder_surface_00001.vtk")==sha(base/"flow02/cylinder_surface_00001.vtk")
    assert len(volume)==73800 and set(seed_names).issubset(keys)
    # SU2 转捩读取器在两个主量之后顺序读取两列附加量；按列名组装，检查逐值往返。
    out.mkdir(parents=True,exist_ok=False)
    target=out/"restart_with_lm_00001.csv"
    ordered=volume[:,[keys.index(k) for k in seed_names]]
    np.savetxt(target,ordered,delimiter=",",comments="",fmt=["%d"]+["%.17g"]*(len(seed_names)-1),
               header=",".join('"'+k+'"' for k in seed_names))
    double_names,double_values=table(target)
    assert double_names==seed_names and np.array_equal(double_values,ordered)
    wall=np.load(run/"wall.npz");ids=wall["nodeIDs"]
    x=wall["x"];wall_xyz=volume[ids,1:4]
    np.testing.assert_allclose(wall_xyz,wall["xyz"],rtol=0,atol=1e-15)
    coords,fields=vtk_tokens(run/"cylinder_surface_00001.vtk")
    gap,idx=cKDTree(coords).query(wall_xyz)
    assert len(np.unique(idx))==len(coords)==1800 and gap.max()<1e-7
    # 标记交线同时属于等温圆柱和绝热安装壁，节点通用 Heat_Flux 在交线上被后者覆盖。
    native=volume[ids,keys.index("Heat_Flux")]
    np.testing.assert_allclose(native,fields["Heat_Flux"][:,0][idx],atol=.0001,rtol=0)
    temp=volume[:,keys.index("Temperature")]
    np.testing.assert_allclose(temp,volume[:,keys.index("Enthalpy")]/c["cpJkgK"]+298.15,rtol=0,atol=2e-12)
    inner=ids+72
    delta=volume[inner,1:4]-wall_xyz
    normals=wall_xyz.copy();normals[:,:,0]=0;normals/=np.linalg.norm(normals,axis=2)[:,:,None]
    dist=np.einsum("ijk,ijk->ij",delta,normals)
    np.testing.assert_allclose(dist,8e-6,rtol=0,atol=1e-14)
    flux=c["conductivityWmK"]*(temp[inner]-c["wallTemperatureK"])/dist
    error=float(np.max(abs(flux[1:-1]-native[1:-1])))
    assert error<1e-6 and np.all(native[[0,-1]]==0)
    widths=np.diff(x)
    node_span=(np.r_[0,widths]+np.r_[widths,0])/2
    dual_power=float(node_span@flux.sum(axis=1)*.033*np.sin(2*np.pi/72))
    monitor=hist[-1,histnames.index("HF")]
    assert abs(dual_power-monitor)<1e-7
    # 将出口数值零点转为新研究工况的上游静压；实际入口中心点单独核对。
    ref_id=ids[12,0]+40*72
    np.testing.assert_allclose(volume[ref_id,1:4],[0,.33,0],rtol=0,atol=1e-14)
    relative=volume[ids,keys.index("Pressure")]-volume[ref_id,keys.index("Pressure")]
    pressure=c["referencePressurePa"]+project(x,relative);heat=project(x,flux)
    assert min(pressure)>0 and np.max(abs(volume[ids,5:8]))<1e-10
    original_page=ROOT/"outputs/v278_validation/gui01/workbench/case-data.js"
    model=json.loads(original_page.read_text("utf8").strip().split("=",1)[1][:-1])
    tri=np.asarray(model["vertices"])[model["faces"]]
    angle=np.rad2deg(np.unwrap(np.arctan2(tri[:,:,2],tri[:,:,1]),axis=1))
    sector=np.rint(np.min(angle,axis=1)).astype(int)%360
    assert np.array_equal(sector,np.tile(np.arange(360),6))
    for i,(lo,hi) in enumerate(zip([-.076,-.033,.033],[-.033,.033,.076])):
        band=tri[720*i:720*(i+1)]
        np.testing.assert_allclose(band[:,:,0].min(axis=1),lo,atol=1e-14,rtol=0)
        np.testing.assert_allclose(band[:,:,0].max(axis=1),hi,atol=1e-14,rtol=0)
        np.testing.assert_allclose(band[:360,:,0].mean(axis=1),(2*lo+hi)/3,atol=1e-14,rtol=0)
        np.testing.assert_allclose(band[360:,:,0].mean(axis=1),(lo+2*hi)/3,atol=1e-14,rtol=0)
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    projected=float(area@heat)
    expected=float(node_span@flux.sum(axis=1)*360*2*.033*np.sin(np.pi/360)/72)
    assert abs(projected-expected)<1e-10
    np.savetxt(out/"trial_faces.csv",np.column_stack((np.arange(1,2161),tri.mean(axis=1),area,pressure,heat)),
               fmt=["%d"]+["%.16g"]*6,header="faceID,xM,yM,zM,areaM2,trialPressurePa,trialHeatFluxWm2",comments="",delimiter=",")
    residuals=dict(zip(histnames[3:12],hist[-1,3:12]))
    positive={k:[float(volume[:,keys.index(k)].min()),float(volume[:,keys.index(k)].max())]
              for k in ("Temperature","Turb_Kin_Energy","Omega","LM_gamma","LM_Re_t")}
    assert all(a>0 for a,b in positive.values())
    result=dict(sourceFrameSHA256=sha(run/"restart_00001.csv"),sourceHistorySHA256=sha(run/"history_00001.csv"),
        seedSHA256=sha(run/"seed_00000.csv"),extendedRestartSHA256=sha(target),extendedNames=seed_names,
        primaryReplayIdentical=True,physicalStepsExecuted=1,physicalDtS=.05*.066/c["velocityMps"],
        innerIterations=120,terminalResiduals=residuals,passesNineResiduals=bool(max(residuals.values()) < -8),
        savedVolumeRanges=positive,nonSharedNativeHeatMaxDifferenceWm2=error,sharedNodes=144,
        sharedCylinderHeatRangeWm2=[float(flux[[0,-1]].min()),float(flux[[0,-1]].max())],
        nativeWallPowerW=monitor,rebuiltDualWallPowerW=dual_power,projectedTrialPowerW=projected,
        upstreamReferenceNode=int(ref_id),upstreamNativeGaugePa=float(volume[ref_id,4]),
        trialFaceCount=len(heat),trialPressureRangePa=[float(pressure.min()),float(pressure.max())],
        trialHeatRangeWm2=[float(heat.min()),float(heat.max())],qualifiedAreaM2=0,productionQualified=False,
        originalGeometryPageSHA256=sha(original_page),faceSHA256=sha(out/"trial_faces.csv"),condition=c)
    (out/"trial.json").write_text(json.dumps(result,indent=2),encoding="utf8")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    print(audit(args.base,args.new_output)["terminalResiduals"])
