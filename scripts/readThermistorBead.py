"""从热敏珠本工况的原生终态恢复有符号热流和压力，投影到完整原面。"""
import argparse
from pathlib import Path

import numpy as np

from readCM011TunnelSector import frame_window, wall_fields
from thermistorBead import AREA, CP, D, DT, ENGINE_SHA, GRIDS, K, P0, RHO, T0, U, body, read, sha, write_json
from thermistorBeadSector import COPIES, sector_rays


def frame(folder):
    rec=read(folder/"engine_process.json")
    assert rec["status"]=="terminal" and rec["exitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["timedOut"] or rec["cancelled"] or rec.get("supervisorError"))
    assert rec["inputSHA256"][rec["command"][0]]==ENGINE_SHA and rec["runtimeFlavor"]=="official"
    for name in ("body.json","case.json","mesh.su2","case.cfg","nodes.npz","initial_flow.csv","geometry.json"):
        assert sha(folder/name)==rec["inputSHA256"][str(folder/name)]
    assert not read(folder/"resources.json")["resourceExceeded"]
    return frame_window(folder,rec)


def evaluate(folder):
    f=frame(folder)
    assert f["normalTerminal"] and f["passesResidualAndHeatWindow"]
    shape=read(folder/"body.json")
    assert shape==body()
    spec=read(folder/"case.json")
    grid=spec["grid"]
    m=GRIDS[grid][0]+1
    geom=read(folder/"geometry.json")
    rays=sector_rays(len(geom["raysRadians"])-1)
    assert geom["rotationCopies"]==COPIES and geom["raysRadians"]==rays.tolist()
    gauge,q,t,t1,sym,wall_speed,chord_area,chord_power,volume=wall_fields(folder,f["savedIndex"],m,rays,
                pressure_scale=RHO*U**2,heat_scale=RHO*CP*U*DT,diameter=D,ambient_pressure=P0)
    assert abs(chord_area-geom["chordWallAreaM2"])<1e-16
    # 先按原模型的子午位置插值，再按当前扇区的实际射线取周向值。
    vertices=np.asarray(shape["V"])
    faces=np.asarray(shape["F"])-1
    angles=np.mod(np.arctan2(vertices[:,2],vertices[:,1]),2*rays[-1])
    h=np.minimum(angles,2*rays[-1]-angles)/rays[1]
    lower=np.minimum(np.floor(h).astype(int),len(rays)-2)
    idx=np.asarray(shape["nodeRow"])*(m-1)/64
    def original(field):
        by_ray=np.array([np.interp(idx,np.arange(m),ray) for ray in field])
        ids=np.arange(len(vertices))
        val=by_ray[lower,ids]*(1-h+lower)+by_ray[lower+1,ids]*(h-lower)
        return val[faces].mean(axis=1)
    dp,flux,tw,near=[original(x) for x in (gauge,q,T0+DT*(t-1),T0+DT*(t1-1))]
    tri=vertices[faces]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=.5*np.linalg.norm(cross,axis=1)
    edges=np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    _,inv,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(counts==2) and np.all(np.bincount(inv,weights=np.sign(edges[:,1]-edges[:,0]))==0)
    assert np.all(np.sum(cross*tri.mean(axis=1),axis=1)>0) and abs(area.sum()/AREA-1)<.005
    tempvol=T0+DT*(volume[:,8]+297.15)
    assert np.min(volume[:,4])*RHO*U**2+P0>0 and np.isfinite(volume).all()
    assert max(sym)<1e-8 and wall_speed<1e-8 and np.all(q<0) and np.all((t1>1)&(t1<2))
    monitor=f["historySectorWallMonitor"]*COPIES*RHO*CP*U*DT*D**2
    powers=[float(area@flux),float(COPIES*chord_power),float(monitor)]
    # 不在任一节点裁温。容差只覆盖 ASCII 终场存储与焓到温度的舍入。
    temperature_ok=bool(tempvol.min()>=T0-1e-7 and tempvol.max()<=T0+DT+1e-7)
    stats=dict(grid=grid, sourceFrameSHA256=sha(folder/"restart.csv"), sourceBodySHA256=sha(folder/"body.json"),
               caseSHA256=sha(folder/"case.json"), configSHA256=sha(folder/"case.cfg"), meshSHA256=sha(folder/"mesh.su2"),
               frame=f, engineSeconds=read(folder/"engine_process.json")["elapsedSeconds"],
               savedVolumeTemperatureRangeK=[tempvol.min(),tempvol.max()],temperatureBounded=temperature_ok,
               originalFaces=len(area),areaM2=area.sum(),analyticAreaM2=AREA,outputOrigin="research",
               firstNodeRangeK=[near.min(),near.max()],wallFluidRangeK=[tw.min(),tw.max()],
               pressurePaRange=[P0+dp.min(),P0+dp.max()],heatFluxWm2Range=[flux.min(),flux.max()],
               wallSpeedOverU=wall_speed,maxSymmetryNormalVelocityOverU=max(sym),
               originalWallPowerW=powers[0],chordWallPowerW=powers[1],fullNativeWallMonitorW=powers[2],
               meanCoefficientWm2K=-powers[0]/(area.sum()*DT),meanNuD=-powers[0]*D/(area.sum()*DT*K),
               areaAndWindowAccepted=temperature_ok,independentlyVerifiedAreaM2=0)
    return stats,np.column_stack((np.arange(len(area)),area,P0+dp,dp,flux,tw,near))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    stats,rows=evaluate(args.run.resolve())
    out=args.new_output
    out.mkdir(parents=True,exist_ok=False)
    np.savetxt(out/"wall.csv",rows,fmt=["%d"]+["%.17g"]*6,delimiter=",",comments="",
               header="faceID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,wallFluidTemperatureK,nearWallCellTemperatureK")
    write_json(out/"summary.json",stats)
