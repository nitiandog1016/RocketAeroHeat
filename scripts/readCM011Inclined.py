"""读取 CM-011 +10° 半域原生保存帧，输出完整原面分布和同一圆周的实际三维节点值。"""
import argparse
import json
from pathlib import Path

import numpy as np

from cm011Inclined import ALPHA, BODY_SHA, FLOW, RAYS, config
from cm011Tunnel import CP, D, DT, P0, RHO, T0, U, meridian, write_json
from readCM011TunnelSector import checked_process, frame_window, read_energy, wall_fields
from runEmbeddedSU2 import sha

SCHEMA = "cm011-inclined-half-energy-v1"
MESH_SHA = "afa998bd0f0496153a70f49d44abb92dc3580472ce27088a18de4da1cc0e8ace"
COLD_SHA = "c53c990a74166997b730fe967271f82f38312a0bfdd89edd514b530b50248fc6"


def state_at_frame(folder):
    rec = checked_process(folder)
    for name,expected in (("body.json",BODY_SHA),("mesh.su2",MESH_SHA)):
        if sha(folder/name)!=expected or rec["inputSHA256"].get(str(folder/name))!=expected:
            raise ValueError("三维工况的网格或原面不匹配")
    geom = json.loads((folder/"geometry.json").read_bytes())
    if (geom["alphaDeg"]!=ALPHA or geom["betaDeg"]!=0 or geom["reflection"]!="Y=0" or geom["copies"]!=2 or
            geom["originalBodySHA256"]!=BODY_SHA or not np.array_equal(geom["velocityMps"],U*FLOW) or
            not np.array_equal(geom["rayRadians"],RAYS)):
        raise ValueError("半域边界或来流方向不匹配")
    if (sha(folder/"case.cfg")!="f93b190bf00a68979c02d56c59ccd2fdb51a3856356ad0c00971f79bd5660b34" or
            (folder/"case.cfg").read_text("ascii")!=config(12000) or sha(folder/"initial_flow.csv")!=COLD_SHA):
        raise ValueError("当前运行与登记的 +10° 配置或初态不匹配")
    return frame_window(folder,rec)


def evaluate(folder):
    frame = state_at_frame(folder)
    body = json.loads((folder/"body.json").read_bytes())
    geom = json.loads((folder/"geometry.json").read_bytes())
    coarse,_ = meridian(1)
    fine,_ = meridian(2)
    parsed = wall_fields(folder,frame["savedIndex"],len(coarse),RAYS)
    p,q,temp,near,normal,no_slip,chord_area,chord_power,state = parsed
    assert np.isclose(geom["chordWallAreaM2"],chord_area,rtol=1e-12,atol=0)
    # 方位角从 +Y 转向 +Z。镜像仅翻转 Y，保留 Z 正负和迎背流差异。
    v = np.asarray(body["V"])
    f = np.asarray(body["F"])-1
    row = np.asarray(body["nodeRow"])
    np.testing.assert_allclose(fine[row],np.column_stack((v[:,0],np.hypot(v[:,1],v[:,2])))/D,rtol=0,atol=1e-12)
    phi = np.arctan2(v[:,2],abs(v[:,1]))
    j = np.minimum(np.searchsorted(RAYS,phi,side="right")-1,len(RAYS)-2)
    assert np.min(j)>=0
    azimuth = (phi-RAYS[j])/(RAYS[j+1]-RAYS[j])
    i = np.minimum(row//2,len(coarse)-2)
    meridional = .5*row-i
    assert np.all((meridional>=0)&(meridional<=1)) and np.all((azimuth>=0)&(azimuth<=1+1e-13))

    def map_to_faces(field):
        a = (1-meridional)*field[j,i]+meridional*field[j,i+1]
        b = (1-meridional)*field[j+1,i]+meridional*field[j+1,i+1]
        return ((1-azimuth)*a+azimuth*b)[f].mean(axis=1)

    gauge,heat,near_k,wall_k = [map_to_faces(field) for field in (p,q,T0+DT*(near-1),T0+DT*(temp-1))]
    pts = v[f]
    area = .5*np.linalg.norm(np.cross(pts[:,1]-pts[:,0],pts[:,2]-pts[:,0]),axis=1)
    valid = (frame["normalTerminal"] and frame["passesResidualAndHeatWindow"] and max(normal)<1e-8 and
             no_slip<1e-8 and np.all(q<0) and np.all((near>1)&(near<2)))
    # 预定截面恰为真实求解站。保留 13 条射线的节点值，原三角面常值另由工作台取截面。
    station = 14
    curve = dict(coarseRowZeroBased=station,fineRowZeroBased=2*station,xM=D*coarse[station,0],
                 radiusM=D*coarse[station,1],azimuthDeg=np.rad2deg(RAYS).tolist(),
                 gaugePressurePa=p[:,station].tolist(),heatFluxWm2=q[:,station].tolist(),
                 windwardMinusLeewardPressurePa=p[0,station]-p[-1,station],
                 windwardMinusLeewardHeatFluxWm2=q[0,station]-q[-1,station])
    regions = []
    for k,name in enumerate(body["regionNames"]):
        use = np.asarray(body["regionCode"])==k
        regions.append(dict(name=name,faces=int(use.sum()),areaM2=area[use].sum(),powerIntoWallW=area[use]@heat[use],
                            gaugePressureRangePa=[min(gauge[use]),max(gauge[use])],
                            heatFluxRangeWm2=[min(heat[use]),max(heat[use])]))
    stats = dict(frame,grid="inclined",alphaDeg=ALPHA,betaDeg=0,velocityMps=geom["velocityMps"],
                 originalFaces=len(f),originalBodySHA256=BODY_SHA,areaM2=area.sum(),triangularPowerW=area@heat,
                 fullChordPowerW=2*chord_power,fullWallMonitorW=2*frame["historySectorWallMonitor"]*RHO*U*D**2*CP*DT,
                 maxWallSpeedOverU=no_slip,maxMirrorNormalSpeedsOverU=normal,coneRing=curve,byRegion=regions,
                 pressureGaugeRangePa=[min(gauge),max(gauge)],heatFluxRangeWm2=[min(heat),max(heat)],
                 firstNodeRangeK=[min(near_k),max(near_k)],wallNodeRangeK=[min(wall_k),max(wall_k)],
                 volumeRangeK=[T0+DT*(h+297.15) for h in (state[:,8].min(),state[:,8].max())],
                 acceptedNumericalWall=bool(valid),numericalWallAccepted=bool(valid),physicalAccuracyValidated=False,
                 sourceVTKSHA256=[frame["frameSHA256"][f'surface_{frame["savedIndex"]:06d}.vtk']])
    return stats,np.column_stack((np.arange(len(f)),P0+gauge,heat,near_k,wall_k,gauge))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("new_output",type=Path)
    parser.add_argument("--native-balance",type=Path)
    args=parser.parse_args()
    stats,rows=evaluate(args.folder.resolve())
    if args.native_balance:
        stats["nativeEnergyBalance"]=read_energy(args.native_balance,stats,copies=2,symmetry=("mirror",),schema=SCHEMA)
    args.new_output.mkdir(parents=True,exist_ok=False)
    write_json(args.new_output/"summary.json",stats)
    if stats["acceptedNumericalWall"]:
        np.savetxt(args.new_output/"wall.csv",rows,delimiter=",",comments="",
                   header="faceID,pressurePa,heatFluxWm2,nearWallCellTemperatureK,wallFluidTemperatureK,gaugePressurePa",
                   fmt=["%d"]+["%.17g"]*5)
