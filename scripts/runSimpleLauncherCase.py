"""单个低速风洞工况的冷启动、三网核对与全表面结果打包。"""
import argparse
import json
import os
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED
import uuid

for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[variable]="1"

import numpy as np

from simpleLauncherTemperature import evaluate
from simpleLauncherTunnel import CASE, CP, K, MU, P_REF, RHO, ROOT, SPEED, T_AIR, T_WALL, prepare
from thermistorBead import read, sha, solve


def compare(a,b):
    left,x=a;right,y=b
    area=y[:,1];dp=x[:,2]-y[:,2];dq=x[:,3]-y[:,3]
    return dict(powerDifference=abs(left['outwardW']/right['outwardW']-1),
                pressureRmsOverDynamic=float(np.sqrt(np.average(dp**2,weights=area))/(.5*RHO*SPEED**2)),
                relativeFluxRms=float(np.sqrt((area@dq**2)/(area@y[:,3]**2))),
                peakOutwardFluxDifference=float(abs(x[:,3].min()/y[:,3].min()-1)),
                regionPowerDifferences={k:abs(v['outwardW']/right['regions'][k]['outwardW']-1)
                                        for k,v in left['regions'].items()})


def seal(flow_dirs,heat_dirs,cross_heat,out):
    """只接收同一原几何及固定工况下经过独立回读、守恒和空间对照的面表。"""
    condition=read(CASE/'tunnel.json')
    assert condition['velocityMps']==[-SPEED,0,0] and condition['ambientPressurePa']==P_REF
    for key,value in [('ambientTemperatureK',T_AIR),('wallTemperatureK',T_WALL),('densityKgM3',RHO),
                      ('dynamicViscosityPaS',MU),('specificHeatJkgK',CP),('conductivityWmK',K)]:
        assert condition[key]==value
    geometry=read(CASE/'geometry.json')
    assert sha(CASE/geometry['file'])==geometry['fileSHA256']==condition['surfaceSHA256']
    frames=[];runs={}
    def checked(flow,heat,grid,spacing):
        s=read(heat/'summary.json');r=read(flow/'engine_process.json');t=s['thermal']
        assert r['status']=='terminal' and r['exitCode']==0 and not(r['cancelled'] or r['timedOut']) and all(r['inputsStable'].values())
        assert s['grid']==grid and s['surfaceSHA256']==geometry['fileSHA256'] and s['heatSpacingD']==spacing
        assert sha(flow/'restart.csv')==r['outputSHA256']['restart.csv']==s['sourceRestartSHA256']
        assert s['wallSHA256']==sha(heat/'wall.csv') and max(s['residuals'].values())< -8
        assert t['temperatureMin']>=-1e-12 and t['temperatureMax']<=1+1e-12 and t['linearResidual']<1e-8
        assert s['heatStability']<.002 and t['balance']<.005
        wall=np.loadtxt(heat/'wall.csv',delimiter=',',skiprows=1)
        assert np.array_equal(wall[:,0],np.arange(1,geometry['faces']+1))
        assert np.isfinite(wall).all() and min(wall[:,2])>0 and min(wall[:,1])>0
        assert abs(wall[:,1].sum()-geometry['areaM2'])<1e-10 and abs(wall[:,1]@wall[:,3]+s['outwardW'])<1e-11
        return s,wall
    for grid,spacing,flow,heat in zip(('coarse','fine','outer'),(.004,.002,.002),flow_dirs,heat_dirs):
        s,wall=checked(flow,heat,grid,spacing)
        frames.append((s,wall))
        receipt=read(flow/'engine_process.json')
        runs[grid]=dict(step=s['finalStep'],flowSeconds=receipt['elapsedSeconds'],rmsMaximum=s['residuals'],
                        wallSHA256=s['wallSHA256'],restartSHA256=s['sourceRestartSHA256'],
                        heatNodes=s['heatNodes'],heatStability=s['heatStability'],thermal=s['thermal'],
                        areaM2=s['areaM2'],outwardW=s['outwardW'])
    cross=checked(flow_dirs[1],cross_heat,'fine',.004)
    grids=compare(frames[0],frames[1]);far=compare(frames[1],frames[2])
    thermal_only=compare(cross,frames[1]);flow_only=compare(frames[0],cross)
    if any(test['powerDifference']>=.05 for test in (grids,far,thermal_only,flow_only)):
        raise ValueError("网格或外域的热功率变化达到 5%，不发布最终结果")
    main,wall=frames[-1]
    for _,other in frames:np.testing.assert_array_equal(wall[:,0:2],other[:,0:2])
    s=dict(accepted=True,condition=condition,surfaceSHA256=geometry['fileSHA256'],wallSHA256=main['wallSHA256'],
           originalFaces=len(wall),areaM2=main['areaM2'],outwardW=main['outwardW'],regions=main['regions'],
           pressureRangePa=[wall[:,2].min(),wall[:,2].max()],heatFluxRangeWm2=[wall[:,3].min(),wall[:,3].max()],
           grids=grids,outerDomain=far,fixedFlowHeatGrid=thermal_only,fixedHeatSpacingFlowGrid=flow_only,
           maxOutwardFluxFaceID=int(wall[np.argmin(wall[:,3]),0]),runs=runs,verifiedByExperimentM2=0,
           method='SU2 8.6.0 不可压恒密度层流动量；散度自由流函数投影、Gmsh 轴对称热网、守恒低阶上风有限元及强等温壁反力')
    out.mkdir(parents=True,exist_ok=False)
    with ZipFile(out/'wind_tunnel.zip','x',ZIP_DEFLATED) as z:
        z.writestr('summary.json',json.dumps(s,ensure_ascii=False,allow_nan=False,indent=2))
        z.write(heat_dirs[-1]/'wall.csv','wall.csv')
    (out/'receipt.json').write_text(json.dumps(dict(sha256=sha(out/'wind_tunnel.zip'),faces=len(wall),
              outwardW=main['outwardW']),indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(faces=len(wall),outwardW=main['outwardW'],meshError=grids['powerDifference'],farError=far['powerDifference'])))


def fresh(job):
    job=job.resolve();job.relative_to(ROOT/'outputs');job.mkdir(parents=True,exist_ok=False)
    flows=[];heats=[]
    for grid,spacing in [('coarse',.004),('fine',.002),('outer',.002)]:
        flow=job/grid;heat=job/('heat_'+grid)
        prepare(flow,grid)
        request=read(flow/'engine_request.json');request['jobToken']=str(uuid.uuid4())
        (flow/'engine_request.json').write_text(json.dumps(request),encoding='utf8')
        solve(flow,memory_bytes=4*1024**3,output_bytes=2*1024**3)
        evaluate(flow,heat,spacing);flows.append(flow);heats.append(heat)
        if grid=='fine':evaluate(flow,job/'cross_heat',.004)
    seal(flows,heats,job/'cross_heat',job/'result')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('job',type=Path)
    args=parser.parse_args();fresh(args.job)
