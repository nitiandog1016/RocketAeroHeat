# -*- coding: utf-8 -*-
"""一次只读定位收敛场的壁邻熵/温度偏差；不移动采样面授资格、不运行CFD。"""
from pathlib import Path
import argparse
import csv
import json
import shutil
import time
import numpy as np
from readEilmerNative import read_field,sha
from eilmerAirThermo import evaluate


def inspect(base,out):
    base,out=map(Path,(base,out));out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    shutil.copyfile(__file__,out/'executed_inspect.py')
    reference=json.loads((base.parent/'v191_validation/audit01/reference.json').read_text())
    rows=[];summary=[];bindings={}
    for name,nx,ny,nz,snapshot in [('coarse',60,2,20,3),('fine',120,2,40,8)]:
        case=base/(name+'01');rawpath=case/'lmrsim/snapshots'/('%04d'%snapshot)/'fluid-0000';bindings[str(rawpath)]=sha(rawpath)
        f=read_field(case,snapshot,nx*ny*nz);thermo=evaluate(f['p'],f['T'])
        for station,ref in [(.45,reference['states'][1]),(.9,reference['states'][2])]:
            i=min(nx-1,int(station*nx));ids=(i*ny)*nz+np.arange(nz)
            s=-f['pos.x'][ids];wall=np.clip(s-.2,0,.4)*np.tan(np.radians(10));offset=f['pos.z'][ids]-wall
            # 保留整条剖面，再单列近壁固定1cm/3cm/5cm；不择优挑一层当外缘。
            ent=thermo['entropyJkgK'][ids]-ref['entropyJkgK'];terr=f['T'][ids]/ref['temperatureK']-1
            perr=f['p'][ids]/ref['pressurePa']-1
            for k,id_ in enumerate(ids):
                rows.append([name,station,int(id_),float(s[k]),float(offset[k]),float(f['p'][id_]),float(f['T'][id_]),
                    float(ent[k]),float(perr[k]),float(terr[k]),float(f['vel.x'][id_]),float(f['vel.z'][id_])])
            selected=[]
            for limit in [.01,.03,.05]:
                mask=offset<=limit
                selected.append(dict(maxWallOffsetM=limit,count=int(mask.sum()),
                    meanTemperatureError=float(np.mean(terr[mask])) if np.any(mask) else None,
                    meanEntropyExcessJkgK=float(np.mean(ent[mask])) if np.any(mask) else None))
            summary.append(dict(grid=name,stationM=station,actualStreamwiseRangeM=[float(s.min()),float(s.max())],
                nearestCellOffsetM=float(offset[0]),nearestTemperatureError=float(terr[0]),
                nearestPressureError=float(perr[0]),nearestEntropyExcessJkgK=float(ent[0]),fixedBands=selected))
    with (out/'normal_profiles.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['grid','station_m','native_cell','s_m','wall_vertical_offset_m','p_Pa','T_K',
            'entropy_minus_reference_JkgK','relative_pressure_error','relative_temperature_error','Ux_mps','Uz_mps']);w.writerows(rows)
    result=dict(passed=True,profileSummary=summary,sourceSHA256=bindings,
                coordinateMeaning='vertical z minus local wall z; not shortest normal distance on ramp',
                causality='localization only; not a unique cause or permission to bypass near-wall values',
                physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();inspect(a.base,a.out)
