"""模块用途：复核三维层流保存态与多marker热流；退出0绝不等同实验通过。"""
import argparse
import csv
from pathlib import Path
import time

import numpy as np

from idealAirProperties import IdealAir
from nacaWallFlux import geometry,evaluate
from su2NACALaminarInterface import mesh,load_state,inventory
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save,paper_cp


def audit(folder,out):
    """模块一：逐帧同阶段能量/EOS、弱热边界与独立面积热功率核对。"""
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_audit.py').write_bytes(Path(__file__).read_bytes())
    receipt=read(folder/'engine_process.json')
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and all(receipt['inputsStable'].values())
    assert not receipt.get('supervisorError') and receipt['checkpoints']['normalScheduleComplete']
    locks={str(folder/name):sha(folder/name) for name in
        ('engine_process.json','mesh.su2','case.cfg','case_manifest.json','geometry.json','history.csv')}
    locks[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
    for name in ('nacaWallFlux.py','su2NACALaminarInterface.py'):
        file=ROOT/'scripts'/name;locks[str(file)]=sha(file)
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    geom=geometry(p,cells,markers);wall_values={r['marker']:r['givenTw_K'] for r in meta['walls']}
    with (folder/'history.csv').open() as f:
        history=list(csv.DictReader(f,skipinitialspace=True))
    core=IdealAir();locks.update(core.locks);frames=[]
    for pair in receipt['checkpoints']['complete']:
        path=folder/pair['volume'];assert sha(path)==pair['volumeSHA256']
        locks[str(path)]=sha(path);a,v=load_state(path);rho=v['Density']
        u=np.column_stack([v['Momentum_'+axis]/rho for axis in 'xyz'])
        thermo=core.state('PT',v['Pressure'],v['Temperature'])
        thermal=v['Energy']/rho-.5*np.sum(u*u,axis=1)
        energy_error=float(max(abs(thermal-thermo['e'])/np.maximum(1,abs(thermo['e']))))
        eos_error=float(max(abs(rho/thermo['rho']-1)))
        cp_error=float(max(abs(paper_cp(v['Temperature'])/thermo['cp']-1)))
        conductivity=thermo['cp']*v['Laminar_Viscosity']/.72
        per_marker=evaluate(geom,v['Temperature'],conductivity,wall_values)
        rows=[];last=np.full(len(p),np.nan);power=0.;counts=np.zeros(len(p),int)
        for tag,group in per_marker.items():
            g=geom[tag];power+=group['powerIntoWall'];ids=group['ids']
            last[ids]=group['qIntoWall'];counts[ids]+=1
            for i,near,distance,area,q in zip(ids,g['normalNeighbor'],g['distance'],g['area'],group['qIntoWall']):
                rows.append([tag,int(i),int(near),*p[i],wall_values[tag],v['Temperature'][i],
                    distance,area,q,area*q,'native-multiple-isothermal-weak-wall-reconstruction',False])
        wall=np.where(counts)[0]
        # 当前CFlowOutput按marker顺序写共享节点的唯一显示列；独立功率始终保留两侧。
        q_error=float(max(abs(last[wall]-v['Heat_Flux'][wall])/np.maximum(1,abs(last[wall]))))
        native_power=float(history[pair['localIteration']]['HF'].strip())
        power_error=abs(power-native_power)/max(1,abs(power))
        filename='wall_marker_%06d.csv'%pair['localIteration']
        with (out/filename).open('w',newline='',encoding='utf8') as f:
            writer=csv.writer(f);writer.writerow(['marker','PointID','normalNeighbor','X_m','Y_m','Z_m',
                'givenTw_K','wallFluidT_K','normalDistance_m','dualArea_m2','qIntoWall_Wm2','powerIntoWall_W','method','physicalQualified'])
            writer.writerows(rows)
        errors=[abs(row[7]-row[6]) for row in rows]
        frame=dict(localIteration=pair['localIteration'],energyMaximumRelativeError=energy_error,
            eosMaximumRelativeError=eos_error,paperCpMaximumRelativeError=cp_error,
            displayedHeatFluxMaximumScaledError=q_error,markerPowerW=power,nativeHistoryPowerW=native_power,
            powerRelativeDifference=power_error,sharedWallNodes=int(sum(counts>1)),
            wallFluidT_RangeK=[float(min(v['Temperature'][wall])),float(max(v['Temperature'][wall]))],
            givenTwMaximumDifferenceK=max(errors),yPlusMaximum=float(max(v['Y_Plus'][wall])),
            qIntoWallRangeWm2=[min(r[10] for r in rows),max(r[10] for r in rows)],
            wallCSV=filename,wallCSV_SHA256=sha(out/filename),
            interfacePassed=bool(max(energy_error,eos_error)<2e-11 and cp_error<2e-8 and q_error<2e-10 and power_error<2e-9))
        frames.append(frame)
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    passed=all(f['interfacePassed'] for f in frames)
    result=dict(interfacePassed=passed,frames=frames,inputSHA256=locks,inputsStable=stable,
        exactStage='laminar-same-flow-state-no-SST',givenWallTemperatureIsSolved=False,
        experimentalQualification=False,newQualifiedArea_m2=0,threeUpdatesAreNotConvergence=True,
        elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',result);print(dict(interfacePassed=passed,frames=len(frames)))
    assert passed,'原严格接口门槛失败，保留原报告，不重算CFD'


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();audit(Path(args.folder).resolve(),Path(args.output).resolve())
