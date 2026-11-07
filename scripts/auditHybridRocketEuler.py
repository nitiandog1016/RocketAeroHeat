# -*- coding: utf-8 -*-
"""真实Euler失败/发展独立收尾；异常场只作原始诊断，不绕过生产导出拒绝。"""
from pathlib import Path
import argparse
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import sha,read_grid,read_field,compare_source
from hybridRocketEuler import local,check_inputs,read_wall,validate_terminal,save
from eilmerAirThermo import evaluate
from hybridEdgeSampling import map_original_wall


def audit(case, preparation, run_receipt, out):
    start=time.monotonic();case,preparation,rp,out=map(Path,(case,preparation,run_receipt,out))
    out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_audit.py')
    r=check_inputs(case,preparation);run=json.loads(rp.read_text());process=json.loads((rp.parent/'native/report.json').read_text())
    assert not process['liveGroupAfterCleanup'] and process['inputsStable']
    for p,d in run['outputSHA256'].items():assert sha(local(p))==d
    text=(rp.parent/'native/stdout.txt').read_text();g=r['geometry'];n=g['cells'];ng=g['ghostCells'];total=n+ng
    storage=re.findall(r'^RAH-STORAGE cells=(\d+) ghosts=(\d+) neq=(\d+) nftl=(\d+) scalarBytes=(\d+) stageLength=(\d+) payloadBytes=(\d+)',text,re.M)
    assert storage==[tuple(map(str,[n,ng,5,4,8,total*5,total*5*4*8*2]))]
    assert ('RAH-WORKSPACE jacobian=0 sourceBackup=0 viscousBackup=0 interpBackup=0 perturbBackup=0 '
            'viscous=%d interp=%d interpWorkspace=%d'%(total,total,total)) in text
    assert 'RAH-WORKSPACE-CELLS internal=64982 implicitObjects=0' in text
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,_=compare_source(grid,case/'mesh.su2')
    surface=json.loads((case/'surface.json').read_text());V=np.array(surface['vertices']);F=np.array(surface['facesOneBased'])-1;tri=V[F]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert abs(area.sum()-1.9646292565715551)<2e-13
    meta=(case/'lmrsim/snapshots/snapshot-times-metadata').read_text()
    entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)\s*\n\s+dt:\s*(\S+)",meta,re.M)
    fields=[];records=[]
    for index,step,t,dt in entries:
        raw=read_field(case,int(index),n);gas=evaluate(raw['p'],raw['T'])
        loc=map_original_wall(surface,grid,mapping['wall'],np.column_stack([raw['pos.'+a] for a in 'xyz']))
        owners=np.array(loc['nativeOwnerCellsZeroBased']);U=np.column_stack([raw['vel.'+k] for k in 'xyz'])
        errors={key:float(np.max(abs(raw[key]-gas[name])/np.maximum(1,abs(gas[name]))))
                for key,name in [('rho','densityKgm3'),('e','internalEnergyJkg'),('a','soundSpeedMps')]}
        assert max(errors.values())<2e-10 and np.all(raw['rho']>0) and np.all(raw['p']>0)
        h0=gas['enthalpyJkg']+np.sum(U*U,axis=1)/2
        records.append(dict(snapshot=int(index),step=int(step),pseudoTime=float(t),reportedDt=float(dt),
            temperatureRangeK=[float(raw['T'].min()),float(raw['T'].max())],pressureRangePa=[float(raw['p'].min()),float(raw['p'].max())],
            fullVolumeEOSErrors=errors,totalEnthalpyRangeJkg=[float(h0.min()),float(h0.max())],
            localDtRangeS=[float(raw['dt_local'].min()),float(raw['dt_local'].max())] if int(index)>0 else None))
        fields.append((raw,owners))
    assert len(fields)>=1
    history=np.atleast_2d(np.loadtxt(case/'lmrsim/transient-residuals.dat'));assert np.isfinite(history).all() and history.shape[1]==10
    free=json.loads((case/'native-input-state.json').read_text());speed=np.linalg.norm([free[k] for k in ['velx','vely','velz']])
    L=float(np.ptp(V[:,0]));rho=free['rho'];scale=np.array([rho*speed/L]+[rho*speed**2/L]*3+[rho*speed**3/L])
    normalized=history[:,3:8]/scale
    np.savetxt(out/'normalized_residuals.csv',np.c_[history[:,0],normalized],delimiter=',',comments='',
               header='step,mass,x_momentum,y_momentum,z_momentum,energy',fmt='%.17g')
    # 独立给出原几何的原生公式长度和有限时间步界；并非运行时对象的实测长度。
    pts=grid['points'][grid['cells'][:,:4]];center=pts.mean(1)
    volume=np.linalg.det((pts[:,1:]-pts[:,:1]).transpose(0,2,1))/6
    lengths=np.cbrt(volume)
    for face in [[0,1,2],[0,1,3],[1,2,3],[2,0,3]]:
        lengths=np.minimum(lengths,2*np.linalg.norm(pts[:,face].mean(1)-center,axis=1))
    assert np.all(lengths>0)
    dt_bounds=[float(.3*lengths.min()/(speed+free['a'])),float(.3*lengths.min()/free['a'])]
    # 生产入口必须拒绝本轮异常任务；初场与异常诊断分开。
    rejected=False
    try:validate_terminal(case,preparation,rp)
    except ValueError:rejected=True
    stop_path=rp.parent/'invalid_time_stop.json'
    stop=json.loads(stop_path.read_text()) if stop_path.exists() else None
    actual_steps=[e['step'] for e in records]
    planned=bool(process['passed'] and actual_steps[0]==0 and 200 in actual_steps and actual_steps[-1]==400)
    reason='normal-completed' if planned else 'incomplete-or-failed-development'
    if stop:
        assert process['nativeExitCode']==-15 and stop['notUserCancellation'] and rejected
        assert len(history)==1 and actual_steps==[0,1] and records[1]['reportedDt']>1e90
        reason='manual-scientific-stop-after-first-invalid-pseudo-time; not temperature/resource/user stop'
    raw,owners=fields[-1]
    np.savetxt(out/'rejected_original_wall_diagnostic.csv',np.c_[np.arange(1,19753),owners,area,tri.mean(1),
        raw['p'][owners],raw['T'][owners],raw['vel.x'][owners],raw['vel.y'][owners],raw['vel.z'][owners]],
        delimiter=',',comments='',header='face_id,owner_zero_based,area_m2,x_m,y_m,z_m,p_raw_Pa,T_raw_fluid_K,Ux_mps,Uy_mps,Uz_mps',fmt='%.17g')
    result=dict(passed=True,scope='read-only audit of terminal task; diagnostic success is not candidate acceptance',
        candidateAccepted=planned,productionExportRejected=rejected,reason=reason,frames=records,
        savedSteps=actual_steps,historyRows=len(history),originalFaces=19752,originalAreaM2=float(area.sum()),
        normalization=dict(lengthM=L,densityKgm3=rho,speedMps=float(speed),scales=scale.tolist()),
        finalNormalizedComponents=normalized[-1].tolist(),numericalConvergencePassed=False,physicalQualification=False,
        rawSourceLengthRangeM=[float(lengths.min()),float(lengths.max())],initialFiniteDtBoundsFromGeometry=dt_bounds,
        runtimeCellLengthsObserved=False,uniqueTimestepCauseEstablished=False,
        processNativeExitCode=process['nativeExitCode'],nativeSupervisedSeconds=process['elapsedSeconds'],
        peakGroupRSSBytes=process['peakObservedGroupRSSBytes'],inputCount=len(process['inputSHA256']),inputsStable=True,
        newSolvedHeatAreaM2=0,sourceSHA256={str(p.resolve()):sha(p) for p in [Path(__file__),rp,preparation,rp.parent/'native/report.json']},
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    print(json.dumps({k:result[k] for k in ['passed','candidateAccepted','reason','peakGroupRSSBytes','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','preparation','run-receipt','out']:p.add_argument('--'+k,required=True)
    a=p.parse_args();audit(a.case,a.preparation,a.run_receipt,a.out)
