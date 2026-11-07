# -*- coding: utf-8 -*-
"""原生Euler有限控制收尾：完成、发展、参照和资格分层，输出科研对照图。"""
from pathlib import Path
import argparse
import csv
import json
import re
import shutil
import time
import numpy as np
from readHybridEulerNative import read_wall_frame
from readEilmerNative import read_grid,read_field,sha
from eilmerAirThermo import evaluate


def audit(case, preparation, run_receipt, reference, out):
    start=time.monotonic();case,run_receipt,out=map(Path,(case,run_receipt,out))
    out.mkdir(parents=True,exist_ok=False)
    files=[Path(__file__)]+[Path(__file__).with_name(n) for n in
        ['readHybridEulerNative.py','readEilmerNative.py','eilmerAirThermo.py','hybridExternalState.py']]
    for p in files:shutil.copyfile(p,out/p.name)
    inputs=files+[Path(preparation),run_receipt,Path(reference)]
    before={str(p):sha(p) for p in inputs};result=dict(passed=False,physicalQualification=False)
    try:
        record=json.loads(run_receipt.read_text());process=json.loads((run_receipt.parent/'native/report.json').read_text())
        assert record['processPassed'] and process['passed'] and not process['liveGroupAfterCleanup']
        text=(run_receipt.parent/'native/stdout.txt').read_text()
        assert not re.search('OutOfMemoryError|integrateInTime failed|An exception was caught|Begin exception message',text)
        assert re.findall(r'^FINAL-STEP: (\d+)',text,re.M)==['400']
        assert re.findall(r'^STOP-REASON: (.+)',text,re.M)==['maximum-steps']
        assert re.findall(r'^Number-type: (.+)',text,re.M)==['real']
        assert list(map(int,re.findall(r'^Step=\s*(\d+)',text,re.M)))==list(range(1,401))
        storage=re.findall(r'^RAH-STORAGE cells=(\d+) ghosts=(\d+) neq=(\d+) nftl=(\d+) scalarBytes=(\d+) stageLength=(\d+) payloadBytes=(\d+)',text,re.M)
        assert storage==[('2400','2720','5','4','8','25600','1638400')]
        assert 'RAH-WORKSPACE jacobian=0 sourceBackup=0 viscousBackup=0 interpBackup=0 perturbBackup=0 viscous=5120 interp=5120 interpWorkspace=5120' in text
        metadata=(case/'lmrsim/snapshots/snapshot-times-metadata').read_text()
        entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)",metadata,re.M)
        assert [(int(i),int(s)) for i,s,t in entries]==[(0,0),(1,200),(2,400)]
        refs=json.loads(Path(reference).read_text());h0=refs['states'][0]['totalEnthalpyJkg']
        frames=[];volume_rows=[];grid=read_grid(case/'lmrsim/grid/grid-0000')
        for snapshot,step,t in entries:
            snapshot=int(snapshot);f,info=read_wall_frame(case,snapshot,preparation,run_receipt)
            frames.append(f);(out/('wall_%04d.json'%snapshot)).write_text(json.dumps(f,indent=2)+'\n')
            raw=read_field(case,snapshot,2400);gas=evaluate(raw['p'],raw['T'])
            assert all(np.all(raw[k]>0) for k in ['p','rho','T','vol'])
            v=np.column_stack([raw['vel.'+k] for k in 'xyz']);hh=gas['enthalpyJkg']+np.sum(v*v,axis=1)/2
            boundary={}
            for name,b in grid['boundaries'].items():
                idx=b['owners'];speed=np.linalg.norm(v[idx],axis=1)
                boundary[name]=dict(minMach=float(np.min(speed/gas['soundSpeedMps'][idx])),
                    pRangePa=[float(raw['p'][idx].min()),float(raw['p'][idx].max())],
                    maxAbsVyMps=float(abs(v[idx,1]).max()))
            volume_rows.append(dict(snapshot=snapshot,step=int(step),temperatureRangeK=[float(raw['T'].min()),float(raw['T'].max())],
                maximumTotalEnthalpyRelativeError=float(np.max(abs(hh/h0-1))),
                volumeMeanTotalEnthalpyRelativeError=float(np.average(abs(hh/h0-1),weights=raw['vol'])),
                fullVolumeEOSErrors=f['fullVolumeEOSErrors'],boundaries=boundary))
        last=frames[-1];prev=frames[-2];areas=np.array(last['faceAreasM2']);centers=np.array(last['faceCentersM'])
        s=-centers[:,0];L1={}
        for key in ['pressurePa','temperatureK','velocityMps']:
            x=np.array(last['fields'][key]);old=np.array(prev['fields'][key])
            if x.ndim==2:d=np.linalg.norm(x-old,axis=1);scale=np.linalg.norm(x,axis=1)
            else:d=abs(x-old);scale=abs(x)
            L1[key]=float(np.sum(areas*d)/np.sum(areas*scale))
        comparisons=[]
        for name,lo,hi,ref in [('compression',.35,.5,refs['states'][1]),('expansion',.85,.95,refs['states'][2])]:
            mask=(s>=lo)&(s<=hi);w=areas[mask];g=last['fields']
            p=float(np.average(np.array(g['pressurePa'])[mask],weights=w));T=float(np.average(np.array(g['temperatureK'])[mask],weights=w))
            v=np.average(np.array(g['velocityMps'])[mask],weights=w,axis=0);ref_v=np.array(ref['velocityMps'])
            angle=float(np.degrees(np.arccos(np.clip(np.dot(v,ref_v)/np.linalg.norm(v)/np.linalg.norm(ref_v),-1,1))))
            comparisons.append(dict(region=name,streamwiseIntervalM=[lo,hi],faceCount=int(mask.sum()),
                 pressurePa=p,referencePressurePa=ref['pressurePa'],pressureRelativeError=p/ref['pressurePa']-1,
                 temperatureK=T,referenceTemperatureK=ref['temperatureK'],temperatureRelativeError=T/ref['temperatureK']-1,
                 velocityMps=v.tolist(),referenceVelocityMps=ref_v.tolist(),directionErrorDeg=angle))
        history=np.loadtxt(case/'lmrsim/transient-residuals.dat');assert history.shape==(400,10)
        assert np.isfinite(history).all() and np.array_equal(history[:,0],np.arange(1,401))
        improvement={name:float(np.mean(np.log10(history[:50,col]))-np.mean(np.log10(history[-50:,col])))
                     for name,col in [('mass',3),('energy',7)]}
        # 原生L2有量纲；同时给原值和首步归一化值，禁止与其他引擎−8直接等同。
        nativeL2=float(history[-1,8]);normalizedL2=float(history[-1,8]/history[0,8])
        comparison_pass=all(abs(d['pressureRelativeError'])<=.05 and abs(d['temperatureRelativeError'])<=.05
                            and d['directionErrorDeg']<=2 for d in comparisons)
        result.update(passed=True,completedUpdates=400,savedSteps=[0,200,400],storageRows=storage,
            nativeSeconds=process['elapsedSeconds'],peakGroupRSSBytes=process['peakObservedGroupRSSBytes'],
            processInputsStable=process['inputsStable'],physicalQualification=False,newSolvedHeatAreaM2=0,
            frames=volume_rows,wallAreaM2=last['originalAreaM2'],wallChanges200To400=L1,
            platformComparisons=comparisons,coarseReferenceScreenPassed=comparison_pass,
            residualLogImprovement=improvement,developmentScreenPassed=all(v>0 for v in improvement.values()),
            finalNativeL2=nativeL2,finalL2RelativeToFirst=normalizedL2,
            convergencePassed=bool(nativeL2<=1e-8 and normalizedL2<=1e-8 and max(L1.values())<=.01),
            convergenceScope='RUN_PLAN global residual did not define normalization; neither native nor first-step normalized L2 passes; no qualification',
            qAvailable=False,edgeQualified=False)
        with (out/'wall_comparison.csv').open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['faceID','s_m','y_m','area_m2','p_Pa','T_K','Ux_mps','Uy_mps','Uz_mps','h0_Jkg','q_status'])
            for i in range(240):writer.writerow([i+1,s[i],centers[i,1],areas[i],last['fields']['pressurePa'][i],
                last['fields']['temperatureK'][i]]+last['fields']['velocityMps'][i]+[last['fields']['totalEnthalpyJkg'][i],'not-solved-by-Euler'])
        # 模块三：独立科研图，不改主UI；场值与坐标/参考身份一同呈现。
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(10,6.5),constrained_layout=True)
        for ax,key,label in [(axes[0,0],'pressurePa','Wall-adjacent pressure [Pa]'),(axes[0,1],'temperatureK','Temperature [K]')]:
            for frame,step in [(frames[1],200),(frames[2],400)]:
                ax.scatter(s,frame['fields'][key],s=6,label='step %d'%step)
            for lo,hi,ref in [(.35,.5,refs['states'][1]),(.85,.95,refs['states'][2])]:
                ax.plot([lo,hi],[ref[key]]*2,'k--',lw=1.4)
            ax.set(xlabel='s = -X [m]',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=8)
        axes[1,0].plot(s,np.degrees(np.arctan2(np.array(last['fields']['velocityMps'])[:,2],-np.array(last['fields']['velocityMps'])[:,0])),'.',ms=3)
        axes[1,0].set(xlabel='s = -X [m]',ylabel='Velocity turn [deg]');axes[1,0].grid(alpha=.2)
        for name,col in [('mass',3),('energy',7),('native L2',8)]:
            axes[1,1].semilogy(history[:,0],history[:,col]/history[0,col],label=name)
        axes[1,1].axhline(1e-8,color='k',ls='--',lw=.8)
        axes[1,1].set(xlabel='Local pseudo-time updates',ylabel='Residual / first update');axes[1,1].legend(fontsize=8);axes[1,1].grid(alpha=.2)
        fig.suptitle('M3 compression / expansion: CEA air Euler control\nUnconverged research states; q and solid wall T not solved',fontsize=11)
        fig.savefig(out/'control_comparison.png',dpi=160);plt.close(fig)
    except Exception as error:result['error']=repr(error)
    result['elapsedSeconds']=time.monotonic()-start
    result['inputSHA256']=before;result['inputsStable']=all(sha(p)==d for p,d in before.items())
    (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','developmentScreenPassed','coarseReferenceScreenPassed','convergencePassed','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','preparation','run-receipt','reference','out']:p.add_argument('--'+k,required=True)
    a=p.parse_args();r=audit(a.case,a.preparation,a.run_receipt,a.reference,a.out)
    raise SystemExit(0 if r['passed'] and r['inputsStable'] else 1)
