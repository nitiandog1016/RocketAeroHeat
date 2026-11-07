# -*- coding: utf-8 -*-
"""将Euler迭代收敛与网格/理论差别分开；不把数值收敛授为物理/热资格。"""
from pathlib import Path
import argparse
import json
import re
import shutil
import time
import numpy as np
from readHybridEulerNative import read_wall_frame
from readEilmerNative import sha,read_field
from eilmerAirThermo import evaluate


def audit(case, preparation, run_receipt, reference, out):
    case,prep,rp,reference,out=map(Path,(case,preparation,run_receipt,reference,out))
    out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    src=[Path(__file__)]+[Path(__file__).with_name(n) for n in ['hybridEulerEvolution.py',
         'readHybridEulerNative.py','readEilmerNative.py','eilmerAirThermo.py','hybridExternalState.py']]
    for p in src:shutil.copyfile(p,out/p.name)
    inputs=src+[prep,rp,reference];hashes={str(p):sha(p) for p in inputs}
    result=dict(passed=False,physicalQualification=False,newSolvedHeatAreaM2=0)
    try:
        r=json.loads(prep.read_text());run=json.loads(rp.read_text());process=json.loads((rp.parent/'native/report.json').read_text())
        assert run['processPassed'] and process['passed'] and not process['liveGroupAfterCleanup']
        text=(rp.parent/'native/stdout.txt').read_text();final=r['expectedControl']['max_step'];first=r.get('sourceStep',0)
        assert re.findall(r'^FINAL-STEP: (\d+)',text,re.M)==[str(final)]
        assert re.findall(r'^STOP-REASON: (.+)',text,re.M)==['maximum-steps']
        assert re.findall(r'^Number-type: (.+)',text,re.M)==['real']
        assert list(map(int,re.findall(r'^Step=\s*(\d+)',text,re.M)))==list(range(first+1,final+1))
        if first:
            assert re.findall(r'^RESTART-STEP: (\d+)',text,re.M)==[str(first)]
            assert re.findall(r'^RESTART-SNAPSHOT: (\d+)',text,re.M)==[str(r['sourceSnapshot'])]
            from hybridEulerEvolution import local_path
            origin=local_path(r['baselineConfig']).parents[1]
            prefix=(origin/'lmrsim/transient-residuals.dat').read_bytes()
            assert (case/'lmrsim/transient-residuals.dat').read_bytes().startswith(prefix)
        else:assert not re.findall(r'^RESTART-STEP:',text,re.M)
        history=np.loadtxt(case/'lmrsim/transient-residuals.dat');assert history.shape==(final,10)
        assert np.isfinite(history).all() and np.array_equal(history[:,0],np.arange(1,final+1))
        ref=json.loads(reference.read_text());free=ref['states'][0]
        rho=free['densityKgm3'];U=np.linalg.norm(free['velocityMps']);L=1.
        scale=np.array([rho*U/L]+[rho*U**2/L]*3+[rho*U**3/L]);normal=history[:,3:8]/scale
        maxima=normal[-50:].max(axis=0);massratio=float(history[-1,9]/(rho*U*.024))
        meta=(case/'lmrsim/snapshots/snapshot-times-metadata').read_text()
        entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)",meta,re.M)
        expected=([0,200,400]+list(range(800,final+1,400))) if first else list(range(0,final+1,400))
        assert [int(e[1]) for e in entries]==expected
        frames=[];allstates=[]
        for i,step,t in entries:
            frame,info=read_wall_frame(case,int(i),prep,rp);frames.append(frame)
            f=read_field(case,int(i),info['cells']);gas=evaluate(f['p'],f['T'])
            assert all(np.all(f[k]>0) for k in ['p','rho','T','vol'])
            assert np.all((f['T']>=200)&(f['T']<800))
            v=np.column_stack([f['vel.'+a] for a in 'xyz']);h0=gas['enthalpyJkg']+np.sum(v*v,axis=1)/2
            allstates.append(dict(snapshot=int(i),step=int(step),temperatureRangeK=[float(f['T'].min()),float(f['T'].max())],
                fullVolumeEOSErrors=frame['fullVolumeEOSErrors'],maximumTotalEnthalpyRelativeError=float(np.max(abs(h0/free['totalEnthalpyJkg']-1)))))
        last=frames[-1];prev=frames[-2];areas=np.array(last['faceAreasM2']);s=-np.array(last['faceCentersM'])[:,0]
        changes={}
        for key in ['pressurePa','temperatureK','velocityMps']:
            a=np.array(last['fields'][key]);b=np.array(prev['fields'][key])
            err=abs(a-b) if a.ndim==1 else np.linalg.norm(a-b,axis=1)
            den=abs(a) if a.ndim==1 else np.linalg.norm(a,axis=1)
            changes[key]=float(np.sum(err*areas)/np.sum(den*areas))
        comparisons=[]
        for name,lo,hi,refstate in [('compression',.35,.5,ref['states'][1]),('expansion',.85,.95,ref['states'][2])]:
            mask=(s>=lo)&(s<=hi);values={key:np.average(np.array(last['fields'][key])[mask],weights=areas[mask],axis=0)
                                      for key in ['pressurePa','temperatureK','velocityMps']}
            v=values['velocityMps'];vr=np.array(refstate['velocityMps'])
            comparisons.append(dict(region=name,faceCount=int(mask.sum()),pressurePa=float(values['pressurePa']),temperatureK=float(values['temperatureK']),
                velocityMps=v.tolist(),pressureRelativeError=float(values['pressurePa']/refstate['pressurePa']-1),
                temperatureRelativeError=float(values['temperatureK']/refstate['temperatureK']-1),
                directionErrorDeg=float(np.degrees(np.arccos(np.clip(np.dot(v,vr)/np.linalg.norm(v)/np.linalg.norm(vr),-1,1))))))
        storage=re.findall(r'^RAH-STORAGE cells=(\d+) ghosts=(\d+) neq=(\d+) nftl=(\d+) scalarBytes=(\d+) stageLength=(\d+) payloadBytes=(\d+)',text,re.M)
        nc=info['cells'];ng=info['ghostCells'];expected_storage=tuple(map(str,[nc,ng,5,4,8,(nc+ng)*5,(nc+ng)*5*4*8*2]))
        assert storage==[expected_storage]
        result.update(passed=True,profile=r['profile'],sourceStep=first,finalStep=final,newUpdates=final-first,
            savedSteps=expected,frames=allstates,normalization=dict(lengthM=L,densityKgm3=rho,speedMps=float(U),
            equationOrder=['mass','x-momentum','y-momentum','z-momentum','total-energy'],scales=scale.tolist(),
            last50MaximumNormalizedResiduals=maxima.tolist(),limit=1e-8),
            residualConverged=bool(np.all(maxima<=1e-8)),nativeL2MixedUnits=float(history[-1,8]),
            massFluxImbalanceRelativeToFreestream=massratio,wallLastPairAreaL1=changes,
            numericalConvergencePassed=bool(np.all(maxima<=1e-8) and max(changes.values())<=.01 and massratio<=1e-6),
            platformComparisons=comparisons,coarseReferenceScreenPassed=all(abs(c['pressureRelativeError'])<=.05
                    and abs(c['temperatureRelativeError'])<=.05 and c['directionErrorDeg']<=2 for c in comparisons),
            nativeSeconds=process['elapsedSeconds'],peakGroupRSSBytes=process['peakObservedGroupRSSBytes'],
            storage=storage,wallAreaM2=last['originalAreaM2'],physicalQualification=False,edgeQualified=False)
        (out/'final_wall.json').write_text(json.dumps(last,indent=2)+'\n')
        np.savetxt(out/'normalized_residuals.csv',np.column_stack([history[:,0],normal]),delimiter=',',
                   header='step,mass,xmom,ymom,zmom,energy',comments='',fmt='%.17g')
    except Exception as error:result['error']=repr(error)
    result.update(inputSHA256=hashes,inputsStable=all(sha(p)==d for p,d in hashes.items()),elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','numericalConvergencePassed','platformComparisons','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','preparation','run-receipt','reference','out']:p.add_argument('--'+k,required=True)
    a=p.parse_args();r=audit(a.case,a.preparation,a.run_receipt,a.reference,a.out)
    raise SystemExit(0 if r['passed'] else 1)
