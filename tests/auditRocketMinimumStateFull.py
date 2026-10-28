"""逐个真实整箭检查点核对原NS方向、连续状态、残差、温域及共享面守恒。"""
from pathlib import Path
import json, hashlib, math, sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
def rows(x):return [x] if isinstance(x,dict) else list(x)
def array_sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def norm(a):return np.linalg.norm(np.asarray(a).ravel())
def audit(folder):
    report=json.loads((folder/'full_comparison.json').read_text('utf8'))
    d=loadmat(ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat',simplify_cells=True,variable_names=['domain'])['domain']
    f=d['openFaces'];w=d['walls'];vol=d['cells']['volume'];n=len(vol)
    first=loadmat(ROOT/'outputs/v83_validation/krylov24/wall_velocity_krylov_evidence.mat',simplify_cells=True,variable_names=['result'])['result']
    start_sha=array_sha(first['Q']);start_rate_sha=array_sha(first['rate']);del first
    results=[]
    for branch in report['branches']:
        previous_sha=start_sha;previous_rate_sha=start_rate_sha;history=[];main_calls=0;probe_calls=0
        for step in branch['steps']:
            if step['probe']:probe_calls+=step['probe']['newActualNSEvaluations']
            if not step['main']:continue
            file=ROOT/step['resultFolder']/'wall_velocity_krylov_evidence.mat'
            s=loadmat(file,simplify_cells=True,variable_names=['result','gas','boundary']);r=s['result'];g=s['gas'];scale=r['residualScale']
            assert array_sha(r['initialQ'])==previous_sha and array_sha(r['initialRate'])==previous_rate_sha
            previous_sha=array_sha(r['Q']);previous_rate_sha=array_sha(r['rate'])
            Q=r['Q'];rho=Q[:,0];u=Q[:,1:4]/rho[:,None];T=(Q[:,4]/rho-.5*np.sum(u*u,axis=1))/(g['R']/(g['gamma']-1))
            assert np.all(rho>0) and np.all(np.isfinite(Q)) and np.min(T)>=80 and np.max(T)<=2000
            assert g['temperatureRange'].tolist()==[80,2000] and s['boundary']['solidWall']['temperature']==300
            logs=rows(r['linearLog']);trace=rows(r['trace']);trials=rows(r['lineLog'])
            assert len(trace)==r['acceptedNonlinearIterations']==1
            e=logs[0]['steadyDirectionEvidence'];jv=(e['ratePlus']-e['rateMinusOrBase'])/scale/e['actualScalarStep'];base=r['initialRate']/scale
            eta=norm(jv+base)/norm(base);slope=np.sum(jv*base)/np.sum(base*base)
            assert eta<=.1 and abs(eta-logs[0]['fullSystemRelativeDefect'])<1e-10 and slope<=-1+eta+1e-10
            t=trace[0];candidate=r['initialQ']+t['alpha']*t['linearUpdateU']*scale
            np.testing.assert_array_equal(candidate,Q)
            selected=loadmat(file.parent.parent/f"step{step['iteration']}"/'selected_direction.mat',simplify_cells=True,variable_names=['direction'])['direction']
            np.testing.assert_array_equal(selected,t['linearUpdateU'])
            scaled=r['rate']/scale;rms=float(np.sqrt(np.mean(scaled*scaled)));maximum=float(np.max(np.abs(scaled)))
            assert abs(rms-r['rmsMerit'])/rms<1e-12 and maximum==r['maximumResidual']
            assert rms<=(1-1e-4*t['alpha'])*r['initialRmsMerit']
            # 全部共享面五守恒量，以独立bincount装配；不复算一个相似的简化残差。
            info=r['info'];integrated=np.zeros_like(Q)
            for c in range(5):
                value=f['area']*info['faceFlux'][:,c]
                mask=f['cellMinus']>0;integrated[:,c]-=np.bincount(f['cellMinus'][mask].astype(int)-1,weights=value[mask],minlength=n)
                mask=f['cellPlus']>0;integrated[:,c]+=np.bincount(f['cellPlus'][mask].astype(int)-1,weights=value[mask],minlength=n)
                mask=w['ownerCell']>0;integrated[:,c]-=np.bincount(w['ownerCell'][mask].astype(int)-1,weights=w['area'][mask]*info['wallFlux'][mask,c],minlength=n)
            divergence_error=norm((integrated/vol[:,None]-r['rate'])/scale)/norm(scaled);assert divergence_error<1e-12
            mask=f['cellPlus']==0
            denominator=np.maximum(1,np.sum(np.abs(info['faceFlux'][mask])*f['area'][mask,None],axis=0)+np.sum(np.abs(info['wallFlux'])*w['area'][:,None],axis=0))
            spatial=float(np.max(np.abs(np.sum(integrated,axis=0)+info['netBoundaryFlux'])/denominator))
            steady=float(np.max(np.abs(info['netBoundaryFlux'])/denominator));assert spatial<=1e-9
            assert abs(steady-r['steadyBoundaryFluxBalance'])<1e-12 and r['options']['steadyResidualTolerance']==2e-6 and r['options']['steadyFluxTolerance']==5e-10
            # 所有候选的单元护栏及实际原残差工作量。
            candidate_calls=0
            for line in trials:
                q=r['initialQ']+line['alpha']*t['linearUpdateU']*scale;rh=q[:,0]
                if np.any(rh<=0):continue
                vel=q[:,1:4]/rh[:,None];ei=q[:,4]-.5*rh*np.sum(vel*vel,axis=1)
                if np.any(ei<=0):continue
                temp=ei/(rh*g['R']/(g['gamma']-1))
                if np.any(temp<80)|np.any(temp>2000):continue
                candidate_calls+=1
            jvp_calls=sum(x['actualOriginalResidualEvaluations'] for x in rows(r['jvpLog']))
            assert r['residualEvaluations']==1+jvp_calls+candidate_calls
            main_calls+=r['residualEvaluations']
            history.append({'iteration':step['iteration'],'alpha':t['alpha'],'RMS':rms,'maximum':maximum,
                'internalC':step['internalCRelativeResidual'],'originalNSDirection':eta,
                'minimumDensity':float(np.min(rho)),'temperatureRange':[float(np.min(T)),float(np.max(T))],
                'spatialBalance':spatial,'steadyBoundaryImbalance':steady,'divergenceRelativeError':divergence_error,
                'mainNSCalls':int(r['residualEvaluations']),'stepTotalSeconds':step['elapsedSeconds'],
                'stateSHA256':previous_sha,'steadyConverged':bool(r['steadyConverged'])})
            del s,r,info,selected,logs,trace,trials,e
        results.append({'method':branch['method'],'history':history,'mainNSCalls':main_calls,'probeNSCalls':probe_calls,
            'branchElapsedSeconds':branch['elapsedSeconds'],'acceptedSteps':len(history)})
    out={'status':'passed-independent-complete-domain-continuous-state-direction-gates-and-all-face-divergence',
        'source':str(folder.relative_to(ROOT)),'branches':results,'initialStateSHA256':start_sha,
        'physicalSolution':False,'newProductionSolvedAreaM2':0}
    (folder/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps(out))
if __name__=='__main__':audit(OUT/(sys.argv[1] if len(sys.argv)>1 else 'full02'))
