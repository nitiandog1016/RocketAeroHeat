# -*- coding: utf-8 -*-
"""独立重建全部合法候选和旧特征分支，并按公开曲线复算被拒面的Riemann候选。"""
import json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT,VAL,sha,exact,PUBLIC_BINDINGS

case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_independent_boundary_path.json');assert not target.exists()
file=VAL/case/'newton_boundary_path_evidence.mat';s=loadmat(file,simplify_cells=True);gas=s['gas'];mapping=s['map']
bindings={str(file.relative_to(ROOT)):sha(file),**PUBLIC_BINDINGS}
source=s['report']['sourceFile'];assert sha(source)==s['report']['sourceSHA256'];bindings[str(__import__('pathlib').Path(source).relative_to(ROOT))]=sha(source)
owner=np.asarray(mapping['ownerCellID'],int)-1;n=np.asarray(mapping['normal']);res=s['reservoir'];g=gas['gamma'];R=gas['R'];cv=R/(g-1)
rows=lambda x:[x] if isinstance(x,dict) else list(x)
details=[];worst_flux=0.;worst_state=0.;exact_count=0
for trial in rows(s['trials']):
    q=s['Q']+trial['alpha']*s['direction']*s['scale'];rho=q[:,0];u=q[:,1:4]/rho[:,None]
    T=(q[:,4]-.5*np.sum(q[:,1:4]**2,axis=1)/rho)/(rho*cv)
    p=np.column_stack((rho,u,T))[owner];np.testing.assert_allclose(p,trial['primitive'],rtol=4e-15,atol=1e-13)
    un=np.sum(p[:,1:4]*n,axis=1);a=np.sqrt(g*R*p[:,4]);ur=np.sum(res[:,1:4]*n,axis=1);ar=np.sqrt(g*R*res[:,4])
    ub=.5*(un+2*a/(g-1)+ur-2*ar/(g-1));ab=(g-1)*.25*(un+2*a/(g-1)-ur+2*ar/(g-1))
    sub=abs(un)<a
    failed=(un < -a)&(ur >= -ar)
    failed|=sub&((ab<=0)|(abs(ub)>=ab-128*np.spacing(np.maximum.reduce((np.ones_like(ab),ab,abs(ub))))))
    np.testing.assert_array_equal(failed,np.asarray(trial['failed'],bool))
    np.testing.assert_allclose(trial['cellNormalMach'],un/a,rtol=4e-15,atol=1e-14)
    np.testing.assert_allclose(trial['acousticSubsonicMarginMs'],ab-abs(ub),rtol=2e-11,atol=3e-12)
    ids=np.flatnonzero(failed);assert not trial['acceptedByOriginalSolver'] or len(ids)==0
    if len(ids):
        assert int(trial['firstFailedIndex'])==int(ids[0]+1)
        e=trial['exactReservoirDiagnostic'];assert e['succeeded']
        fluxes=np.atleast_2d(e['flux']);states=np.atleast_2d(e['samplePrimitive']);regions=np.atleast_1d(e['sampleRegion'])
        for j,k in enumerate(ids):
            state,flux,region=exact(p[k],res[k],n[k],gas)
            worst_flux=max(worst_flux,float(np.max(abs(flux-fluxes[j]))/max(1,np.max(abs(flux)))))
            worst_state=max(worst_state,float(np.max(abs(state-states[j]))/max(1,np.max(abs(state)))))
            assert region==regions[j] and 80<=state[4]<=2000
            exact_count+=1
    details.append(dict(alpha=float(trial['alpha']),failedBoundaryFaces=int(len(ids)),
                        firstOpenFaceID=int(mapping['openFaceID'][ids[0]]) if len(ids) else None))
assert exact_count==18 and worst_flux<3e-12 and worst_state<3e-12
assert all(sha(ROOT/p)==digest for p,digest in bindings.items())
out=dict(status='passed-independent-full-legal-boundary-path-and-18-public-Riemann-reference-faces',processID=os.getpid(),
    trialDetails=details,exactReferenceFaces=exact_count,maximumFluxRelativeDifference=worst_flux,maximumPrimitiveRelativeDifference=worst_state,
    oldCharacteristicGuardPreserved=True,fullNSRecomputed=False,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__),referenceHelperSHA256=sha(__import__('riemannReservoirReference').__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
