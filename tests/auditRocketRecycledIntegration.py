# -*- coding: utf-8 -*-
"""独立核对增强主链的模型/预条件预算、周期正交/当前像及总方向身份。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];kind=sys.argv[2]
assert case.replace('_','').isalnum() and kind in ('unit','real');target=VAL/(case+'_independent_recycled_integration.json');assert not target.exists()
file=VAL/case/('steady_newton_unit_evidence.mat' if kind=='unit' else 'wall_velocity_krylov_evidence.mat');s=loadmat(file,simplify_cells=True)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

rows=lambda x:[x] if isinstance(x,dict) else list(np.atleast_1d(x))
bindings={str(file.relative_to(ROOT)):sha(file)};checks=[]
if kind=='unit':
    kernelFile=VAL/'recycledkernel01/recycled_kernel_evidence.mat';kernel=loadmat(kernelFile,simplify_cells=True,variable_names=['unitA','Q'])
    A=kernel['unitA'];bindings[str(kernelFile.relative_to(ROOT))]=sha(kernelFile)
for variant in (['result','limited'] if kind=='unit' else ['result']):
    r=s[variant];assert r['krylovMethod']=='SVD-selected-augmented-restarted-right-GMRES'
    total_steps=0;total_actions=0;total_preconditions=0;maximum_retained=0
    for a in rows(r['linearLog']):
        if 'recycledKrylovReport' not in a:continue
        rr=a['recycledKrylovReport'];cycles=rows(rr['cycles']);expected=0;steps=0
        assert rr['basisStorageUpperBoundBytes']<=rr['maximumBasisBytes'] and not rr['originalNSAcceptanceClaimed']
        assert rr['maximumBasisOrthogonalityError']<1e-8 and rr['maximumRetainedImageProbeRelativeDifference']<1e-8
        for c in cycles:
            j=int(c['arnoldiSteps']);k=int(c['retainedAfter']);steps+=j;expected+=j+1+(k+1 if k else 0)
            assert c['modelResidualNonincreaseVerified'] and c['actualResidualNorm']<=c['startResidualNorm']+1e-8*max(1,c['startResidualNorm'])
            maximum_retained=max(maximum_retained,int(c['retainedBefore']),k)
            if kind=='unit':
                np.testing.assert_array_equal(kernel['Q'],r['initialQ']);n=A.shape[0]
                Z=np.asarray(c['Z']).reshape((n,j),order='F');V=np.asarray(c['V']).reshape((n,j+1),order='F');H=np.asarray(c['smallH']).reshape((j+1,j),order='F')
                assert np.linalg.norm(A@Z-V@H)/max(1,np.linalg.norm(A@Z))<1e-8
        assert expected==rr['modelApplications']<=rr['maximumModelApplications'] and steps==rr['arnoldiSteps']==rr['preconditionerApplications']
        assert steps<=r['options']['maximumKrylovIterations']
        if a.get('initialGuessUsed',False):
            g=a['initialGuessEvidence'];np.testing.assert_array_equal(a['steadyDirectionEvidence']['directionU'],g['directionU']+g['correctionDirectionU'])
            np.testing.assert_allclose(rr['absoluteResidual'],np.atleast_1d(rr['residualHistory'])[-1],rtol=1e-13,atol=0)
            np.testing.assert_allclose(a['gmresResidualHistory'],rr['residualHistory'],rtol=0,atol=0)
            if kind=='unit':
                correction=g['correctionDirectionU'].ravel(order='F');direct=g['correctionRhs'].ravel(order='F')-A@correction
                np.testing.assert_allclose(np.linalg.norm(direct),rr['absoluteResidual'],rtol=1e-8,atol=1e-8)
        assert a['actualCombinedDirectionCheckedWithOriginalNS'] and a['actualNSResidualDenominator']=='original-base-RHS'
        total_steps+=steps;total_actions+=rr['modelApplications'];total_preconditions+=rr['preconditionerApplications']
    action_rows=rows(r['jvpLog']);recorded_models=sum(bool(x.get('recycledKrylovModelAction',False)) for x in action_rows)
    assert total_actions==r['recycledKrylovModelApplications']==recorded_models and total_steps==r['recycledKrylovArnoldiSteps']
    assert total_preconditions==r['recycledKrylovPreconditionerApplications']==r['preconditionerApplications']==r['coarsePreconditionerModelApplications']
    assert r['selectedAnalyticApplications']==total_actions+r.get('steadyInitialGuessModelApplications',0)
    checks.append(dict(variant=variant,arnoldiSteps=total_steps,modelApplications=total_actions,preconditionerApplications=total_preconditions,
        maximumRetainedDirectionsUsed=maximum_retained,acceptedIterations=int(r['acceptedNonlinearIterations']),actualNSCalls=int(r['residualEvaluations'])))
if kind=='unit':
    assert s['report']['oldSingleGuessCompleteResultReplayed'] and s['report']['basisBudgetRejected']
    assert s['result']['acceptedNonlinearIterations']==1 and s['limited']['candidateTrials']==0 and s['limited']['acceptedNonlinearIterations']==0
out=dict(status='passed-independent-recycled-main-current-images-work-budgets-and-total-direction-identity',case=case,kind=kind,processID=os.getpid(),
    checks=checks,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
