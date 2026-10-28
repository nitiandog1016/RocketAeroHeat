# -*- coding: utf-8 -*-
"""核对混合接入的缓存/新构造身份、新状态像与原NS工作分账。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_mixed_integration.json');assert not target.exists()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'steady_newton_unit_evidence.mat';digest=sha(file);s=loadmat(file,simplify_cells=True);r=s['result'];fresh=s['fresh'];rebuilt=s['rebuilt'];report=s['report']
assert r['coarseBasis']=='mixed-conservative-primitive' and r['coarseUnknowns']==147 and r['coarseBaseConservativeColumns']==140 and r['coarseAdditionalPrimitiveColumns']==7
assert r['coarseOperatorReuses']==1 and r['coarseOperatorBuilds']==0 and r['coarseConstructionModelApplications']==0
assert fresh['coarseOperatorReuses']==0 and fresh['coarseOperatorBuilds']==1 and fresh['coarseConstructionModelApplications']==147
for a in [r,fresh]:
    assert a['coarsePreconditionerModelApplications']==a['preconditionerApplications'] and a['acceptedNonlinearIterations']==1 and a['residualEvaluations']==3
    assert a['originalNSCombinedDirectionCheckPreserved'] and a['linearLog']['actualCombinedDirectionCheckedWithOriginalNS']
assert s['limited']['candidateTrials']==0 and s['limited']['acceptedNonlinearIterations']==0 and s['limited']['residualEvaluations']==2
np.testing.assert_array_equal(rebuilt['baseQ'],s['newQ']);assert not np.array_equal(s['newQ'],r['initialQ'])
newQ=s['newQ'];gas=s['gas'];rho=newQ[:,0];u=newQ[:,1:4]/rho[:,None];T=(newQ[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*gas['R']/(gas['gamma']-1))
np.testing.assert_allclose(s['newInfo']['primitive'],np.column_stack((rho,u,T)),rtol=5e-14,atol=0)
G=(rebuilt['image'].T@rebuilt['image']).toarray();np.testing.assert_allclose(G,rebuilt['normalMatrix'],rtol=5e-11,atol=5e-13)
np.testing.assert_allclose(rebuilt['upper'].T@rebuilt['upper'],G,rtol=5e-11,atol=5e-13);assert rebuilt['constructionModelApplications']==147
source=Path(r['options']['coarseOperatorFile']);before=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse']
difference=rebuilt['image']-before['image'];assert np.linalg.norm(difference.data)>0
for name,digest2 in rebuilt['liftSourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest2
assert report['invalidCacheOrLiftRejections']==4 and report['separateNewStateValidationNSEvaluations']==1
assert sha(file)==digest
out=dict(status='passed-independent-mixed-cache-fresh-build-new-state-identity-and-work-accounting',case=case,processID=os.getpid(),
    reusedInitialBasis=True,freshConstructionModelApplications=147,newStateConstructionModelApplications=147,newStateImageDifferenceNorm=float(np.linalg.norm(difference.data)),
    actualOriginalNSCallsCached=3,actualOriginalNSCallsFresh=3,actualOriginalNSCallsStrictLimited=2,separateNewStateNSEvaluations=1,
    invalidCacheOrLiftRejections=4,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):digest,str(source.relative_to(ROOT)):sha(source)},scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
