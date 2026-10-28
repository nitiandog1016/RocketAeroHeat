# -*- coding: utf-8 -*-
"""独立核对局部策略缓存/新构造、原NS接受与拒绝、新Q重建和模型工作分账。"""
from pathlib import Path
import hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]; VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1]; assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_local_integration.json'); assert not target.exists()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

file=VAL/case/'steady_newton_unit_evidence.mat'; s=loadmat(file,simplify_cells=True)
r=s['result']; fresh=s['fresh']; limited=s['limited']; rebuilt=s['rebuilt']; report=s['report']
assert r['coarseBasis']=='local-conservative-enriched-mixed' and r['coarseUnknowns']==149
assert r['coarseRetainedMixedColumns']==147 and r['coarseLocalConservativeColumns']==2
assert r['coarseOperatorReuses']==1 and r['coarseOperatorBuilds']==0 and r['coarseConstructionModelApplications']==0
assert r['coarseReuseVerificationModelApplications']==1 and r['coarseReuseReport']['sourceAndStateVerified']
assert r['coarseReuseReport']['freshFineModelProbeRelativeDifference']<1e-9
assert fresh['coarseOperatorReuses']==0 and fresh['coarseOperatorBuilds']==1 and fresh['coarseConstructionModelApplications']==149
for a in (r,fresh):
    assert a['acceptedNonlinearIterations']==1 and a['originalNSCombinedDirectionCheckPreserved']
    assert a['coarsePreconditionerModelApplications']==a['preconditionerApplications']
    assert a['linearLog']['fullSystemTargetMet'] and a['linearLog']['target']==.1
    assert a['linearLog']['actualCombinedDirectionCheckedWithOriginalNS'] and a['linearLog']['fullSystemRelativeDefect']<=.1
assert limited['candidateTrials']==0 and limited['acceptedNonlinearIterations']==0 and limited['residualEvaluations']==2
assert limited['options']['linearRelativeTolerance']==1e-7
np.testing.assert_array_equal(limited['Q'],r['initialQ'])
np.testing.assert_array_equal(rebuilt['baseQ'],s['newQ']); np.testing.assert_array_equal(rebuilt['mixedLift']['baseQ'],s['newQ'])
assert not np.array_equal(s['newQ'],r['initialQ']) and rebuilt['constructionModelApplications']==149
G=(rebuilt['image'].T@rebuilt['image']).toarray()
np.testing.assert_allclose(G,rebuilt['normalMatrix'],rtol=5e-10,atol=5e-13)
np.testing.assert_allclose(rebuilt['upper'].T@rebuilt['upper'],G,rtol=5e-10,atol=5e-13)
source=Path(r['options']['coarseOperatorFile']); parent=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse']
for key in ('retainedMixedColumns','localCellIDs','localEquationIDs'):
    np.testing.assert_array_equal(rebuilt[key],parent[key])
difference=rebuilt['image']-parent['image']; difference_norm=float(np.linalg.norm(difference.data)); assert difference_norm>0
np.testing.assert_array_equal(rebuilt['mixedLift']['selectedPrimitiveColumns'],parent['mixedLift']['selectedPrimitiveColumns'])
assert report['oldMixedCompleteResultReplayed'] and report['invalidCacheOrLiftRejections']==8
assert report['integerConfigurationInputTested'] and report['separateNewStateValidationNSEvaluations']==1
assert report['newQProbeRelativeDifference']<1e-9
for name,digest in rebuilt['liftSourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest
bindings={str(p.relative_to(ROOT)):sha(p) for p in (file,source)}
out=dict(status='passed-independent-local-cache-fresh-current-Q-rebuild-and-original-NS-gate-accounting',case=case,
    processID=os.getpid(),cachedActualNSCalls=int(r['residualEvaluations']),freshActualNSCalls=int(fresh['residualEvaluations']),
    strictLimitedActualNSCalls=2,separateNewStateNSEvaluations=1,newStateImageDifferenceNorm=difference_norm,
    freshConstructionModelApplications=149,newStateConstructionModelApplications=149,invalidCacheOrLiftRejections=8,
    oldMixedCompleteResultReplayed=True,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8'); print(json.dumps(out))
