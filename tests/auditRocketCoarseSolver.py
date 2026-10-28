# -*- coding: utf-8 -*-
"""独立核对两层接入的模型/原NS调用分账及缓存来源；实际方向由另一审计复算。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];kind=sys.argv[2];assert case.replace('_','').isalnum() and kind in ('unit','real')
target=VAL/(case+'_independent_coarse_solver.json');assert not target.exists()
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
file=VAL/case/('steady_newton_unit_evidence.mat' if kind=='unit' else 'wall_velocity_krylov_evidence.mat')
bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True);r=s['result'];rows=lambda x:[x] if isinstance(x,dict) else list(x)
assert r['rightPreconditioner']=='two-level-minres' and r['jacobianAction']=='selected-analytic' and r['linearEquation']=='steady-newton'
assert r['fixedLinearPreconditioner'] and r['originalNSCombinedDirectionCheckPreserved'] and not r['coarseCorrectionChangesPhysicalMesh']
assert r['coarseOperatorBuilds']==0 and r['coarseOperatorReuses']==1 and r['coarseConstructionModelApplications']==0
assert r['coarsePreconditionerModelApplications']==r['preconditionerApplications'] and r['coarseReuseVerificationModelApplications']==1
reuse=r['coarseReuseReport'];assert reuse['sourceAndStateVerified'] and reuse['freshFineModelProbeRelativeDifference']<1e-9
source=Path(r['options']['coarseOperatorFile']);source=source if source.is_absolute() else ROOT/source;bindings[str(source.relative_to(ROOT))]=sha(source)
coarse=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse']
np.testing.assert_array_equal(coarse['baseQ'],r['initialQ']);np.testing.assert_array_equal(coarse['scale'],r['residualScale'])
additional=len(np.atleast_1d(coarse['selectedPrimitiveColumns'])) if 'selectedPrimitiveColumns' in coarse else 0
local='localCellIDs' in coarse
count=len(np.atleast_1d(coarse['retainedMixedColumns']))+len(np.atleast_1d(coarse['localCellIDs'])) if local else 5*coarse['space']['groupCount']+additional
assert coarse['imageBytes']==r['coarseImageBytes']<=r['options']['coarseMaximumImageBytes'] and count==r['coarseUnknowns']
if local:
    np.testing.assert_array_equal(coarse['mixedLift']['selectedPrimitiveColumns'],r['options']['coarseMixedColumns'])
    for saved,key in [('retainedMixedColumns','retainedMixedColumns'),('localCellIDs','cellIDs'),('localEquationIDs','equationIDs')]:
        np.testing.assert_array_equal(coarse[saved],r['options']['coarseLocalColumns'][key])
    for name,digest in coarse['liftSourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest
if additional:
    np.testing.assert_array_equal(np.atleast_1d(coarse['selectedPrimitiveColumns']),np.atleast_1d(r['options']['coarseMixedColumns']))
    for name,digest in coarse['liftSourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest
for name,digest in coarse['sourceIdentity'].items():
    path=ROOT/'internalFlow'/('' if name=='resolveNSCharacteristicBoundary3D' else 'cutCell')/(name+'.m');assert sha(path)==digest
del coarse
logs=rows(r['linearLog']);actions=rows(r['jvpLog']);models=[x for x in actions if x['scheme']=='selected-analytic'];checks=[x for x in actions if x['scheme']!='selected-analytic']
assert len(models)==r['selectedAnalyticApplications'] and len(checks)==len(logs)
assert sum(x['actualOriginalResidualEvaluations'] for x in actions)==r['jacobianResidualEvaluations']
for a in logs:
    assert a['actualCombinedDirectionCheckedWithOriginalNS'] and a['coarseAlgebraicRelativeResidual']<1e-8
    assert a['fineModelResidualNonincreaseVerified'] and a['coarseUnknowns']==r['coarseUnknowns']
    assert a['preconditionerResidualIsInitialGraphScan'] and not a['originalNSLinearAccuracyClaimed']
    assert a['coarsePredictedFineRelativeResidual']<=1+1e-9
details=dict(coarseUnknowns=int(r['coarseUnknowns']),coarseImageBytes=int(r['coarseImageBytes']),
    outerAnalyticApplications=len(models),preconditionerModelApplications=int(r['coarsePreconditionerModelApplications']),
    constructionModelApplications=0,reuseVerificationModelApplications=1,originalNSDirectionChecks=len(checks),
    totalSegmentNSEvaluations=int(r['residualEvaluations']),reuseProbeRelativeDifference=float(reuse['freshFineModelProbeRelativeDifference']),
    actualLinearDefects=[float(a['fullSystemRelativeDefect']) for a in logs])
if kind=='unit':
    fresh=s['fresh'];assert fresh['coarseOperatorBuilds']==1 and fresh['coarseOperatorReuses']==0 and fresh['coarseConstructionModelApplications']==r['coarseUnknowns']
    np.testing.assert_array_equal(fresh['Q'],r['Q']);np.testing.assert_array_equal(fresh['rate'],r['rate'])
    limited=s['limited'];assert limited['candidateTrials']==0 and limited['acceptedNonlinearIterations']==0 and limited['residualEvaluations']==2
    assert r['acceptedNonlinearIterations']==1 and r['residualEvaluations']==3 and s['report']['invalidCacheOrPolicyRejections']==11
    details.update(freshAndReusedEndpointIdentical=True,invalidCacheOrPolicyRejections=11)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-two-level-cache-identity-and-model-NS-work-accounting',case=case,kind=kind,processID=os.getpid(),details=details,
    fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
