# -*- coding: utf-8 -*-
"""独立核对同Q初猜、当前模型剩余、原分母容差和总方向原NS证据。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];kind=sys.argv[2]
assert case.replace('_','').isalnum() and kind in ('unit','real')
revision=sys.argv[3] if len(sys.argv)>3 else '';assert not revision or revision.isdigit()
target=VAL/(case+'_independent_steady_guess'+revision+'.json');assert not target.exists()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

file=VAL/case/('steady_newton_unit_evidence.mat' if kind=='unit' else 'wall_velocity_krylov_evidence.mat')
s=loadmat(file,simplify_cells=True);bindings={str(file.relative_to(ROOT)):sha(file)};checks=[]
rows=lambda x:[x] if isinstance(x,dict) else list(np.atleast_1d(x))
variants=['result','ready','limited'] if kind=='unit' else ['result']
for variant in variants:
    r=s[variant];warm=[a for a in rows(r['linearLog']) if a.get('initialGuessUsed',False)]
    assert warm and len(warm)==r['steadyInitialGuessModelApplications']
    model_rows=[a for a in rows(r['jvpLog']) if a.get('initialGuessModelAction',False)]
    assert len(model_rows)==len(warm) and all(a['scheme']=='selected-analytic' and a['actualOriginalResidualEvaluations']==0 for a in model_rows)
    for a in warm:
        g=a['initialGuessEvidence'];source=g['sourceReport'];path=Path(source['sourceFile']);coarse=Path(source['sourceCoarseFile'])
        if not path.is_absolute():path=ROOT/path
        if not coarse.is_absolute():coarse=ROOT/coarse
        assert sha(path)==source['sourceFileSHA256'] and sha(coarse)==source['sourceCoarseFileSHA256']
        assert source['sourceStateAndGeometryVerified'] and not source['oldDerivativeReused'] and not source['originalNSAcceptanceClaimed']
        bindings[str(path.relative_to(ROOT))]=sha(path);bindings[str(coarse.relative_to(ROOT))]=sha(coarse)
        parent=loadmat(path,simplify_cells=True,variable_names=[source['sourceVariant']])[source['sourceVariant']]
        np.testing.assert_array_equal(parent['initialQ'],r['initialQ']);np.testing.assert_array_equal(parent['Q'],r['initialQ'])
        np.testing.assert_array_equal(parent['residualScale'],r['residualScale']);assert parent['acceptedNonlinearIterations']==0
        old=rows(parent['linearLog'])[-1]['steadyDirectionEvidence']
        np.testing.assert_array_equal(g['directionU'],old['directionU']*source['directionScale']);del parent
        rhs=r['initialRate']/r['residualScale'];np.testing.assert_array_equal(g['correctionRhs'],rhs+g['currentModelValue'])
        np.testing.assert_array_equal(a['steadyDirectionEvidence']['directionU'],g['directionU']+g['correctionDirectionU'])
        defect=float(np.linalg.norm(g['correctionRhs'])/np.linalg.norm(rhs))
        np.testing.assert_allclose(a['initialGuessModelRelativeDefect'],defect,rtol=4e-13,atol=0)
        expected=min(1-np.finfo(float).eps,r['options']['linearRelativeTolerance']/max(defect,np.finfo(float).tiny))
        np.testing.assert_allclose(a['effectiveCorrectionTolerance'],expected,rtol=4e-13,atol=0)
        assert a['actualNSResidualDenominator']=='original-base-RHS' and a['gmresResidualDenominator']=='current-initial-guess-correction-RHS'
        history=np.atleast_1d(a['gmresResidualHistory']);np.testing.assert_allclose(history[0],np.linalg.norm(g['correctionRhs']),rtol=4e-13,atol=0)
        actual=float(np.linalg.norm(rhs+a['steadyDirectionEvidence']['value'])/np.linalg.norm(rhs))
        np.testing.assert_allclose(actual,a['fullSystemRelativeDefect'],rtol=4e-13,atol=0)
        assert a['target']==r['options']['linearRelativeTolerance']
        if a['gmresSkipped']:
            assert np.isnan(a['gmresReportedResidual']) and len(history)==1 and defect<=a['target']
        else:
            np.testing.assert_allclose(a['gmresReportedResidual'],history[-1]/history[0],rtol=2e-10,atol=1e-12)
        checks.append(dict(variant=variant,initialModelDefect=defect,effectiveCorrectionTolerance=expected,actualTotalNSDefect=actual,
            originalTarget=float(a['target']),gmresSkipped=bool(a['gmresSkipped']),acceptedFlowUpdates=int(r['acceptedNonlinearIterations'])))
if kind=='unit':
    assert s['report']['oldLocalCompleteResultReplayed'] and s['report']['invalidSourceRejections']==8
    assert not s['result']['linearLog']['gmresSkipped'] and s['ready']['linearLog']['gmresSkipped']
    assert s['limited']['candidateTrials']==0 and s['limited']['acceptedNonlinearIterations']==0
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-same-Q-guess-current-remainder-total-direction-and-original-normalization',case=case,kind=kind,
    processID=os.getpid(),checks=checks,oldDerivativeReused=False,incrementOnlyAcceptanceUsed=False,
    completeIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
