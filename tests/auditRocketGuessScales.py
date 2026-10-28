# -*- coding: utf-8 -*-
"""独立复算预条件尺度探针的最佳标量、当前剩余、模型投影及未更新身份。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1]
assert case.replace('_','').isalnum();target=VAL/(case+'_independent_guess_scales.json');assert not target.exists()
file=VAL/case/'guess_scale_evidence.mat';s=loadmat(file,simplify_cells=True);report=s['report'];H=s['history'];remaining=s['remaining'];norm=np.linalg.norm

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

source=Path(s['sourceFile']);source=source if source.is_absolute() else ROOT/source;r=loadmat(source,simplify_cells=True,variable_names=['result'])['result'];a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
np.testing.assert_array_equal(remaining,r['initialRate']/r['residualScale']+a['steadyDirectionEvidence']['value']);del r
checks=[]
for j,row in enumerate(H):
    if not np.isfinite(row[3]):continue
    action=s['products'][j];direction=s['directions'][j]
    scalar=float(np.sum(remaining*action)/np.sum(action*action));prediction=remaining-scalar*action
    np.testing.assert_allclose(row[2],scalar,rtol=1e-10,atol=1e-14)
    np.testing.assert_allclose(row[3],norm(prediction)/norm(s['rhs']),rtol=4e-13,atol=0)
    np.testing.assert_allclose(row[4],1-norm(prediction)**2/norm(remaining)**2,rtol=1e-6,atol=2e-14)
    assert np.isfinite(direction).all() and norm(prediction)<=norm(remaining)+1e-8
    checks.append(dict(pseudoDt=float(row[0]),graphRelativeResidual=float(row[1]),oneDirectionPredictedDefect=float(row[3]),squaredNormFractionRemoved=float(row[4])))
assert report['newActualNSEvaluations']==0 and report['newFlowAcceptances']==0 and report['newCoarseConstructionApplications']==0
assert report['modelApplications']==2*len(checks) and len(checks)==6
assert abs(report['bestOneDirectionPredictedDefect']-min(x['oneDirectionPredictedDefect'] for x in checks))<1e-14
out=dict(status='passed-independent-current-remaining-preconditioner-scale-projection-screen',case=case,processID=os.getpid(),checks=checks,
    modelApplications=int(report['modelApplications']),newActualNSEvaluations=0,newFlowAcceptances=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,source)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
