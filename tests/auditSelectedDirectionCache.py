# -*- coding: utf-8 -*-
"""绑定原已核对方向，逐值检查缓存模型、线性叠加与零新NS调用身份。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_direction_cache.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'selected_direction_cache_evidence.mat';s=loadmat(file,simplify_cells=True);source=Path(s['sourceFile']);bindings={str(p.relative_to(ROOT)):sha(p) for p in (file,source)}
original=loadmat(source,simplify_cells=True,variable_names=['allEvidence'])['allEvidence'];items=[original] if isinstance(original,dict) else list(original)
def matrices(value,count):return [value] if count==1 else list(value)
directions=matrices(s['directions'],len(items));expected=matrices(s['expected'],len(items));values=matrices(s['values'],len(items));errors=[]
for item,d,reference,value in zip(items,directions,expected,values):
    np.testing.assert_array_equal(d,item['directionU']);np.testing.assert_array_equal(reference,item['viscous']+item['euler'])
    error=float(np.linalg.norm(value-reference)/max(1,np.linalg.norm(reference)));errors.append(error);assert error<1e-10
linear=.7*values[0]-1.3*s['second'];linearity=float(np.linalg.norm(s['combined']-linear)/max(1,np.linalg.norm(linear)));assert linearity<1e-10
assert s['report']['newOriginalNSEvaluations']==0 and s['report']['newAcceptedIterations']==0 and not s['report']['physicalSolution']
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-selected-direction-cache-source-values-and-linear-superposition',case=case,processID=os.getpid(),
    coldActionRelativeDifferences=errors,linearSuperpositionRelativeDifference=linearity,newOriginalNSEvaluations=0,newAcceptedIterations=0,
    physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
