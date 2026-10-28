# -*- coding: utf-8 -*-
"""独立重建小域图矩阵/逐色解，或复核真实图探针的既有作用、着色及叠加。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];kind=sys.argv[2];assert case.replace('_','').isalnum() and kind in ('unit','probe')
target=VAL/(case+'_independent_graph_preconditioner.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
def read(path,names=None):
    bindings[str(path.relative_to(ROOT))]=sha(path);return loadmat(path,simplify_cells=True,variable_names=names)
def rel(a,b):return float(np.max(abs(a-b))/max(1,np.max(abs(b))))
def groups(value):return [value] if isinstance(value,dict) else list(value)
if kind=='unit':
    s=read(VAL/case/'graph_preconditioner_evidence.mat');op=s['operator'];n=int(op['sourceFaces']*op['layers']);nf=int(op['sourceFaces']);M=np.zeros((5*n,5*n))
    for i in range(n):
        rows=np.arange(5*i,5*i+5);M[np.ix_(rows,rows)]=op['diagonal'][:,:,i]
        if i>=nf:M[np.ix_(rows,rows-5*nf)]=op['lower'][:,:,i]
        if i<n-nf:M[np.ix_(rows,rows+5*nf)]=op['upper'][:,:,i]
    m=op['sideMinus'].astype(int)-1;p=op['sidePlus'].astype(int)-1
    for j,(mi,pi) in enumerate(zip(m,p)):
        rows=np.arange(5*mi,5*mi+5);cols=np.arange(5*pi,5*pi+5)
        M[np.ix_(rows,cols)]+=op['minusCoupling'][:,:,j];M[np.ix_(cols,rows)]+=op['plusCoupling'][:,:,j]
    assert rel(M,s['M'])<1e-14 and rel(M,s['fdMatrix'])<5e-8
    color=op['color'].astype(int);assert np.all(color[m%nf]!=color[p%nf])
    # 独立5n矩阵按颜色子系统直接解；不使用MATLAB缓存逆或扫描函数。
    rhs=s['rhs'];u=np.zeros((n,5));gg=groups(op['colorGroups']);order=list(range(len(gg)))+list(range(len(gg)-2,-1,-1))
    for sweep in range(int(op['fixedSymmetricSweeps'])):
        for k in order:
            cells=np.atleast_1d(gg[k]['cellIDs']).astype(int)-1;rows=(5*cells[:,None]+np.arange(5)).ravel();others=np.setdiff1d(np.arange(5*n),rows)
            vector=u.ravel().copy();vector[rows]=np.linalg.solve(M[np.ix_(rows,rows)],rhs.ravel()[rows]-M[np.ix_(rows,others)]@vector[others]);u=vector.reshape(n,5)
    scan_error=rel(u,s['update']);assert scan_error<1e-11
    factor_error=0.;upper_error=0.
    for layer in range(int(op['layers'])):
        for face in range(nf):
            i=layer*nf+face;D=op['diagonal'][:,:,i].copy()
            if layer:D-=op['lower'][:,:,i]@op['lineUpperFactor'][:,:,i-nf]
            factor_error=max(factor_error,rel(D@op['lineInverse'][:,:,i],np.eye(5)))
            upper_error=max(upper_error,rel(np.linalg.solve(D,op['upper'][:,:,i]),op['lineUpperFactor'][:,:,i]))
    assert factor_error<1e-10 and upper_error<1e-10
    linearity=rel(s['combined'],.7*s['update']-1.3*s['u2']);assert linearity<1e-12
    prior=read(VAL/'lateralunit01/lateral_coupling_evidence.mat',['exactSideGraphAction']);assert rel(s['graphAction'],prior['exactSideGraphAction'])<1e-11
    details=dict(cellCount=n,sideFaceCount=len(m),fullMatrixRelativeDifference=rel(M,s['M']),independentMATLABFDRelativeDifference=rel(M,s['fdMatrix']),
        denseColoredScanRelativeDifference=scan_error,cachedInverseRelativeDifference=factor_error,cachedUpperFactorRelativeDifference=upper_error,
        linearSuperpositionRelativeDifference=linearity)
else:
    s=read(VAL/case/'graph_probe_evidence.mat');report=s['report'];prior=read(VAL/'lateralrocket01/lateral_coupling_evidence.mat',['exactSideGraphAction'])
    domain=read(ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat',['domain'])['domain'];nf=len(domain['inputGeometry']['F'])
    f=domain['openFaces'];mask=(f['cellPlus']>0)&(f['kind']==3);m=f['cellMinus'][mask].astype(int)-1;p=f['cellPlus'][mask].astype(int)-1;color=s['color'].astype(int)
    assert np.all(color[m%nf]!=color[p%nf]) and len(np.unique(color))==report['colorCount']
    action_error=rel(s['graphAction'],prior['exactSideGraphAction']);linearity=rel(s['combined'],.7*s['update']-1.3*s['other'])
    assert action_error<2e-11 and linearity<1e-11 and np.all(np.isfinite(s['update']))
    assert not report['fullSparseGraphAllocated'] and not report['peakMemoryMeasured'] and report['newOriginalNSEvaluations']==0
    details=dict(cellCount=len(domain['cells']['volume']),sideFaceCount=len(m),colorCount=int(report['colorCount']),
        allSideColorConflictCount=int(np.count_nonzero(color[m%nf]==color[p%nf])),priorIndependentActionRelativeDifference=action_error,
        linearSuperpositionRelativeDifference=linearity,operatorCacheBytes=float(report['operatorCacheBytes']),fullSparseGraphAllocated=False,peakMemoryMeasured=False)
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-graph-preconditioner-'+kind,case=case,processID=os.getpid(),details=details,
    fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
