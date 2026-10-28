# -*- coding: utf-8 -*-
"""只读实际可变块方向，按精确恒等式分开旧粗近似缺陷与新F映射不一致贡献。"""
from pathlib import Path
import os,json,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha,arrays,local

name=sys.argv[1];assert name.replace('_','').isalnum();folder=VAL/name
target=folder/'flexible_schur_defect_diagnostic.json';assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rect(g,x):
    owner=np.atleast_1d(g['owner']).astype(np.int64)-1;col=np.atleast_1d(g['column']).astype(np.int64)-1
    blocks=np.asarray(g['values']).reshape(5,5,len(col),order='F').transpose(2,0,1)
    output=np.zeros((int(g['rowCount']),5))
    for first in range(0,len(col),4096):
        ii=slice(first,first+4096);np.add.at(output,owner[ii],np.einsum('kab,kb->ka',blocks[ii],x[col[ii]]))
    return output


def schur(g,x):
    n=int(g['cellCount']);ptr=np.atleast_1d(g['rowPointer']).astype(np.int64)-1
    col=np.atleast_1d(g['column']).astype(np.int64)-1;order=np.atleast_1d(g['order']).astype(np.int64)-1
    owner=np.repeat(np.arange(n),np.diff(ptr));blocks=np.asarray(g['values']).reshape(5,5,len(col),order='F').transpose(2,0,1)
    output=np.zeros((n,5));ordered=x[order]
    for first in range(0,len(col),4096):
        ii=slice(first,first+4096);np.add.at(output,owner[ii],np.einsum('kab,kb->ka',blocks[ii],ordered[col[ii]]))
    value=np.empty_like(output);value[order]=output;return value


s=read(folder/'flexible_reduction_real_evidence.mat',['assemblyFile','detail'])
assert sha(s['assemblyFile'])=='b04c5c8d604c52c9289ce5d08127270a41c31322eef5ac44cd0fda4ebc0c76ca'
r=read(s['assemblyFile'],['reduction'])['reduction'];plan=arrays(r['Ffactor'])
f=np.atleast_1d(r['split']['fineOriginalCVs']).astype(np.int64)-1;c=np.atleast_1d(r['split']['coarseOriginalCVs']).astype(np.int64)-1
reports=s['detail']['preconditionerReports'];reports=[reports] if isinstance(reports,dict) else list(np.atleast_1d(reports).ravel())
results=[]
for j,p in enumerate(reports):
    d=read(p['options']['evidenceFile'],['rhs','firstF','secondF','xC','coarseRHS'])
    bx=rect(r['B'],d['xC']);mf_bx=local(plan,bx.ravel()).reshape(len(f),5)
    coarse_defect=d['coarseRHS']-schur(r['schur'],d['xC'])
    mismatch=rect(r['E'],d['firstF']-d['secondF']-mf_bx)
    c_remaining=d['rhs'][c]-rect(r['E'],d['secondF'])-rect(r['D'],d['xC'])
    denominator=max(np.finfo(float).tiny,np.linalg.norm(d['rhs']))
    relation=float(np.linalg.norm(c_remaining-coarse_defect-mismatch)/denominator)
    assert relation<1e-8
    results.append(dict(application=j+1,oldSchurLocalEquationDefectRelativeToFullRHS=float(np.linalg.norm(coarse_defect)/denominator),
        changedFMapContributionRelativeToFullRHS=float(np.linalg.norm(mismatch)/denominator),
        actualParentCRowsRelativeToFullRHS=float(np.linalg.norm(c_remaining)/denominator),
        exactDecompositionRelativeDifference=relation))
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-read-only-exact-coarse-defect-and-changed-F-map-decomposition-for-actual-flexible-directions',processID=os.getpid(),applications=results,
    identity='b_C-E*x_F-D*x_C = (b_C-E*f_first-S_old*x_C) + E*(f_first-x_F-M_old*B*x_C)',
    oldSNotClaimedToMatchNewFMap=True,uniqueGlobalFailureCauseProven=False,
    separateReadOnlyOldSApplications=len(results),separateReadOnlyOldFLocalApplications=len(results),
    separateReadOnlyRectangularApplications=4*len(results),newOriginalCApplications=0,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
