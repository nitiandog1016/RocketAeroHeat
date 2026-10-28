# -*- coding: utf-8 -*-
"""独立核对完整五量CSR、块ILU(0)消元与三角解；不运行原NS。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1]
assert case.replace('_','').isalnum();target=VAL/(case+'_independent_block_ilu.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def rows(value):return [value] if isinstance(value,dict) else list(np.atleast_1d(value).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/case/'block_ilu_unit_evidence.mat';s=loadmat(file,simplify_cells=True);n=len(s['Q']);checks=[]
for item in rows(s['cases']):
    packed,factor=item['packed'],item['factor'];ptr=packed['rowPointer'].astype(int)-1;col=packed['column'].astype(int)-1;order=packed['order'].astype(int)-1
    np.testing.assert_array_equal(np.sort(order),np.arange(n));np.testing.assert_array_equal(packed['inverseOrder'].astype(int)[order]-1,np.arange(n))
    indices=(order[:,None]+np.arange(5)*n).ravel();expected=s['C'][np.ix_(indices,indices)];assembled=np.zeros((5*n,5*n));pattern=np.zeros((n,n),bool)
    for i in range(n):
        cc=col[ptr[i]:ptr[i+1]];assert np.all(cc[1:]>cc[:-1]) and i in cc
        for p in range(ptr[i],ptr[i+1]):assembled[5*i:5*i+5,5*col[p]:5*col[p]+5]=packed['values'][:,:,p];pattern[i,col[p]]=True
    assembly_error=rel(assembled,expected);assert assembly_error<2e-13 and np.array_equal(pattern,pattern.T)
    independent=assembled.copy();retained=discarded=0
    for i in range(n):
        ii=slice(5*i,5*i+5)
        for j in np.flatnonzero(pattern[i,:i]):
            jj=slice(5*j,5*j+5);multiplier=np.linalg.solve(independent[jj,jj].T,independent[ii,jj].T).T;independent[ii,jj]=multiplier
            for k in np.flatnonzero(pattern[j,j+1:])+j+1:
                if pattern[i,k]:independent[ii,5*k:5*k+5]-=multiplier@independent[jj,5*k:5*k+5];retained+=1
                else:discarded+=1
    observed=np.zeros_like(assembled);L=np.eye(5*n);U=np.zeros_like(assembled);inverse_error=0.
    for i in range(n):
        ii=slice(5*i,5*i+5)
        for p in range(ptr[i],ptr[i+1]):
            j=col[p];jj=slice(5*j,5*j+5);observed[ii,jj]=factor['values'][:,:,p]
            if j<i:L[ii,jj]=observed[ii,jj]
            else:U[ii,jj]=observed[ii,jj]
        inverse_error=max(inverse_error,rel(factor['inverseDiagonal'][:,:,i],np.linalg.solve(independent[ii,ii],np.eye(5))))
    factor_error=rel(observed,independent);assert factor_error<2e-10 and inverse_error<2e-10
    np.testing.assert_allclose(L,item['L'],rtol=2e-11,atol=1e-10);np.testing.assert_allclose(U,item['U'],rtol=2e-11,atol=1e-10)
    rhs=item['rhs'][order].ravel();actual=item['value'][order].ravel();direct=np.linalg.solve(U,np.linalg.solve(L,rhs));solve_error=rel(actual,direct);assert solve_error<2e-10
    assert rel(item['combined'],.7*item['value']-1.3*item['other'])<2e-10
    assert retained==factor['factorizationReport']['retainedBlockProducts'] and discarded==factor['factorizationReport']['discardedOutOfPatternFillProducts']
    assert factor['packingReport']['persistentFactorStorageUpperBoundBytes']<=factor['options']['maximumPersistentBytes']
    assert not factor['factorizationReport']['diagonalCompensationApplied'] and not factor['factorizationReport']['originalGraphConnectionsDropped']
    checks.append(dict(ordering=packed['options']['ordering'],cells=n,blocks=len(col),assemblyRelativeDifference=assembly_error,
        independentFactorRelativeDifference=factor_error,inverseDiagonalRelativeDifference=inverse_error,triangularSolveRelativeDifference=solve_error,
        retainedBlockProducts=retained,discardedOutOfPatternFillProducts=discarded))
assert s['report']['budgetRejected'] and s['report']['doubleFactorRejected'] and s['report']['newActualNSEvaluations']==0
source=Path(s['sourceFile']);out=dict(status='passed-independent-full-C-CSR-block-ILU0-elimination-and-triangular-reference',case=case,processID=os.getpid(),cases=checks,
    newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,source)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
