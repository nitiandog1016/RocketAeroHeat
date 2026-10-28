# -*- coding: utf-8 -*-
"""独立稠密一级模式、零新增块、ILU算术及三角解核对；不重构物理Jacobian。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_fill1_unit.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def pattern(p):
    n=int(p['cellCount']);mask=np.zeros((n,n),bool);ptr=p['rowPointer'].astype(int)-1;cols=p['column'].astype(int)-1
    for i in range(n):mask[i,cols[ptr[i]:ptr[i+1]]]=True
    return mask,ptr,cols
file=VAL/case/'block_fill1_unit_evidence.mat';s=loadmat(file,simplify_cells=True);results=[]
for c in np.atleast_1d(s['cases']).ravel():
    old=c['original'];p=c['enriched'];f=c['factor'];base,_,_=pattern(old);mask,ptr,cols=pattern(p);n=len(base)
    expected=(base.astype(int)+np.tril(base,-1).astype(int)@np.triu(base,1).astype(int))>0;np.testing.assert_array_equal(mask,expected)
    mapping=p['enrichmentReport']['originalToEnrichedPosition'].astype(int)-1
    np.testing.assert_array_equal(p['values'][:,:,mapping],old['values']);added=np.ones(len(cols),bool);added[mapping]=False;assert np.all(p['values'][:,:,added]==0)
    order=p['order'].astype(int)-1;indices=(order[:,None]+np.arange(5)*n).ravel();C=s['C'][np.ix_(indices,indices)];reference=C.copy()
    for i in range(n):
        rows=slice(5*i,5*i+5)
        for j in np.flatnonzero(mask[i,:i]):
            js=slice(5*j,5*j+5);multiplier=np.linalg.solve(reference[js,js].T,reference[rows,js].T).T;reference[rows,js]=multiplier
            for k in np.flatnonzero(mask[j,j+1:])+j+1:
                if mask[i,k]:reference[rows,5*k:5*k+5]-=multiplier@reference[js,5*k:5*k+5]
    stored=np.zeros_like(C);L=np.eye(5*n);U=np.zeros_like(C)
    for i in range(n):
        for pidx in range(ptr[i],ptr[i+1]):
            j=cols[pidx];block=f['values'][:,:,pidx];stored[5*i:5*i+5,5*j:5*j+5]=block
            if j<i:L[5*i:5*i+5,5*j:5*j+5]=block
            else:U[5*i:5*i+5,5*j:5*j+5]=block
    error=float(np.linalg.norm(stored-reference)/max(1,np.linalg.norm(reference)));assert error<1e-8
    rhs=c['rhs'].ravel(order='F')[indices];expected_x=np.linalg.solve(U,np.linalg.solve(L,rhs));observed=c['solution'].ravel(order='F')[indices]
    solve_error=float(np.linalg.norm(observed-expected_x)/max(1,np.linalg.norm(expected_x)));assert solve_error<1e-8
    product=L@U;backward=float(np.linalg.norm((product-C)*np.repeat(np.repeat(mask,5,axis=0),5,axis=1))/max(1,np.linalg.norm(C)));assert backward<1e-8
    results.append(dict(ordering=old['options']['ordering'],originalBlocks=int(base.sum()),expandedBlocks=int(mask.sum()),addedInitiallyZeroBlocks=int(added.sum()),
                        independentFactorRelativeDifference=error,denseTriangularSolveDifference=solve_error,retainedPatternLUProductDifference=backward))
assert s['report']['budgetRejected'] and s['report']['doubleEnrichmentRejected']
out=dict(status='passed-independent-small-symbolic-level1-all-C-values-factors-and-triangular-solves',processID=os.getpid(),cases=results,
         newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
         inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
