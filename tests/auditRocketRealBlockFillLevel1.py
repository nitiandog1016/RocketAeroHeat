# -*- coding: utf-8 -*-
"""独立重建真实一级符号模式，逐行核对全部保留LU位置及舍弃填充误差。"""
from pathlib import Path
import hashlib,json,os,sys,time
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix,tril,triu
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_real_fill1_factors.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)
factorfile=VAL/case/'block_fill1_factor_evidence.mat';s=read(factorfile,['factor','packedFile','sourceBinding']);f=s['factor'];binding=s['sourceBinding']
p=read(s['packedFile'],['packed','sourceBinding']);assert p['sourceBinding']==binding;p=p['packed'];del s
n=int(p['cellCount']);oldptr=p['rowPointer'].astype(np.int64)-1;oldcol=p['column'].astype(np.int64)-1;pv=p['values'].transpose(2,0,1)
ptr=f['rowPointer'].astype(np.int64)-1;col=f['column'].astype(np.int64)-1;diag=f['diagonalPosition'].astype(np.int64)-1;fv=f['values'].transpose(2,0,1)
P=csr_matrix((np.ones(len(oldcol)),oldcol,oldptr),shape=(n,n));expected=P+tril(P,-1,format='csr')@triu(P,1,format='csr');expected.sort_indices();expected.data[:]=1
np.testing.assert_array_equal(expected.indptr,ptr);np.testing.assert_array_equal(expected.indices,col);del P,expected
np.testing.assert_array_equal(p['order'],f['order']);mapping=f['enrichmentReport']['originalToEnrichedPosition'].astype(np.int64)-1
np.testing.assert_array_equal(col[mapping],oldcol);owners=np.repeat(np.arange(n),np.diff(oldptr));assert np.all((mapping>=ptr[owners])&(mapping<ptr[owners+1]));del owners
e=f['enrichmentReport'];assert len(col)-len(oldcol)==e['addedInitiallyZeroBlocks'] and e['persistentFactorStorageUpperBoundBytes']<=p['options']['maximumPersistentBytes']
assert e['algorithmArrayConstructionEstimateBytes']<=p['options']['maximumConstructionBytes']
retained_error2=outside_error2=den2=0.;retained_count=discarded_count=outside_blocks=0;maximum_backward=inverse_defect=0.;eq_error2=np.zeros(5);eq_norm2=np.zeros(5);timer=time.monotonic()
for i in range(n):
    begin,end=ptr[i],ptr[i+1];columns=col[begin:end];locations={int(c):j for j,c in enumerate(columns)}
    C=np.zeros((end-begin,5,5));old=slice(oldptr[i],oldptr[i+1]);C[mapping[old]-begin]=pv[old]
    product=np.zeros_like(C);bounds=np.linalg.norm(C,axis=(1,2));outside={}
    for lower in range(begin,diag[i]):
        j=int(col[lower]);multiplier=fv[lower]
        for upper in range(diag[j],ptr[j+1]):
            k=int(col[upper]);value=multiplier@fv[upper];where=locations.get(k)
            if where is not None:
                product[where]+=value;bounds[where]+=np.linalg.norm(multiplier)*np.linalg.norm(fv[upper])
                if k>j:retained_count+=1
            else:
                if k>j:discarded_count+=1
                if k in outside:outside[k]+=value
                else:outside[k]=value.copy()
    for upper in range(diag[i],end):
        where=locations[int(col[upper])];product[where]+=fv[upper];bounds[where]+=np.linalg.norm(fv[upper])
    error=product-C;retained_error2+=float(np.sum(error*error));den2+=float(np.sum(C*C));eq_norm2+=np.sum(C*C,axis=(0,2))
    maximum_backward=max(maximum_backward,float(np.max(np.linalg.norm(error,axis=(1,2))/np.maximum(bounds,np.finfo(float).tiny))));eq_error2+=np.sum(error**2,axis=(0,2))
    for value in outside.values():outside_error2+=float(np.sum(value*value));eq_error2+=np.sum(value**2,axis=1)
    outside_blocks+=len(outside);inverse_defect=max(inverse_defect,float(np.max(abs(fv[diag[i]]@f['inverseDiagonal'][:,:,i]-np.eye(5)))))
    if (i+1)%50000==0:print('INDEPENDENT_FILL1_LU rows=%d/%d'%(i+1,n),flush=True)
relative=float(np.sqrt(retained_error2/max(1,den2)));assert relative<1e-8 and maximum_backward<1e-8 and inverse_defect<1e-8
assert retained_count==f['factorizationReport']['retainedBlockProducts'] and discarded_count==f['factorizationReport']['discardedOutOfPatternFillProducts']
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-real-level1-symbolic-pattern-all-retained-LU-products-and-all-pivot-inverses',processID=os.getpid(),
    cellCount=n,originalBlocks=len(oldcol),expandedBlocks=len(col),addedInitiallyZeroBlocks=len(col)-len(oldcol),
    retainedPatternLUProductRelativeDifference=relative,maximumCancellationAwareRetainedBlockBackwardError=maximum_backward,allPivotInverseDefect=inverse_defect,
    fullLUApproximationRelativeDifference=float(np.sqrt((retained_error2+outside_error2)/max(1,den2))),LUApproximationDifferenceByConservationRow=np.sqrt(eq_error2/np.maximum(eq_norm2,np.finfo(float).tiny)).tolist(),
    retainedBlockProducts=retained_count,discardedFillProducts=discarded_count,outsidePatternLUBlocks=outside_blocks,auditSeconds=time.monotonic()-timer,
    correctFactorArithmeticDoesNotCertifyStableInverse=True,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
