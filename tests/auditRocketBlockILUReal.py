# -*- coding: utf-8 -*-
"""逐块核对真实完整C与ILU0的全部原模式系数，并计算舍弃填充误差。"""
from pathlib import Path
import hashlib,json,os,sys,time
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_real_C_and_factors.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)
def rows(value):return [value] if isinstance(value,dict) else list(np.atleast_1d(value).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
folder=VAL/case;p=read(folder/'steady_graph_packed_evidence.mat');packed=p['packed'];f=read(folder/'block_ilu_factor_evidence.mat');factor=f['factor']
assert p['sourceBinding']==f['sourceBinding'] and p['sourceIdentity']==f['sourceIdentity']
source=Path(p['sourceBinding']['sourceFile']);assert sha(source)==p['sourceBinding']['sourceFileSHA256'];r=read(source,['result'])['result'];op=r['firstOperator']
np.testing.assert_array_equal(f['Q'],r['Q']);np.testing.assert_array_equal(f['scale'],r['residualScale']);del r
op.pop('wlsNormalFactor',None);op.pop('colorGroups',None)
for name,digest in f['sourceIdentity'].items():
    current=ROOT/'internalFlow/cutCell'/(name+'.m')
    if sha(current)!=digest:
        epoch=VAL/'source_epochs'/('pre_blockilu_failure_'+name+'.m');assert epoch.is_file() and sha(epoch)==digest
n=int(packed['cellCount']);nf=int(packed['sourceFaces']);ptr=packed['rowPointer'].astype(np.int64)-1;col=packed['column'].astype(np.int64)-1;diag=packed['diagonalPosition'].astype(np.int64)-1
order=packed['order'].astype(np.int64)-1;inv=packed['inverseOrder'].astype(np.int64)-1;nb=len(col)
np.testing.assert_array_equal(np.sort(order),np.arange(n));np.testing.assert_array_equal(inv[order],np.arange(n))
for key in ('rowPointer','column','diagonalPosition','order','inverseOrder'):np.testing.assert_array_equal(packed[key],factor[key])
assert ptr[0]==0 and ptr[-1]==nb and np.all(col[diag]==np.arange(n))
position=csr_matrix((np.arange(nb,dtype=np.int64)+1,col.astype(np.int32),ptr.astype(np.int32)),shape=(n,n))
expected=np.zeros((nb,5,5));covered=np.zeros(nb,bool);raw=0
def add(physical_rows,physical_cols,pages,page_ids=None):
    global raw
    raw+=len(physical_rows)
    for first in range(0,len(physical_rows),4096):
        rr=physical_rows[first:first+4096];cc=physical_cols[first:first+4096]
        where=np.asarray(position[inv[rr],inv[cc]]).ravel()-1;assert np.all(where>=0)
        ids=slice(first,first+len(rr)) if page_ids is None else page_ids[first:first+4096]
        np.add.at(expected,where,pages[:,:,ids].transpose(2,0,1));covered[where]=True
self_ids=np.arange(n);normal=np.arange(n-nf);second=np.arange(n-2*nf);D=op['diagonal'].copy()
for a in range(5):D[a,a]-=1/op['pseudoTimeStepSeconds']
add(self_ids,self_ids,D);del D
add(normal,normal+nf,op['upper'],normal);add(normal+nf,normal,op['lower'],normal+nf)
add(second,second+2*nf,op['secondUpper'],second);add(second+2*nf,second,op['secondLower'],second+2*nf)
m=op['sideMinus'].astype(np.int64)-1;other=op['sidePlus'].astype(np.int64)-1
add(m,other,op['minusCoupling']);add(other,m,op['plusCoupling']);assert np.all(covered)
assert raw==packed['packingReport']['rawBlockContributions'] and nb==packed['packingReport']['uniqueStoredBlocks']
pv=packed['values'].transpose(2,0,1);fv=factor['values'].transpose(2,0,1);err2=den2=0.;equation_norm2=np.zeros(5)
for first in range(0,nb,4096):
    target_values=expected[first:first+4096];difference=target_values-pv[first:first+4096]
    err2+=float(np.sum(difference*difference));den2+=float(np.sum(target_values*target_values));equation_norm2+=np.sum(target_values**2,axis=(0,2))
assembly_error=float(np.sqrt(err2/max(1,den2)));assert assembly_error<2e-11
del expected,position,covered,op
owners=np.repeat(np.arange(n),np.diff(ptr));probe_errors=[]
for probe in rows(p['probes']):
    vector=probe['probe'][order];product=np.zeros((n,5))
    for first in range(0,nb,8192):
        selected=slice(first,first+8192);np.add.at(product,owners[selected],np.einsum('kab,kb->ka',pv[selected],vector[col[selected]]))
    result=np.zeros((n,5));result[order]=product;error=rel(result,probe['expected']);assert error<2e-10 and rel(result,probe['packedValue'])<2e-10;probe_errors.append(error)
del owners
inverse_defect=0.
for first in range(0,n,4096):
    ids=slice(first,min(first+4096,n));D=fv[diag[ids]];inverse=factor['inverseDiagonal'][:,:,ids].transpose(2,0,1)
    inverse_defect=max(inverse_defect,float(np.max(abs(D@inverse-np.eye(5)))))
assert inverse_defect<1e-8
retained_error2=outside_error2=0.;retained_count=discarded_count=outside_blocks=0;maximum_backward=0.;equation_error2=np.zeros(5);timer=time.monotonic()
for i in range(n):
    begin,end=ptr[i],ptr[i+1];columns=col[begin:end];assert np.all(columns[1:]>columns[:-1]);locations={int(c):j for j,c in enumerate(columns)}
    product=np.zeros((end-begin,5,5));bounds=np.linalg.norm(pv[begin:end],axis=(1,2));outside={}
    for lower in range(begin,diag[i]):
        j=int(col[lower]);assert j<i;multiplier=fv[lower]
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
    error=product-pv[begin:end];retained_error2+=float(np.sum(error*error))
    maximum_backward=max(maximum_backward,float(np.max(np.linalg.norm(error,axis=(1,2))/np.maximum(bounds,np.finfo(float).tiny))))
    equation_error2+=np.sum(error**2,axis=(0,2))
    for value in outside.values():outside_error2+=float(np.sum(value*value));equation_error2+=np.sum(value**2,axis=1)
    outside_blocks+=len(outside)
    if (i+1)%50000==0:print('INDEPENDENT_BLOCK_LU rows=%d/%d'%(i+1,n),flush=True)
retained_error=float(np.sqrt(retained_error2/max(1,den2)));assert retained_error<1e-8 and maximum_backward<1e-8
assert retained_count==factor['factorizationReport']['retainedBlockProducts'] and discarded_count==factor['factorizationReport']['discardedOutOfPatternFillProducts']
assert packed['packingReport']['persistentFactorStorageUpperBoundBytes']<=packed['options']['maximumPersistentBytes']
assert packed['packingReport']['algorithmArrayConstructionEstimateBytes']<=packed['options']['maximumConstructionBytes']
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-real-complete-C-all-blocks-and-all-retained-ILU0-product-audit',case=case,processID=os.getpid(),cellCount=n,unknowns=5*n,
    rawBlockContributions=raw,uniqueStoredBlocks=nb,allCBlockAssemblyRelativeDifference=assembly_error,fullVectorProbeRelativeDifferences=probe_errors,
    allPivotInverseDefect=inverse_defect,allRetainedBlockLUProductRelativeDifference=retained_error,maximumCancellationAwareRetainedBlockBackwardError=maximum_backward,
    outsidePatternLUBlocks=outside_blocks,fullLUApproximationRelativeDifference=float(np.sqrt((retained_error2+outside_error2)/max(1,den2))),
    LUApproximationDifferenceByConservationRow=np.sqrt(equation_error2/np.maximum(equation_norm2,np.finfo(float).tiny)).tolist(),
    retainedBlockProducts=retained_count,discardedFillProducts=discarded_count,auditProductSeconds=time.monotonic()-timer,
    numericFactorStorageUpperBoundBytes=int(packed['packingReport']['persistentFactorStorageUpperBoundBytes']),
    matrixAndFactorIdentityDoesNotProveStableInverse=True,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
