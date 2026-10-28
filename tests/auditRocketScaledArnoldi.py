# -*- coding: utf-8 -*-
"""独立验证整方向缩放后的物理Z/完整像、小域稠密解和未改的默认/拒绝合同。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_scaled_arnoldi.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
def verify(A,b,value,krylov,expected_Z):
    c=krylov['cycles'];j=int(c['arnoldiSteps']);V=c['V'];Z=c['Z'];H=c['smallH'];scales=np.atleast_1d(c['imageNormalizationScales'])
    assert c['basisOrthogonalityError']<1e-8 and c['arnoldiRelationRelativeDifference']<1e-8
    assert krylov['originalOrthogonalityThreshold']==krylov['originalModelRelationThreshold']==1e-8
    ze=rel(Z,expected_Z(V,scales));ae=rel(A@Z,V@H);assert max(ze,ae)<2e-10
    weights=np.linalg.lstsq(H,np.r_[np.linalg.norm(b),np.zeros(j)],rcond=1e-13)[0]
    assert rel(H@weights,H@np.atleast_1d(c['weights']))<2e-10
    residual=np.linalg.norm(b-A@value)/np.linalg.norm(b);assert residual<1.1e-10
    assert rel(c['actualRemaining'],b-A@c['candidate'])<2e-10
    assert krylov['modelApplications']==j+1 and krylov['preconditionerApplications']==j
    assert np.max(c['reorthogonalizationPasses'])<=8 and krylov['basisStorageUpperBoundBytes']<=krylov['maximumBasisBytes']
    return dict(relativeResidual=float(residual),physicalZRelativeDifference=ze,completeArnoldiRelationDifference=ae,arnoldiSteps=j,
                imageScaleRange=[float(min(scales)),float(max(scales))],maximumReorthogonalizationPasses=int(max(c['reorthogonalizationPasses'])))
file=VAL/case/'scaled_arnoldi_unit_evidence.mat';s=loadmat(file,simplify_cells=True)
toy=verify(s['A'],s['b'],s['value'],s['detail'],lambda V,d:V[:,:len(d)]*s['gains'][:len(d)][None,:]/d[None,:])
f=s['factor'];n=int(f['cellCount']);L=np.eye(5*n);U=np.zeros((5*n,5*n));ptr=f['rowPointer'].astype(int)-1;cols=f['column'].astype(int)-1
order=f['order'].astype(int)-1;indices=(order[:,None]+np.arange(5)*n).ravel()
for i in range(n):
    for p in range(ptr[i],ptr[i+1]):
        j=cols[p]
        if j<i:L[5*i:5*i+5,5*j:5*j+5]=f['values'][:,:,p]
        else:U[5*i:5*i+5,5*j:5*j+5]=f['values'][:,:,p]
def P(V,d):
    out=np.zeros((5*n,len(d)));out[indices]=np.linalg.solve(U,np.linalg.solve(L,V[indices,:len(d)]))/d[None,:];return out
small=verify(s['C'],s['rhs'].ravel(order='F'),s['smallValue'].ravel(order='F'),s['smallReport']['krylov'],P)
default_file=VAL/case/'default_replay/failure_observer_unit_summary.json';default=json.loads(default_file.read_text(encoding='utf-8-sig'))
assert default['oldCompleteKernelOutputExactlyReplayed'] and default['observerSuccessOutputExactlyReplayed'] and default['intentionalNonlinearActionRejected']
if 'nullValue' in s:
    assert np.all(s['nullValue']==0) and s['nullReport']['flag']!=0 and s['report']['zeroImageRejectedAsUnresolved']
    assert s['nullReport']['cycles']['imageNormalizationScales']==1
out=dict(status='passed-independent-scaled-Arnoldi-physical-Z-full-small-C-and-default-contract',processID=os.getpid(),toy=toy,smallC=small,
         newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
         inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,default_file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
