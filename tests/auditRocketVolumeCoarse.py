# -*- coding: utf-8 -*-
"""独立检查正面邻接聚合、真实体积传递、粗层LU及模型投影；不另写NS求解器。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix,eye
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,kind=sys.argv[1:3]
assert case.replace('_','').isalnum() and kind in ('unit','space','real');target=VAL/(case+'_independent_volume_coarse.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
filename={'unit':'volume_coarse_unit_evidence.mat','space':'volume_coarse_space_evidence.mat','real':'volume_coarse_real_evidence.mat'}[kind]
file=VAL/case/filename;bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True)
if kind=='space':
    source=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat';bindings[str(source.relative_to(ROOT))]=sha(source)
    cut=loadmat(source,simplify_cells=True,variable_names=['domain'])['domain'];map=s['map'].astype(int)-1;P=s['P'].tocsr();R=s['R'].tocsr();vol=s['groupVolume'];options=s['report']['options']
else:
    space=s['space'];cut=space['geometry'];map=space['cellToGroup'].astype(int)-1;P=space['prolongation'].tocsr();R=space['restriction'].tocsr();vol=space['groupVolume'];options=space['options']
v=cut['cells']['volume'];c=cut['cells']['centroid'];n=len(v);ng=len(vol);assert map.min()==0 and map.max()+1==ng and np.all(v>0)
expected_volume=np.bincount(map,weights=v,minlength=ng);np.testing.assert_allclose(expected_volume,vol,rtol=3e-14,atol=0)
expected_P=coo_matrix((np.ones(n),(np.arange(n),map)),shape=(n,ng)).tocsr();expected_R=coo_matrix((v/expected_volume[map],(map,np.arange(n))),shape=(ng,n)).tocsr()
def sparse_error(a,b):
    difference=a-b;return float(np.max(abs(difference.data),initial=0))
assert sparse_error(P,expected_P)==0 and sparse_error(R,expected_R)<3e-14
identity_error=sparse_error(R@P,eye(ng));assert identity_error<1e-12
bins=np.zeros((n,3),dtype=int)
for d,count in enumerate(np.asarray(options['cartesianBins']).astype(int)):
    indices=np.maximum(1,np.minimum(n,np.floor(np.arange(1,count)*n/count+.5).astype(int)))-1
    thresholds=np.unique(np.sort(c[:,d])[indices]);bins[:,d]=np.sum(c[:,d,None]>thresholds[None,:],axis=1)
nf=len(cut['inputGeometry']['F']);layer=np.arange(n)//nf+1;breaks=np.atleast_1d(options['layerBreaks']);radial=np.sum(layer[:,None]>breaks[None,:],axis=1)
_,bucket=np.unique(np.column_stack((bins,radial)),axis=0,return_inverse=True)
f=cut['openFaces'];inside=(f['area']>0)&(f['cellMinus']>0)&(f['cellPlus']>0);m=f['cellMinus'][inside].astype(int)-1;p=f['cellPlus'][inside].astype(int)-1;same=bucket[m]==bucket[p]
assert np.all(map[m[same]]==map[p[same]])
adj=coo_matrix((np.ones(2*sum(same)),(np.r_[m[same],p[same]],np.r_[p[same],m[same]])),shape=(n,n)).tocsr()
count,labels=connected_components(adj,directed=False);assert count==ng
for group in range(ng):
    assert len(np.unique(bucket[map==group]))==1 and len(np.unique(labels[map==group]))==1
details=dict(fineCells=n,coarseGroups=ng,coarseUnknowns=5*ng,connectedComponents=int(count),restrictionProlongationIdentityError=identity_error,
    totalVolumeRelativeDifference=float(abs(sum(expected_volume)-sum(v))/sum(v)),fineCellDiscardCount=0)
if kind!='space':
    coarse=s['coarse'];C=coarse['matrix'];permutation=coarse['rowPermutation'].astype(int)-1
    factor=np.linalg.norm(coarse['lower']@coarse['upper']-C[permutation])/max(1,np.linalg.norm(C));assert factor<1e-12
    details['coarseFactorRelativeDifference']=float(factor)
    if kind=='unit':
        fd_error=float(np.linalg.norm(C-s['fdMatrix'])/np.linalg.norm(s['fdMatrix']));assert fd_error<2e-6
        solution=np.linalg.solve(C,s['coarseRhs'].ravel(order='F')).reshape((ng,5),order='F')
        np.testing.assert_allclose(solution,s['coarseSolution'],rtol=2e-10,atol=2e-12)
        linear=.7*s['update']-1.3*s['other'];linearity=float(np.linalg.norm(s['combined']-linear)/max(1,np.linalg.norm(linear)));assert linearity<1e-11
        details.update(measuredOriginalNSCoarseFDRelativeDifference=fd_error,fixedTwoLevelLinearityRelativeDifference=linearity,
            actualOriginalNSEvaluations=int(s['report']['actualOriginalNSEvaluations']))
    else:
        projected=(C@s['coarseVector'].ravel(order='F')).reshape((ng,5),order='F')
        np.testing.assert_allclose(projected,s['modelCoarse'],rtol=2e-10,atol=2e-8)
        np.testing.assert_allclose(R@s['afterCorrection'],s['restricted'],rtol=2e-10,atol=2e-8)
        residual=float(np.linalg.norm(R@s['afterCorrection'])/max(1,np.linalg.norm(R@s['rhs'])));assert residual<1e-8
        details.update(correctedModelCoarseRelativeResidual=residual,
            initialFineModelRelativeResidual=float(np.linalg.norm(s['beforeCorrection'])/np.linalg.norm(s['rhs'])),
            correctedFineModelRelativeResidual=float(np.linalg.norm(s['afterCorrection'])/np.linalg.norm(s['rhs'])),actualOriginalNSEvaluations=0)
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-connected-volume-transfer-and-Galerkin-coarse-'+kind,case=case,processID=os.getpid(),details=details,
    fullIndependentNSSolver=False,physicalMeshChanged=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
