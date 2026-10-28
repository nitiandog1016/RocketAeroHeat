# -*- coding: utf-8 -*-
"""独立核对冻结完整C的全部Arnoldi列、实际剩余及首末两列物理预条件方向。"""
from pathlib import Path
import hashlib,json,os,sys,time
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_scaled_ilu_replay.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)
folder=VAL/case;file=folder/'scaled_ilu_replay_evidence.mat';rejected=False
if file.exists():
    s=read(file);k=s['detail']['krylov'];c=k['cycles'];b=s['rhs'].ravel(order='F');actual=c['actualRemaining'];candidate=c['candidate'];binding=s
else:
    file=folder/'scaled_arnoldi_failure_evidence.mat';s=read(file);c=s['failure'];b=c['b'];actual=c['actualRemaining'];candidate=c['candidate'];binding=read(folder/'input_binding.mat');k=c;rejected=True
j=int(c['arnoldiSteps']);V=c['V'];Z=c['Z'];H=c['smallH'];scales=np.atleast_1d(c['imageNormalizationScales'])
pack=read(binding['packedFile'],['packed'])['packed'];n=int(pack['cellCount']);order=pack['order'].astype(int)-1
ptr=pack['rowPointer'].astype(np.int64)-1;col=pack['column'].astype(np.int64)-1;pv=pack['values'].transpose(2,0,1);owners=np.repeat(np.arange(n),np.diff(ptr))
def action(vector):
    v=vector.reshape((n,5),order='F')[order];out=np.zeros((n,5))
    for first in range(0,len(col),4096):
        ids=slice(first,first+4096);np.add.at(out,owners[ids],np.einsum('kab,kb->ka',pv[ids],v[col[ids]]))
    physical=np.zeros((n,5));physical[order]=out;return physical.ravel(order='F')
maximum_relation=0.
for i in range(j):
    image=action(Z[:,i]);formed=V@H[:,i];error=np.linalg.norm(image-formed)/max(np.finfo(float).tiny,np.linalg.norm(image))
    maximum_relation=max(maximum_relation,float(error));assert error<1e-8
    if (i+1)%8==0:print('INDEPENDENT_SCALED_C columns=%d/%d'%(i+1,j),flush=True)
remaining=b-action(candidate);rhs_norm=np.linalg.norm(b);actual_difference=float(np.linalg.norm(remaining-actual)/rhs_norm);assert actual_difference<1e-8
residual=float(np.linalg.norm(remaining)/rhs_norm)
del pack,pv,owners
f=read(binding['factorFile'],['factor'])['factor'];fv=f['values'].transpose(2,0,1);diag=f['diagonalPosition'].astype(int)-1;inverse=f['inverseDiagonal'].transpose(2,0,1)
np.testing.assert_array_equal(f['order'].astype(int)-1,order)
fptr=f['rowPointer'].astype(np.int64)-1;fcol=f['column'].astype(np.int64)-1
probe_errors=[]
for probe in sorted(set([0,j-1])):
    rhs=V[:,probe].reshape((n,5),order='F')[order];forward=np.zeros((n,5));solution=np.zeros((n,5))
    for i in range(n):
        ids=slice(fptr[i],diag[i]);forward[i]=rhs[i]-np.einsum('kab,kb->a',fv[ids],forward[fcol[ids]])
    for i in range(n-1,-1,-1):
        ids=slice(diag[i]+1,fptr[i+1]);solution[i]=inverse[i]@(forward[i]-np.einsum('kab,kb->a',fv[ids],solution[fcol[ids]]))
    physical=np.zeros((n,5));physical[order]=solution/scales[probe];reference=physical.ravel(order='F')
    error=float(np.linalg.norm(reference-Z[:,probe])/max(np.finfo(float).tiny,np.linalg.norm(Z[:,probe])));assert error<1e-8
    probe_errors.append(dict(column=probe+1,relativePhysicalDirectionDifference=error));print('INDEPENDENT_SCALED_P column=%d/%d'%(probe+1,j),flush=True)
active=V.shape[1]-int(np.linalg.norm(V[:,-1])==0);activeV=V[:,:active]
orth=float(np.linalg.norm(activeV.T@activeV-np.eye(active)));assert orth<1e-8
assert k['modelApplications']==j+1 and k['preconditionerApplications']==j
if 'originalNSAcceptanceClaimed' in k:assert not k['originalNSAcceptanceClaimed']
assert max(c['reorthogonalizationPasses'])<=8
if not rejected:
    assert abs(residual-k['relativeResidual'])<1e-8 and not s['report']['physicalSolution'] and s['report']['newActualNSEvaluations']==0
else:assert c['orthogonalityThreshold']==c['modelRelationThreshold']==1e-8 and not c['usableDirectionClaimed']
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-frozen-real-C-all-Arnoldi-columns-actual-residual-and-two-ILU-direction-probes',processID=os.getpid(),
    originalArnoldiGateRejected=rejected,cellCount=n,unknowns=5*n,arnoldiSteps=j,relativeCResidual=residual,
    basisOrthogonalityError=orth,maximumCompleteCColumnRelationDifference=maximum_relation,actualResidualRecomputationRelativeToRhs=actual_difference,
    physicalILUDirectionProbeErrors=probe_errors,imageScaleRange=[float(min(scales)),float(max(scales))],
    maximumReorthogonalizationPasses=int(max(c['reorthogonalizationPasses'])),CApplications=int(k['modelApplications']),ILUApplications=int(k['preconditionerApplications']),
    newFactorConstructions=0,newCoarseConstructions=0,newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
