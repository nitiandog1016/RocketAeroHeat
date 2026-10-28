# -*- coding: utf-8 -*-
"""用独立稠密矩阵复算增强Arnoldi关系、最小残差、保留像与原NS方向。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
file=VAL/case/'recycled_kernel_evidence.mat';target=VAL/(case+'_independent_recycled_kernel.json');assert not target.exists()
s=loadmat(file,simplify_cells=True)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def rows(x):return [x] if isinstance(x,dict) else list(np.atleast_1d(x))
details=[]
for name,A,b,x,r in [('recycled',s['A'],s['b'],s['x'],s['recycled']),('plain',s['A'],s['b'],s['plainX'],s['plain']),
                    ('identity',np.eye(len(s['b'])),s['b'],s['identityX'],s['identity']),('unit',s['unitA'],s['unitB'],s['unitX'],s['unitReport'])]:
    n=len(b);observed_evidence=0;expected_actions=0;total_steps=0;maximum_relation=0.;maximum_projection=0.
    np.testing.assert_allclose(np.linalg.norm(b-A@x),r['absoluteResidual'],rtol=1e-9,atol=1e-10)
    assert r['basisStorageUpperBoundBytes']<=r['maximumBasisBytes'] and r['arnoldiSteps']<=r['options']['maximumArnoldiSteps']
    for c in rows(r['cycles']):
        j=int(c['arnoldiSteps']);k=int(c['retainedBefore']);next_k=int(c['retainedAfter']);total_steps+=j
        V=np.asarray(c['V']).reshape((n,j+1),order='F');Z=np.asarray(c['Z']).reshape((n,j),order='F')
        U=np.asarray(c['retainedU']).reshape((n,k),order='F');C=np.asarray(c['retainedC']).reshape((n,k),order='F')
        H=np.asarray(c['smallH']).reshape((j+1,j),order='F');weights=np.atleast_1d(c['weights'])
        relation=float(np.linalg.norm(A@Z-V@H)/max(1,np.linalg.norm(A@Z)));assert relation<1e-8;maximum_relation=max(maximum_relation,relation)
        if k:
            assert np.linalg.norm(A@U-C)/max(1,np.linalg.norm(C))<1e-8
            assert np.linalg.norm(C.T@C-np.eye(k))<1e-8
        np.testing.assert_allclose(c['candidate'],c['xStart']+Z@weights,rtol=1e-10,atol=1e-10)
        actual=b-A@c['candidate'];np.testing.assert_allclose(actual,c['actualRemaining'],rtol=1e-7,atol=1e-8)
        assert np.linalg.norm(actual)<=c['startResidualNorm']+1e-8*max(1,c['startResidualNorm'])
        if j:
            image=A@Z;reference=np.linalg.lstsq(image,c['rStart'],rcond=1e-13)[0]
            error=float(np.linalg.norm(image@(reference-weights))/max(1,np.linalg.norm(b)));assert error<2e-8;maximum_projection=max(maximum_projection,error)
        if next_k:
            nextU=np.asarray(c['nextU']).reshape((n,next_k),order='F');nextC=np.asarray(c['nextC']).reshape((n,next_k),order='F')
            assert np.linalg.norm(A@nextU-nextC)/max(1,np.linalg.norm(nextC))<1e-8
            assert np.linalg.norm(nextC.T@nextC-np.eye(next_k))<1e-8
            D=np.zeros((k+j+1,k+j));D[:k,:k]=np.eye(k);D[k:,k:]=H
            Y=np.asarray(c['retentionCoefficients']).reshape((k+j,next_k),order='F');sigma=np.atleast_1d(c['retentionSingularValues'])
            assert np.linalg.norm(Y.T@Y-np.eye(next_k))<1e-8
            assert np.linalg.norm(D.T@D@Y-Y*(sigma*sigma))/max(1,np.linalg.norm(D)**2)<1e-8
            spectrum=np.linalg.svd(D,compute_uv=False);eligible=spectrum[spectrum>r['options']['retentionRelativeThreshold']*spectrum.max()]
            np.testing.assert_allclose(np.sort(sigma),np.sort(eligible[-next_k:]),rtol=1e-8,atol=1e-12)
        expected_actions+=j+1+(next_k+1 if next_k else 0)
        for value in c.values():
            if isinstance(value,np.ndarray):observed_evidence+=value.nbytes
    assert total_steps==r['arnoldiSteps']==r['preconditionerApplications'] and expected_actions==r['modelApplications']<=r['maximumModelApplications']
    assert observed_evidence<=r['basisStorageUpperBoundBytes']
    details.append(dict(case=name,arnoldiSteps=total_steps,modelApplications=expected_actions,relativeResidual=float(r['relativeResidual']),
        maximumArnoldiRelationDifference=maximum_relation,maximumIndependentSVDProjectionDifference=maximum_projection,storedCycleEvidenceBytes=observed_evidence))
errors=[]
for direction,plus,minus,model in zip(s['validationDirections'],s['validationPlus'],s['validationMinus'],s['validationModel']):
    measured=(plus-minus)/s['scale']/(2*s['h']);error=float(np.linalg.norm(measured-model)/max(1,np.linalg.norm(measured)));assert error<2e-6;errors.append(error)
    np.testing.assert_allclose(model.ravel(order='F'),-s['unitA']@direction.ravel(order='F'),rtol=1e-12,atol=1e-7)
assert s['report']['actualOriginalNSEvaluations']==5 and s['report']['unitJacobianConstructionModelApplications']==s['unitA'].shape[1]
assert s['zero']['modelApplications']==0 and s['zero']['flag']==0
out=dict(status='passed-independent-recycled-Arnoldi-SVD-retained-images-budget-and-unit-original-NS-directions',case=case,
    processID=os.getpid(),details=details,unitOriginalNSFDErrors=errors,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
