# -*- coding: utf-8 -*-
"""独立核对有界混合粗像逐列来源、提升、Gram/SVD投影及原NS方向。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_mixed_coarse.json');assert not target.exists()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'mixed_coarse_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True)
coarse=s['coarse'];selected=np.atleast_1d(coarse['selectedPrimitiveColumns']).astype(int)-1;image=coarse['image'].tocsc();nc=int(coarse['baseConservativeColumns']);n=len(coarse['baseQ'])
assert len(set(selected))==len(selected) and image.shape==(5*n,nc+len(selected)) and len(selected)<=coarse['selectionOptions']['maximumAdditionalColumns']
assert int(coarse['imageBytes'])==16*image.nnz+8*(image.shape[1]+1)<=coarse['maximumImageBytes']==268435456
for key,part in [('spatialFile','spatial'),('primitiveFile','primitive')]:
    p=Path(s[key]);bindings[str(p.relative_to(ROOT))]=sha(p);parent=loadmat(p,simplify_cells=True,variable_names=['coarse'])['coarse']
    if part=='spatial':
        diff=image[:,:nc]-parent['image'];assert not diff.nnz or np.max(abs(diff.data))==0
        np.testing.assert_array_equal(coarse['columnLiftScale'][:nc],parent['columnLiftScale']);np.testing.assert_array_equal(coarse['baseQ'],parent['baseQ'])
    else:
        diff=image[:,nc:]-parent['image'][:,selected];assert not diff.nnz or np.max(abs(diff.data))==0
        np.testing.assert_array_equal(coarse['columnLiftScale'][nc:],parent['columnLiftScale'][selected]);np.testing.assert_array_equal(coarse['baseQ'],parent['baseQ'])
G=(image.T@image).toarray();np.testing.assert_allclose(G,coarse['normalMatrix'],rtol=5e-11,atol=5e-13)
np.testing.assert_allclose(coarse['upper'].T@coarse['upper'],G,rtol=5e-11,atol=5e-13)
P=coarse['space']['prolongation'];PP=coarse['primitiveSpace']['prolongation'];ng=P.shape[1];npg=PP.shape[1];Q=coarse['baseQ'];rho=Q[:,0];u=Q[:,1:4]/rho[:,None]
cv=coarse['gas']['R']/(coarse['gas']['gamma']-1);T=(Q[:,4]/rho-.5*np.sum(u*u,axis=1))/cv
def lift(c):
    c=c*coarse['columnLiftScale'];base=P@c[:nc].reshape((ng,5),order='F');primitive=np.zeros(5*npg);primitive[selected]=c[nc:]
    fine=PP@primitive.reshape((npg,5),order='F');dr=fine[:,0];du=fine[:,1:4];dt=fine[:,4]
    dq=np.column_stack((dr,u*dr[:,None]+rho[:,None]*du,(cv*T+.5*np.sum(u*u,axis=1))*dr+rho*(np.sum(u*du,axis=1)+cv*dt)))
    return base+dq/coarse['scale']
remaining=s['remaining'].ravel(order='F');c=np.linalg.solve(G,image.T@remaining);prediction=remaining-image@c
correction_error=float(np.linalg.norm(lift(c)-s['correction'])/max(1,np.linalg.norm(s['correction'])));assert correction_error<5e-8
prediction_error=float(np.linalg.norm(prediction-s['predicted'].ravel(order='F'))/max(1,np.linalg.norm(remaining)));assert prediction_error<3e-8
np.testing.assert_allclose(lift(s['probeCoefficients']),s['probeDirection'],rtol=5e-12,atol=1e-12)
assert np.linalg.norm(image@s['probeCoefficients']-s['probeActual'].ravel(order='F'))/max(1,np.linalg.norm(s['probeActual']))<1e-9
linearity=float(np.linalg.norm(s['combined']-(.7*s['correction']-1.3*s['other']))/max(1,np.linalg.norm(s['combined'])));assert linearity<1e-8
assert np.linalg.norm(prediction)<=np.linalg.norm(s['oldPredicted'])+1e-8*max(1,np.linalg.norm(s['oldPredicted']))
for name,digest in coarse['liftSourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest
details=dict(baseColumns=nc,selectedPrimitiveColumns=(selected+1).tolist(),mixedColumns=image.shape[1],imageBytes=int(coarse['imageBytes']),
    independentLiftCorrectionRelativeDifference=correction_error,independentPredictionRelativeDifference=prediction_error,fixedLinearityRelativeDifference=linearity,
    predictedRemainingRelativeDefect=float(np.linalg.norm(prediction)/np.linalg.norm(s['rhs'])))
if s['report']['kind']=='unit':
    ref=np.linalg.lstsq(image.toarray(),remaining,rcond=None)[0];svd_error=float(np.linalg.norm(image@(ref-c))/max(1,np.linalg.norm(remaining)));assert svd_error<3e-8;details['independentDenseSVDImageDifference']=svd_error
    errors=[]
    for plus,minus,model in zip(s['validationPlus'],s['validationMinus'],s['validationModel']):
        measured=(plus-minus)/coarse['scale']/(2*s['h']);err=float(np.linalg.norm(measured-model)/max(1,np.linalg.norm(measured)));assert err<2e-6;errors.append(err)
    assert len(errors)==2 and s['report']['actualOriginalNSEvaluations']==5;details['actualNSDirectionFDRelativeDifferences']=errors
else:
    p=Path(s['sourceFile']);bindings[str(p.relative_to(ROOT))]=sha(p);r=loadmat(p,simplify_cells=True,variable_names=['result'])['result'];a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
    np.testing.assert_array_equal(s['remaining'],r['initialRate']/r['residualScale']+a['steadyDirectionEvidence']['value']);np.testing.assert_array_equal(Q,r['Q'])
    assert s['report']['actualOriginalNSEvaluations']==0
assert s['report']['constructionModelApplications']==0 and all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-bounded-mixed-column-identity-lift-projection-and-direction-evidence',case=case,processID=os.getpid(),details=details,
    actualFlowUpdates=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
