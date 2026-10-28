# -*- coding: utf-8 -*-
"""独立核对逐层体积嵌套、稀疏像空间与投影收益，不把模型预测当实际NS更新。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix,kron,eye,diags
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_layer_coarse.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'coarse_layer_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True)
parent=Path(s['parentFile']);bindings[str(parent.relative_to(ROOT))]=sha(parent);old=loadmat(parent,simplify_cells=True,variable_names=['coarse'])['coarse']
coarse=s['coarse'];space=coarse['space'];v=space['geometry']['cells']['volume'];mapping=space['cellToGroup'].astype(int)-1;n=len(v);ng=int(space['groupCount'])
vol=np.bincount(mapping,weights=v,minlength=ng);P=coo_matrix((np.ones(n),(np.arange(n),mapping)),shape=(n,ng)).tocsr()
R=coo_matrix((v/vol[mapping],(mapping,np.arange(n))),shape=(ng,n)).tocsr()
def error(a,b):return float(np.max(abs((a-b).data),initial=0))
assert error(P,space['prolongation'])==0 and error(R,space['restriction'])<1e-13
np.testing.assert_allclose(vol,space['groupVolume'],rtol=5e-14,atol=0)
parent_transfer=R@old['space']['prolongation'];assert error(parent_transfer,s['parentTransfer'])<1e-13
assert error(P@parent_transfer,old['space']['prolongation'])<1e-12
image=coarse['image'].tocsc();old_image=old['image'].tocsc();T=kron(eye(5),parent_transfer).tocsc()
T=diags(1/coarse['columnLiftScale'])@T@diags(old['columnLiftScale'])
parent_error=float(np.linalg.norm((image@T-old_image).data)/np.linalg.norm(old_image.data));assert parent_error<1e-9
G=(image.T@image).toarray();np.testing.assert_allclose(G,coarse['normalMatrix'],rtol=5e-11,atol=5e-13)
np.testing.assert_allclose(coarse['upper'].T@coarse['upper'],G,rtol=5e-11,atol=5e-13)
assert coarse['imageBytes']<=coarse['maximumImageBytes'] and np.isfinite(image.data).all()
remaining=s['failureRemaining'].ravel(order='F');rhs_norm=np.linalg.norm(s['rhs']);projections={}
for name,c in [('old',old),('new',coarse)]:
    basis=c['image'];a=np.linalg.solve(c['normalMatrix'],basis.T@remaining);prediction=remaining-basis@a
    expected=c['space']['prolongation']@(a*c['columnLiftScale']).reshape((int(c['space']['groupCount']),5),order='F')
    np.testing.assert_allclose(s[name+'Correction'],expected,rtol=3e-8,atol=1e-9)
    actual=s[name+'Predicted'].ravel(order='F');projection_error=float(np.linalg.norm(actual-prediction)/max(1,np.linalg.norm(prediction)));assert projection_error<3e-8
    projections[name]=dict(norm=float(np.linalg.norm(prediction)),relativeToOriginalRHS=float(np.linalg.norm(prediction)/rhs_norm),independentPredictionRelativeDifference=projection_error)
assert projections['new']['norm']<=projections['old']['norm']+1e-8*max(1,projections['old']['norm'])
fraction=1-(projections['new']['norm']/projections['old']['norm'])**2
assert abs(fraction-s['report']['additionalRemainingSquaredNormFractionRemoved'])<2e-8
details=dict(fineCells=n,coarseUnknowns=5*ng,parentImageRelativeDifference=parent_error,projection=projections,additionalRemainingSquaredNormFractionRemoved=fraction)
if s['report']['kind']=='unit':
    reference=np.linalg.lstsq(image.toarray(),remaining,rcond=None)[0];prediction=remaining-image@reference
    assert np.linalg.norm(prediction-s['newPredicted'].ravel(order='F'))/max(1,np.linalg.norm(prediction))<3e-8
    directions=list(s['validationDirections']);plus=list(s['validationPlus']);minus=list(s['validationMinus']);models=list(s['validationModel']);errors=[]
    for direction,p,m,model in zip(directions,plus,minus,models):
        measured=(p-m)/coarse['scale']/(2*s['h']);err=float(np.linalg.norm(measured-model)/max(1,np.linalg.norm(measured)));assert err<2e-6;errors.append(err)
    assert len(errors)==2 and s['report']['actualOriginalNSEvaluations']==5;details['actualNSCoarseDirectionFDRelativeDifferences']=errors
else:
    source=Path(s['sourceFile']);bindings[str(source.relative_to(ROOT))]=sha(source)
    r=loadmat(source,simplify_cells=True,variable_names=['result'])['result'];a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
    np.testing.assert_array_equal(s['failureRemaining'],r['initialRate']/r['residualScale']+a['steadyDirectionEvidence']['value'])
    np.testing.assert_array_equal(coarse['baseQ'],r['Q']);assert s['report']['actualOriginalNSEvaluations']==0
for name,digest in coarse['sourceIdentity'].items():
    path=ROOT/'internalFlow'/('' if name=='resolveNSCharacteristicBoundary3D' else 'cutCell')/(name+'.m');assert sha(path)==digest
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-nested-layer-coarse-transfer-image-and-projection-evidence',case=case,processID=os.getpid(),details=details,
    actualNonlinearUpdates=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
