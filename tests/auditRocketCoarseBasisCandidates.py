# -*- coding: utf-8 -*-
"""独立重建粗基提升、体积嵌套和细层投影，并保留没有收益的候选结果。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix,kron,eye,diags
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,basis=sys.argv[1:3]
assert case.replace('_','').isalnum() and basis in ('primitive','spatial');target=VAL/(case+'_'+basis+'_independent_basis.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/basis/'coarse_candidate_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True)
coarse=s['coarse'];space=coarse['space'];parent_file=Path(s['parentFile']);bindings[str(parent_file.relative_to(ROOT))]=sha(parent_file)
parent=loadmat(parent_file,simplify_cells=True,variable_names=['coarse'])['coarse'];n=len(coarse['baseQ']);ng=int(space['groupCount']);Q=coarse['baseQ'];scale=coarse['scale']
volume=space['geometry']['cells']['volume'];mapping=space['cellToGroup'].astype(int)-1;vol=np.bincount(mapping,weights=volume,minlength=ng)
P=coo_matrix((np.ones(n),(np.arange(n),mapping)),shape=(n,ng)).tocsr();R=coo_matrix((volume/vol[mapping],(mapping,np.arange(n))),shape=(ng,n)).tocsr()
def sparse_error(a,b):return float(np.max(abs((a-b).data),initial=0))
assert sparse_error(P,space['prolongation'])==0 and sparse_error(R,space['restriction'])<1e-13
np.testing.assert_allclose(vol,space['groupVolume'],rtol=5e-14,atol=0)
rho=Q[:,0];u=Q[:,1:4]/rho[:,None];cv=coarse['gas']['R']/(coarse['gas']['gamma']-1);temperature=(Q[:,4]/rho-.5*np.sum(u*u,axis=1))/cv
def lift(coeff):
    fine=P@coeff
    if basis=='spatial':return fine
    drho=fine[:,0];du=fine[:,1:4];dT=fine[:,4]
    return np.column_stack((drho,drho[:,None]*u+rho[:,None]*du,
        (cv*temperature+.5*np.sum(u*u,axis=1))*drho+rho*(np.sum(u*du,axis=1)+cv*dT)))/scale
if basis=='primitive':assert sha(ROOT/'internalFlow/cutCell/liftRocketPrimitiveCoarse.m')==coarse['primitiveLiftSourceSHA256']
image=coarse['image'].tocsc();G=(image.T@image).toarray();np.testing.assert_allclose(G,coarse['normalMatrix'],rtol=5e-11,atol=5e-13)
np.testing.assert_allclose(coarse['upper'].T@coarse['upper'],G,rtol=5e-11,atol=5e-13);assert coarse['imageBytes']<=268435456
remaining=s['remaining'].ravel(order='F');coeff=np.linalg.solve(G,image.T@remaining)
expected=lift((coeff*coarse['columnLiftScale']).reshape((ng,5),order='F'));correction_error=float(np.linalg.norm(expected-s['correction'])/max(1,np.linalg.norm(expected)));assert correction_error<3e-8
prediction=remaining-image@coeff;prediction_error=float(np.linalg.norm(prediction-s['predicted'].ravel(order='F'))/max(1,np.linalg.norm(prediction)));assert prediction_error<3e-8
assert np.linalg.norm(prediction)<=np.linalg.norm(remaining)+1e-8*max(1,np.linalg.norm(remaining))
linearity=float(np.linalg.norm(s['combined']-(.7*s['correction']-1.3*s['other']))/max(1,np.linalg.norm(s['combined'])));assert linearity<1e-8
# 单位首稿的coefficients/direction变量保留最后FD方向；像重放使用源中明确的固定0.17探针，不能混用。
probe=np.sin(np.arange(1,5*ng+1)*.17).reshape((ng,5),order='F');probe_product=image@(probe.ravel(order='F')/coarse['columnLiftScale'])
np.testing.assert_allclose(probe_product,s['imageProduct'],rtol=3e-10,atol=1e-7)
assert np.linalg.norm(s['actual'].ravel(order='F')-probe_product)/max(1,np.linalg.norm(probe_product))<1e-9
details=dict(coarseUnknowns=5*ng,imageBytes=int(coarse['imageBytes']),independentCorrectionRelativeDifference=correction_error,
    independentPredictionRelativeDifference=prediction_error,fixedLinearityRelativeDifference=linearity,
    predictedRelativeDefect=float(np.linalg.norm(prediction)/np.linalg.norm(s['rhs'])),extraSquaredNormFractionRemoved=float(1-np.linalg.norm(prediction)**2/np.linalg.norm(s['oldPredicted'])**2))
if basis=='spatial':
    transfer=R@parent['space']['prolongation'];assert sparse_error(P@transfer,parent['space']['prolongation'])<1e-12
    refinement=space['targetRefinement'];cut=space['geometry'];V=cut['inputGeometry']['V'];F=cut['inputGeometry']['F'].astype(int)-1;nf=len(F);nl=n//nf
    fraction=(V[F].mean(axis=1)[:,0]-V[:,0].min())/(V[:,0].max()-V[:,0].min());layer=np.arange(n)//nf+1
    targeted=np.tile(fraction>=refinement['options']['axialFractionThreshold'],nl)&np.isin(layer,np.atleast_1d(refinement['options']['targetLayers']))
    old_map=parent['space']['cellToGroup'].astype(int)-1;bucket=np.zeros(n,dtype=int);c=cut['cells']['centroid']
    for group in range(int(parent['space']['groupCount'])):
        ids=np.flatnonzero((old_map==group)&targeted)
        if len(ids):bucket[ids]=1+(c[ids,0]>np.median(c[ids,0]))
    faces=cut['openFaces'];inside=(faces['area']>0)&(faces['cellMinus']>0)&(faces['cellPlus']>0);m=faces['cellMinus'][inside].astype(int)-1;p=faces['cellPlus'][inside].astype(int)-1
    same=(old_map[m]==old_map[p])&(bucket[m]==bucket[p]);a=coo_matrix((np.ones(2*sum(same)),(np.r_[m[same],p[same]],np.r_[p[same],m[same]])),shape=(n,n)).tocsr()
    count,labels=connected_components(a,directed=False);assert count==ng
    for group in range(ng):assert len(np.unique(labels[mapping==group]))==1
    details['independentConnectedGroupCount']=int(count)
    parent_of=np.asarray(transfer.argmax(axis=1)).ravel();counts=np.bincount(parent_of,minlength=int(parent['space']['groupCount']));copied=0
    for b in range(5):
        for g,pg in enumerate(parent_of):
            if counts[pg]!=1:continue
            col=b*ng+g;old_col=b*int(parent['space']['groupCount'])+pg
            assert sparse_error(image[:,col],parent['image'][:,old_col])==0 and coarse['columnLiftScale'][col]==parent['columnLiftScale'][old_col];copied+=1
    assert copied==coarse['copiedParentImageColumns'] and coarse['constructionModelApplications']==5*ng-copied
    details['unchangedParentImageColumnsExactlyReused']=copied
if len(s['validationDirections']):
    errors=[]
    for j,(direction,plus,minus,model) in enumerate(zip(s['validationDirections'],s['validationPlus'],s['validationMinus'],s['validationModel']),1):
        coefficients=np.sin(np.arange(1,5*ng+1)*(.13+.09*j)).reshape((ng,5),order='F');np.testing.assert_allclose(lift(coefficients),direction,rtol=5e-13,atol=1e-14)
        measured=(plus-minus)/scale/(2*s['h']);error=float(np.linalg.norm(measured-model)/max(1,np.linalg.norm(measured)));assert error<2e-6;errors.append(error)
    assert len(errors)==2 and s['result']['validationOriginalNSCalls']==4;details['originalNSDirectionFDRelativeDifferences']=errors
    reference=np.linalg.lstsq(image.toarray(),remaining,rcond=None)[0];assert np.linalg.norm(image@(reference-coeff))/max(1,np.linalg.norm(prediction))<3e-8
    details['legacyUnitProbeVariableNamesClarified']=not bool('imageProbeDirection' in s)
else:
    assert s['result']['validationOriginalNSCalls']==0
    source=Path(s['sourceFile']);bindings[str(source.relative_to(ROOT))]=sha(source);r=loadmat(source,simplify_cells=True,variable_names=['result'])['result']
    a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
    np.testing.assert_array_equal(s['remaining'],r['initialRate']/r['residualScale']+a['steadyDirectionEvidence']['value']);np.testing.assert_array_equal(Q,r['Q'])
for name,digest in coarse['sourceIdentity'].items():
    path=ROOT/'internalFlow'/('' if name=='resolveNSCharacteristicBoundary3D' else 'cutCell')/(name+'.m');assert sha(path)==digest
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-coarse-candidate-lift-projection-transfer-and-direction-evidence',case=case,basis=basis,processID=os.getpid(),details=details,
    actualNonlinearUpdates=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
