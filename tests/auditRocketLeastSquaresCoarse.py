# -*- coding: utf-8 -*-
"""独立核对细层最小二乘粗校正、原体积投影、源码绑定与非增性质。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import kron,eye,diags
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_least_squares_coarse.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'least_squares_coarse_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True);coarse=s['coarse'];space=coarse['space']
source=Path(s['sourceFile']);bindings[str(source.relative_to(ROOT))]=sha(source)
image=coarse['image'].tocsc();P=space['prolongation'].tocsr();R=space['restriction'].tocsr();ng=int(space['groupCount']);n=P.shape[0];G=np.asarray(coarse['normalMatrix'])
assert image.shape==(5*n,5*ng) and np.isfinite(image.data).all() and coarse['imageBytes']<=coarse['maximumImageBytes']
np.testing.assert_allclose(np.asarray(image.multiply(image).sum(axis=0)).ravel(),1,rtol=5e-12,atol=0)
np.testing.assert_allclose((image.T@image).toarray(),G,rtol=5e-11,atol=5e-13)
np.testing.assert_allclose(coarse['upper'].T@coarse['upper'],G,rtol=5e-11,atol=5e-13)
projected=(kron(eye(5),R)@image@diags(1/coarse['columnLiftScale'])).toarray()
projection_error=float(np.linalg.norm(projected-s['priorMatrix'])/np.linalg.norm(s['priorMatrix']));assert projection_error<1e-10
remaining=s['remaining'];vector=remaining.ravel(order='F');right=image.T@vector;coefficients=np.linalg.solve(G,right)
correction=P@(coefficients*coarse['columnLiftScale']).reshape((ng,5),order='F');expected=s['initial']+correction
update_error=float(np.linalg.norm(expected-s['update'])/max(1,np.linalg.norm(s['update'])));assert update_error<2e-8
predicted=(vector-image@coefficients).reshape((n,5),order='F');actual=s['after']
projection_action_error=float(np.linalg.norm(predicted-actual)/max(1,np.linalg.norm(actual)));assert projection_action_error<2e-8
assert np.linalg.norm(actual)<=np.linalg.norm(remaining)+1e-9*max(1,np.linalg.norm(remaining))
linearity=float(np.linalg.norm(s['combined']-(.7*s['update']-1.3*s['other']))/max(1,np.linalg.norm(s['combined'])));assert linearity<1e-9
source_bindings={}
for name,digest in coarse['sourceIdentity'].items():
    path=ROOT/'internalFlow'/('' if name=='resolveNSCharacteristicBoundary3D' else 'cutCell')/(name+'.m')
    assert sha(path)==digest,(name,digest);source_bindings[str(path.relative_to(ROOT))]=digest
details=dict(fineCells=n,coarseUnknowns=5*ng,imageNonzeros=int(image.nnz),imageBytes=int(coarse['imageBytes']),
    volumeGalerkinProjectionRelativeDifference=projection_error,independentUpdateRelativeDifference=update_error,
    predictedVersusFullModelRemainingRelativeDifference=projection_action_error,fixedLinearityRelativeDifference=linearity,
    initialFineModelRelativeResidual=float(np.linalg.norm(remaining)/np.linalg.norm(s['rhs'])),correctedFineModelRelativeResidual=float(np.linalg.norm(actual)/np.linalg.norm(s['rhs'])))
if s['report']['kind']=='unit':
    reference=np.linalg.lstsq(image.toarray(),vector,rcond=None)[0];u=s['initial']+P@(reference*coarse['columnLiftScale']).reshape((ng,5),order='F')
    qr_error=float(np.linalg.norm(u-s['update'])/max(1,np.linalg.norm(u)));assert qr_error<2e-7;details['independentDenseLeastSquaresRelativeDifference']=qr_error
    details['originalNSCoarseFDRelativeDifference']=float(np.linalg.norm(projected-s['fdMatrix'])/np.linalg.norm(s['fdMatrix']));assert details['originalNSCoarseFDRelativeDifference']<2e-6
assert s['report']['actualOriginalNSEvaluations']==0 and s['report']['newAcceptedIterations']==0
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-sparse-fine-L2-coarse-projection-nonincreasing-and-source-binding',case=case,processID=os.getpid(),details=details,
    actualOriginalNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,numericalSourceSHA256=source_bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
