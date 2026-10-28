# -*- coding: utf-8 -*-
"""比较同Q两条保存原NS有限差分方向的张成空间；预测不替代新组合原NS复核。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,previous=sys.argv[1:3]
assert case!=previous and all(x.replace('_','').isalnum() for x in (case,previous))
target=VAL/case/'saved_direction_span_diagnostic.json';assert not target.exists()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

files=[VAL/x/'wall_velocity_krylov_evidence.mat' for x in (previous,case)]
bindings={str(p.relative_to(ROOT)):sha(p) for p in files};base=None;images=[];directions=[];steps=[];defects=[]
for path in files:
    r=loadmat(path,simplify_cells=True,variable_names=['result'])['result']
    logs=r['linearLog'];log=logs if isinstance(logs,dict) else list(logs)[-1];e=log['steadyDirectionEvidence']
    assert e['baseAcceptedIteration']==0 and r['acceptedNonlinearIterations']==0
    if base is None:base=r['initialQ'].copy();rhs=r['initialRate']/r['residualScale'];scale=r['residualScale'].copy()
    else:
        np.testing.assert_array_equal(base,r['initialQ']);np.testing.assert_array_equal(rhs,r['initialRate']/r['residualScale']);np.testing.assert_array_equal(scale,r['residualScale'])
    images.append(-e['value'].ravel(order='F'));directions.append(e['directionU']);steps.append(float(e['actualScalarStep']))
    defects.append(float(log['fullSystemRelativeDefect']));del r
T=np.column_stack(images);vector=rhs.ravel(order='F');weights,_,rank,singular=np.linalg.lstsq(T,vector,rcond=None)
prediction=vector-T@weights;combined=weights[0]*directions[0]+weights[1]*directions[1]
defect=float(np.linalg.norm(prediction)/np.linalg.norm(vector))
assert defect<=min(defects)+1e-10
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='source-bound-two-saved-finite-difference-direction-span-diagnostic-not-new-original-NS-validation',
    processID=os.getpid(),cases=[previous,case],savedActualNSRelativeDefects=defects,finiteDifferenceSteps=steps,
    spanWeights=weights.tolist(),spanRank=int(rank),spanCondition=float(singular[0]/singular[-1]),
    predictedSpanRelativeDefect=defect,originalActualLinearTarget=.1,predictionReachesTarget=bool(defect<=.1),
    combinedScaledDirectionInfinityNorm=float(np.max(abs(combined))),
    finiteDifferenceImagesAreApproximate=True,newCombinedDirectionCheckedWithOriginalNS=False,
    newModelApplications=0,newActualNSEvaluations=0,newFlowUpdates=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
