# -*- coding: utf-8 -*-
"""只投影保存的原NS失败剩余，比较两个已测粗像的互补性；不构造超预算联合算子。"""
from pathlib import Path
import hashlib,json,os,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/case/'alternative_coarse_projection.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'wall_velocity_krylov_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};r=loadmat(file,simplify_cells=True,variable_names=['result'])['result']
assert r['acceptedNonlinearIterations']==0;np.testing.assert_array_equal(r['initialQ'],r['Q']);Q=r['Q'];a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
rhs=r['initialRate']/r['residualScale'];remaining=(rhs+a['steadyDirectionEvidence']['value']).ravel(order='F');reference=np.linalg.norm(rhs);before=np.linalg.norm(remaining);del r,a,rhs;gc.collect()
checks=[];total_bytes=0
for name in ['spatial','primitive']:
    source=VAL/'basiscandidatesrocket01'/name/'coarse_candidate_evidence.mat';bindings[str(source.relative_to(ROOT))]=sha(source)
    coarse=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse'];np.testing.assert_array_equal(Q,coarse['baseQ'])
    image=coarse['image'];G=coarse['normalMatrix'];weights=np.linalg.solve(G,image.T@remaining);prediction=remaining-image@weights;after=np.linalg.norm(prediction)
    checks.append(dict(basis=name,coarseUnknowns=int(5*coarse['space']['groupCount']),sparseImageBytes=int(coarse['imageBytes']),
        predictedRemainingRelativeDefect=float(after/reference),remainingSquaredNormFractionRemoved=float(1-(after/before)**2)))
    total_bytes+=int(coarse['imageBytes']);del coarse,image,G,weights,prediction;gc.collect()
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='diagnosed-existing-coarse-alternative-projections-on-same-original-NS-failed-direction',processID=os.getpid(),
    actualNSRemainingRelativeDefect=float(before/reference),separateAlternativeProjections=checks,
    fullCombinedImageBytesWouldBe=total_bytes,maximumSingleImageBytes=268435456,combinedOperatorBuilt=False,
    newActualNSEvaluations=0,newModelApplications=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
