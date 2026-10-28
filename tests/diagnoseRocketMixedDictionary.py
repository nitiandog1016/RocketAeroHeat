# -*- coding: utf-8 -*-
"""用小联合Gram诊断已有候选字典的投影能力，不构造或使用超预算联合细像算子。"""
from pathlib import Path
import hashlib,json,os,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/case/'mixed_dictionary_diagnostic.json';assert not target.exists()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'wall_velocity_krylov_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)};r=loadmat(file,simplify_cells=True,variable_names=['result'])['result']
assert r['acceptedNonlinearIterations']==0 and r['coarseBasis']=='mixed-conservative-primitive';np.testing.assert_array_equal(r['initialQ'],r['Q'])
Q=r['Q'];selected=np.atleast_1d(r['options']['coarseMixedColumns']).astype(int)-1;a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[-1]
rhs=r['initialRate']/r['residualScale'];remaining=(rhs+a['steadyDirectionEvidence']['value']).ravel(order='F');reference=np.linalg.norm(rhs);before=np.linalg.norm(remaining)
del r,a,rhs;gc.collect();images=[];bytes_total=0
for name in ['spatial','primitive']:
    source=VAL/'basiscandidatesrocket01'/name/'coarse_candidate_evidence.mat';bindings[str(source.relative_to(ROOT))]=sha(source)
    c=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse'];np.testing.assert_array_equal(Q,c['baseQ']);images.append(c['image']);bytes_total+=int(c['imageBytes']);del c;gc.collect()
S,P=images;ns=S.shape[1];G0=(S.T@S).toarray();GP=(P.T@P).toarray();C=(S.T@P).toarray();G=np.block([[G0,C],[C.T,GP]]);G=.5*(G+G.T);b=np.r_[S.T@remaining,P.T@remaining]
indices=np.r_[np.arange(ns),ns+selected];weights=np.linalg.solve(G[np.ix_(indices,indices)],b[indices]);current=remaining-S@weights[:ns]-P[:,selected]@weights[ns:]
eigenvalues,vectors=np.linalg.eigh(G);largest=float(eigenvalues[-1]);fits=[]
for cutoff in [1e-8,1e-10,1e-12]:
    keep=eigenvalues>cutoff*largest;v=vectors[:,keep];coeff=v@((v.T@b)/eigenvalues[keep]);prediction=remaining-S@coeff[:ns]-P@coeff[ns:];after=np.linalg.norm(prediction)
    assert after<=before+1e-8*max(1,before)
    fits.append(dict(relativeEigenvalueCutoff=cutoff,retainedRank=int(sum(keep)),predictedRemainingRelativeDefect=float(after/reference),
                     remainingSquaredNormFractionRemoved=float(1-(after/before)**2)))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='diagnosed-existing-source-dictionary-spectrum-sensitive-projections-without-joint-fine-image-allocation',processID=os.getpid(),
    actualNSRemainingRelativeDefect=float(before/reference),currentSelectedSpaceProjectedRelativeDefect=float(np.linalg.norm(current)/reference),
    sourceDictionaryColumns=G.shape[0],spectrumSensitiveFullDictionaryPredictions=fits,smallestGramEigenvalue=float(eigenvalues[0]),largestGramEigenvalue=largest,
    sourceImageBytesTotal=bytes_total,maximumSingleImageBytes=268435456,fullJointFineImageAllocated=False,fullDictionaryPreconditionerBuilt=False,
    newActualNSEvaluations=0,newModelApplications=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
