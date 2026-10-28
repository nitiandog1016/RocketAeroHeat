# -*- coding: utf-8 -*-
"""小域稠密Schur参考，区分热力学变量时间正则；不把参考解冒充真实求解器。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1]
assert case.replace('_','').isalnum();target=VAL/case/'dense_schur_references.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'inner_graph_unit_evidence.mat';s=loadmat(file,simplify_cells=True);source=Path(s['sourceFile'])
t=loadmat(source,simplify_cells=True,variable_names=['operator']);tau=float(t['operator']['pseudoTimeStepSeconds'])
A,B=s['steadyA'],s['B'];n=len(s['Q']);N=len(A);b=(s['rate']/s['scale']).ravel(order='F');eye=np.eye(N);results=[]
for name,fields in [('allShifted',[]),('densityUnshifted',[0]),('energyUnshifted',[4]),('densityEnergyUnshifted',[0,4]),('allUnshifted',[0,1,2,3,4])]:
    selected=np.concatenate([np.arange(n)+field*n for field in fields]).astype(int) if fields else np.array([],int)
    C=B.copy();C[selected,selected]-=1/tau;sv=np.linalg.svd(C,compute_uv=False);threshold=N*np.finfo(float).eps*sv[0];rank=int(np.count_nonzero(sv>threshold))
    row=dict(name=name,unshiftedConservativeFields=fields,numericalRank=rank,rankThreshold=float(threshold),
        minimumSingularValue=float(sv[-1]),conditionNumber=float(sv[0]/sv[-1]) if sv[-1]>0 else None,
        notImplementedForRealDomain=True,physicalSolution=False)
    if rank<N:
        row['status']='rank-deficient-dense-reference-not-a-usable-inverse';results.append(row);continue
    inverse=np.linalg.solve(C,eye);x=inverse@b;value=A@x;alpha=float(b@value/(value@value));preconditioned=np.linalg.svd(A@inverse,compute_uv=False)
    row.update(status='finite-dense-reference-only',steadyAPreconditionedCondition=float(preconditioned[0]/preconditioned[-1]),
        currentAResidualOfReferenceDirection=float(np.linalg.norm(b-value)/np.linalg.norm(b)),
        bestScalar=float(alpha),bestScalarModelRelativeDefect=float(np.linalg.norm(b-alpha*value)/np.linalg.norm(b)))
    if 0<len(selected)<N:
        other=np.setdiff1d(np.arange(N),selected);D=C[np.ix_(other,other)];F=C[np.ix_(selected,other)];E=C[np.ix_(other,selected)]
        DiE=np.linalg.solve(D,E);S=C[np.ix_(selected,selected)]-F@DiE;yd=np.linalg.solve(D,b[other]);xs=np.linalg.solve(S,b[selected]-F@yd)
        schur=np.zeros(N);schur[selected]=xs;schur[other]=yd-DiE@xs;difference=float(np.linalg.norm(schur-x)/max(1,np.linalg.norm(x)));assert difference<2e-8
        singular=np.linalg.svd(S,compute_uv=False);row.update(schurScalarUnknowns=len(selected),eliminatedUnknowns=len(other),
            schurCondition=float(singular[0]/singular[-1]),denseSchurReferenceDifference=difference)
    results.append(row)
out=dict(status='source-bound-small-domain-dense-Schur-time-regularization-reference-only',processID=os.getpid(),cellCount=n,conservativeUnknowns=N,
    pseudoTimeStepSeconds=tau,candidates=results,notARealDomainOrPhysicalValidation=True,preconditionerChangesPhysicalEquation=False,
    newSelectedModelApplications=0,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,source)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
