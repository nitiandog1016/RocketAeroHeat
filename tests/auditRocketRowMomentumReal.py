# -*- coding: utf-8 -*-
"""独立核对冻结因子复用、全部行尺度、PR^{-1}及每条R*A物理方向。"""
from pathlib import Path
import os,json,hashlib,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum();folder=VAL/name;target=VAL/(name+'_independent_row_momentum.json');assert not target.exists();bindings={}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for b in iter(lambda:stream.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path,variables=None):
 path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path);return loadmat(path,simplify_cells=True,variable_names=variables)
def norm(a):
 scale=float(np.max(abs(a))) if np.size(a) else 0.
 return float(np.linalg.norm(a/scale)*scale) if scale else 0.
def rel(a,b):return norm(a-b)/max(1,norm(b))
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
summary_file=folder/'row_momentum_real_summary.json';assert summary_file.stat().st_size<128*1024;summary=json.loads(summary_file.read_text(encoding='utf-8-sig'));assert summary['numericalRunComplete'] and not summary['candidateRejected']
s=read(folder/'row_momentum_real_evidence.mat');binding=s['sourceBinding'];assert summary['sourceBinding']==binding
for key,value in s['sourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==value
previous=json.loads((VAL/'momentumlinereal03_independent_momentum_lines.json').read_text());assert previous['status'].startswith('passed-')
fdata=read(s['factorFile']);f=fdata['factor'];factor_key=str(Path(s['factorFile']).relative_to(ROOT));assert bindings[factor_key]==s['factorSHA256']==previous['inputSHA256'][factor_key]
a=read(s['assemblyFile'],['momentumMatrix']);m=a['momentumMatrix'];n=int(m['cellCount']);nf=int(f['sourceFaces']);nl=int(f['layers']);N=3*nl
ptr=m['rowPointer'].astype(np.int64)-1;col=m['column'].astype(np.int64)-1;owners=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr));values=m['values'].transpose(2,0,1);mapping=f['inverseOrder'].astype(np.int64)-1;inverse=f['inverseLine'].transpose(2,0,1)
row_max=np.zeros((n,3))
for component in range(3):np.maximum.at(row_max[:,component],owners,np.max(abs(values[:,component,:]),axis=1))
R=1/row_max;np.testing.assert_array_equal(R,s['rowScale']);assert np.all(np.isfinite(R)&(R>0))
def action(vector):
 out=np.zeros((n,3))
 for first in range(0,len(col),4096):
  ids=slice(first,first+4096);np.add.at(out,owners[ids],np.einsum('kab,kb->ka',values[ids],vector[col[ids]]))
 return out
def precondition(vector):
 # 先把尺度右端还原，再按原物理CV/层/分量组织各法向列。
 original=(vector/R)[mapping].reshape(nl,nf,3).transpose(1,0,2).reshape(nf,N);solved=np.einsum('fij,fj->fi',inverse,original)
 physical=solved.reshape(nf,nl,3).transpose(1,0,2).reshape(n,3);out=np.empty_like(physical);out[mapping]=physical;return out
parent=read(s['rhsFile'],['rhs','sourceBinding']);rhs=s['rhs'];np.testing.assert_array_equal(rhs,parent['rhs']);assert parent['sourceBinding']==binding and sha(s['rhsFile'])==s['rhsSHA256']
d=s['detail'];k=d['krylov'];cycle=rows(k['cycles'])[0];Z=cycle['Z'];V=cycle['V'];H=cycle['smallH'];scales=cycle['imageNormalizationScales'];images=np.zeros_like(Z);maximum_map=0.
for j in range(Z.shape[1]):
 images[:,j]=(R*action(Z[:,j].reshape((n,3),order='F'))).ravel(order='F')
 expected=precondition(V[:,j].reshape((n,3),order='F')).ravel(order='F')/scales[j];maximum_map=max(maximum_map,rel(expected,Z[:,j]));assert maximum_map<1e-8
 if (j+1)%8==0:print('INDEPENDENT_ROW_MOMENTUM_COLUMNS '+str(j+1),flush=True)
actual_image=action(s['value']);image_error=rel(actual_image,s['image']);assert image_error<1e-8
remaining=rhs-actual_image;raw_residual=norm(remaining)/norm(rhs);scaled_residual=norm(R*remaining)/norm(R*rhs)
assert abs(raw_residual-d['relativeAResidual'])/max(1,raw_residual)<1e-8 and abs(scaled_residual-d['relativeRowScaledResidual'])/max(1,scaled_residual)<1e-8
tolerance=float(d['options']['relativeTolerance']);assert bool(d['originalUnitTargetMet'])==(d['relativeAResidual']<=tolerance) and bool(d['rowScaledTargetMet'])==(k['flag']==0)
relation=rel(images,V@H);orth=float(np.linalg.norm(V.T@V-np.eye(V.shape[1])));assert relation<1e-8 and orth<1e-8
weights=cycle['weights'];candidate=cycle['xStart']+Z@weights;assert rel(candidate,s['value'].ravel(order='F'))<1e-8
predicted=cycle['rStart']-V@(H@weights);combined=norm((R*remaining).ravel(order='F')-predicted)/max(1,norm(R*rhs));assert combined<1e-8
assert d['AApplications']==d['solverAApplications']+1 and d['originalUnitDiagnosticAApplications']==1 and d['lineInverseApplications']==Z.shape[1]
assert all(sha(ROOT/path)==value for path,value in bindings.items())
out=dict(status='passed-independent-real-row-scale-all-physical-directions-PRinverse-and-dual-residuals',processID=os.getpid(),cellCount=n,scalarRows=3*n,
 frozenFactorByteIdentityAndPriorFullInverseAuditReused=True,newFactorConstructions=0,allPhysicalArnoldiColumns=Z.shape[1],maximumCompletePRinverseMapDifference=maximum_map,
 completeRowSystemArnoldiRelativeDifference=relation,basisOrthogonalityError=orth,rowSystemCombinationRelationRelativeDifference=combined,
 originalAImageRelativeDifference=image_error,originalAResidual=raw_residual,rowScaledResidual=scaled_residual,originalUnitTargetMet=bool(d['originalUnitTargetMet']),rowScaledTargetMet=bool(d['rowScaledTargetMet']),
 originalTargetNotReplacedByScaledFlag=True,originalARowNormRange=[float(row_max.min()),float(row_max.max())],originalARowNormSpan=float(row_max.max()/row_max.min()),
 unscaledReferenceOriginalAResidual=summary['unscaledReferenceOriginalAResidual'],originalAResidualImproved=raw_residual<summary['unscaledReferenceOriginalAResidual'],
 solverAApplications=int(d['solverAApplications']),originalUnitDiagnosticAApplications=1,lineInverseApplications=int(d['lineInverseApplications']),
 newSelectedJApplications=0,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,wholeFiveFieldOrOriginalNSConvergenceProven=False,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
