# -*- coding: utf-8 -*-
"""独立核对真实全部法向列逆、完整A作用及保存的每条Arnoldi方向。"""
from pathlib import Path
import os,json,hashlib,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from decimal import Decimal,localcontext
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum();folder=VAL/name;target=VAL/(name+'_independent_momentum_lines.json');assert not target.exists();bindings={}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path,variables=None):
 path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
 return loadmat(path,simplify_cells=True,variable_names=variables)
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def norm(a):
 scale=float(np.max(abs(a))) if np.size(a) else 0.
 return float(np.linalg.norm(a/scale)*scale) if scale else 0.
def rel(a,b):return norm(a-b)/max(1,norm(b))
report_folder=VAL/sys.argv[2] if len(sys.argv)>2 else folder
assert (report_folder/'momentum_lines_real_summary.json').stat().st_size<128*1024,'Use the explicit compact recovery report; never parse oversized basis JSON.'
report=json.loads((report_folder/'momentum_lines_real_summary.json').read_text(encoding='utf-8-sig'));assert report['factorAvailable'] and report['threeActionsComplete']
if report.get('driverReportingInterrupted',False):
 assert report['originalDriverExitCode']==4294967295 and report['summaryRecoveryNewAApplications']==report['summaryRecoveryNewPreconditionerApplications']==0
 stop=json.loads((folder/'agent_reporting_stop_record.json').read_text(encoding='utf-8-sig'))
 assert sha(folder/'momentum_line_krylov_evidence.mat')==stop['savedNumericalEvidenceSHA256'] and not stop['numericalFailureClaimed']
fdata=read(folder/'momentum_line_factor_evidence.mat');f=fdata['factor'];binding=fdata['sourceBinding'];assert binding==report['sourceBinding']
for key,value in fdata['sourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==value
a=read(fdata['assemblyFile'],['momentumMatrix','algebraicOrder','sourceBinding']);m=a['momentumMatrix'];assert a['sourceBinding']==binding
n=int(f['cellCount']);nf=int(f['sourceFaces']);nl=int(f['layers']);N=3*nl;order=a['algebraicOrder'].astype(np.int64)-1;mapping=f['inverseOrder'].astype(np.int64)-1
np.testing.assert_array_equal(order[mapping],np.arange(n));assert nf*nl==n
ptr=m['rowPointer'].astype(np.int64)-1;col=m['column'].astype(np.int64)-1;owners=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr));values=m['values'].transpose(2,0,1)
rowcv=order[owners];columncv=order[col];mask=rowcv%nf==columncv%nf;selected=np.flatnonzero(mask);faces=rowcv[selected]%nf;rl=rowcv[selected]//nf;cl=columncv[selected]//nf
lines=np.zeros((nf,N,N));inverse=f['inverseLine'].transpose(2,0,1)
for i in range(3):
 for j in range(3):lines[faces,3*rl+i,3*cl+j]=values[selected,i,j]
assert len(selected)==f['factorizationReport']['includedNormalLineBlocks'];del rowcv,columncv,faces,rl,cl,mask,selected
fallback_faces=np.asarray(f.get('inverseCertification',{}).get('fallbackFaces',np.zeros(nf))).astype(bool)
for face in np.flatnonzero(fallback_faces):
 original=lines[face].copy();lines[face]=0
 for layer in range(nl):
  ids=slice(3*layer,3*layer+3);lines[face,ids,ids]=original[ids,ids]
assert int(np.count_nonzero(fallback_faces))==f['factorizationReport'].get('fallbackColumnCount',0)
maximum=0.;line_norms=np.zeros(nf);gamma=N*np.finfo(float).eps/(1-N*np.finfo(float).eps);decimal_faces=[]
certificate=f.get('inverseCertification');source_upper=np.asarray(certificate['originalInputResidualUpperBounds']) if certificate else None
if certificate:assert np.all(source_upper<1e-8)
for first in range(0,nf,128):
 ids=slice(first,first+128);product=lines[ids]@inverse[ids];ordinary=np.max(abs(product-np.eye(N)),axis=(1,2));maximum=max(maximum,float(ordinary.max()));line_norms[ids]=np.max(np.sum(abs(inverse[ids]),axis=2),axis=1)
 if certificate:
  error_bound=gamma*np.max(np.sum(abs(lines[ids]),axis=2),axis=1)*np.max(abs(inverse[ids]),axis=(1,2))/(1-gamma)**2
  selected_exact=np.asarray(certificate['accurateCheckedFaces'])[ids].astype(bool)|(ordinary+error_bound>=1e-8)
  decimal_faces.extend((first+np.flatnonzero(selected_exact)).tolist())
 else:assert ordinary.max()<1e-8
np.testing.assert_allclose(line_norms,f['inverseNormInf'],rtol=1e-13,atol=0)
maximum_exact=0.
with localcontext() as context:
 context.prec=80
 for count,face in enumerate(decimal_faces):
  AA=[[Decimal.from_float(float(x)) for x in row] for row in lines[face]];XX=[[Decimal.from_float(float(x)) for x in row] for row in inverse[face]]
  exact=max(abs((Decimal(1) if i==j else Decimal(0))-sum((AA[i][q]*XX[q][j] for q in range(N)),Decimal(0))) for i in range(N) for j in range(N))
  assert exact<Decimal('1e-8') and exact<=Decimal.from_float(float(source_upper[face]));maximum_exact=max(maximum_exact,float(exact))
  if (count+1)%100==0:print('INDEPENDENT_DECIMAL80_NORMAL_LINES '+str(count+1),flush=True)
del lines;gc.collect();print('INDEPENDENT_ALL_NORMAL_LINE_INVERSES_VERIFIED',flush=True)
def precondition(rhs):
 original=rhs[mapping].reshape(nl,nf,3).transpose(1,0,2).reshape(nf,N)
 solved=np.einsum('fij,fj->fi',inverse,original);original=solved.reshape(nf,nl,3).transpose(1,0,2).reshape(n,3);answer=np.empty_like(original);answer[mapping]=original;return answer
def action(vector):
 answer=np.zeros((n,3))
 for first in range(0,len(col),4096):
  ids=slice(first,first+4096);np.add.at(answer,owners[ids],np.einsum('kab,kb->ka',values[ids],vector[col[ids]]))
 return answer
cases=read(folder/'momentum_line_actions.mat',['cases'])['cases'];case_results=[]
for case in rows(cases):
 expected=precondition(case['rhs']);image=action(case['value']);mapping_error=rel(expected,case['value']);action_error=rel(image,case['Avalue']);assert mapping_error<1e-8 and action_error<1e-8
 residual=norm(case['rhs']-image)/norm(case['rhs']);gain=norm(case['value'])/norm(case['rhs']);summary=case['summary']
 assert abs(residual-summary['relativeAResidual'])/max(1,residual)<1e-8 and abs(gain-summary['inverseGainSeconds'])/max(1,gain)<1e-8
 case_results.append(dict(name=summary['name'],completePreconditionerMapRelativeDifference=mapping_error,fullAActionRelativeDifference=action_error,relativeAResidual=residual,inverseGainSeconds=gain))
rhs_data=read(fdata['rhsFile'],['momentumRhs','sourceBinding','algebraicOrder']);assert rhs_data['sourceBinding']==binding;np.testing.assert_array_equal(rhs_data['algebraicOrder'].astype(np.int64)-1,order)
rhs=rhs_data['momentumRhs'];b=rhs.ravel(order='F');inner_result=None
if report['innerSolveComplete']:
 k=read(folder/'momentum_line_krylov_evidence.mat');np.testing.assert_array_equal(rhs,k['rhs']);detail=k['detail'];kr=detail['krylov'];cycle=rows(kr['cycles'])[0]
 candidate=k['value'];remaining=rhs-action(candidate);residual=norm(remaining)/norm(rhs);assert abs(residual-detail['relativeAResidual'])<1e-8
else:
 failure=read(folder/'momentum_arnoldi_failure.mat',['failure'])['failure'];cycle=failure;candidate=failure['candidate'].reshape((n,3),order='F');remaining=rhs-action(candidate);residual=norm(remaining)/norm(rhs);kr=None
Z=cycle['Z'];V=cycle['V'];H=cycle['smallH'];scales=cycle['imageNormalizationScales'];images=np.zeros_like(Z);maximum_map=0.
for j in range(Z.shape[1]):
 x=Z[:,j].reshape((n,3),order='F');images[:,j]=action(x).ravel(order='F');expected=precondition(V[:,j].reshape((n,3),order='F')).ravel(order='F')/scales[j]
 maximum_map=max(maximum_map,rel(expected,Z[:,j]));assert maximum_map<1e-8
 if (j+1)%8==0:print('INDEPENDENT_MOMENTUM_ARNOLDI_COLUMNS '+str(j+1),flush=True)
relation=rel(images,V@H);active=V[:,:min(V.shape[1],Z.shape[1]+1)];orth=float(np.linalg.norm(active.T@active-np.eye(active.shape[1])))
weights=cycle['weights'];combination=cycle['xStart']+Z@weights;assert rel(combination,candidate.ravel(order='F'))<1e-8
predicted=cycle['rStart']-V@(H@weights);actual=remaining.ravel(order='F');combined_error=norm(actual-predicted)/max(1,norm(b))
if report['innerSolveComplete']:assert relation<1e-8 and orth<1e-8 and combined_error<1e-8
else:assert failure['status']=='rejected-by-unchanged-original-Arnoldi-gate' and (failure['basisOrthogonalityError']>=1e-8 or failure['modelRelationRelativeDifference']>=1e-8)
inner_result=dict(completed=bool(report['innerSolveComplete']),originalAResidual=residual,allPhysicalArnoldiColumns=Z.shape[1],allColumnAMapRelativeDifference=relation,maximumCompleteLinePreconditionerMapDifference=maximum_map,
 basisOrthogonalityError=orth,originalACombinationRelationRelativeDifference=combined_error,originalGateThreshold=1e-8,internalTargetMet=bool(report.get('internalMomentumTargetMet',False)),usableAsFixedLinearSchurInverseClaimed=False)
row_max=np.zeros((n,3))
for component in range(3):np.maximum.at(row_max[:,component],owners,np.max(abs(values[:,component,:]),axis=1))
assert np.all(np.isfinite(row_max)&(row_max>0));row_scale=1/row_max
singular=np.linalg.svd(H,compute_uv=False);projection_condition=float(singular[0]/singular[-1]) if singular[-1]>0 else None
inner_result.update(originalARowInfinityNormRange=[float(row_max.min()),float(row_max.max())],originalARowNormSpan=float(row_max.max()/row_max.min()),
 rowScaledResidualDiagnosticOnly=float(norm(row_scale*remaining)/norm(row_scale*rhs)),rowScaledNormUsedAsAcceptance=False,
 savedArnoldiSmallMatrixConditionEstimate=projection_condition,globalConditionOrRankProven=False)
assert all(sha(ROOT/path)==value for path,value in bindings.items())
out=dict(status='passed-independent-all-real-normal-line-inverses-complete-A-and-all-physical-Arnoldi-columns',processID=os.getpid(),cellCount=n,normalLines=nf,unknownsPerLine=N,
 maximumOrdinaryBLASLineProductDefect=maximum,maximumAllOriginalInputResidualUpperBound=float(source_upper.max()) if certificate else maximum,
 decimal80CheckedLines=len(decimal_faces),maximumDecimal80CheckedInverseDefect=maximum_exact,
 originalInverseDefectThreshold=1e-8,ordinaryBLASAndAccurateCertificateKeptDistinct=True,explicitOriginalCVDiagonalFallbackColumns=int(np.count_nonzero(fallback_faces)),
 allOriginalACouplingsKeptInFullAction=True,rejectedFullLineInverseClaimedValid=False,threeFixedActionCases=case_results,innerMomentum=inner_result,newSelectedJApplications=0,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,
 completeFiveFieldOrOriginalNSConvergenceProven=False,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
