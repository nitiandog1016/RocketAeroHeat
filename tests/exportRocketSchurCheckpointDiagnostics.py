# -*- coding: utf-8 -*-
"""把已保存的真实主元、全S复算及新零主元拒绝接回同源科研工作台。"""
from pathlib import Path
import os,json,hashlib
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'block_linear_research_summary_v3.json';assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
out=read(VAL/'block_linear_research_summary_v2.json');receipt=VAL/'schurcheckpointreal01/schur_failure_checkpoint_summary.json';r=read(receipt);f=r['failureSummary'];auditFile=VAL/'schurcheckpointreal01_independent_real_checkpoint.json';audit=read(auditFile)
assert f['sourceBinding']['sourceFileSHA256']==out['sourceEvidenceSHA256'] and f['algebraicRow']==188887 and f['originalControlVolume']==61925
assert f['maximumInverseDefect']>=f['originalInverseDefectThreshold']==1e-8 and f['ledgerExactlyReproducesPartialRow'] and not f['preconditionerUsable']
assert audit['status'].startswith('passed-') and audit['allSchurCoefficientRelativeDifference']<1e-10
details=out['failedCellDetails'];details.update(localInverseDefectValue=f['maximumInverseDefect'],pivotRcond=f['pivotRcond'],pivotGrowthFromOriginalDiagonal=audit['pivotNormGrowthFromOriginalSchurDiagonal'],
 warning='This location belongs to the original fixed 2-field block rejection. The new scalar ILUTP zero-pivot index was not reported by the library.')
out['pivotFailureCheckpoint']=dict(receiptFile=str(receipt),receiptSHA256=sha(receipt),auditFile=str(auditFile),auditSHA256=sha(auditFile),
 allRealSchurInputCoefficientsIndependentlyChecked=True,allSchurCoefficientRelativeDifference=audit['allSchurCoefficientRelativeDifference'],oldBlockStillRejected=True)
old=out['studies'][4];old.update(caseID='schurcheckpointreal01',evidenceFile=str(receipt),evidenceSHA256=sha(receipt),displayState='原分解拒绝（已定位）',
 note='原CV61925 / 原面2669 / 第4层：逆检查误差4.76837e−7，原门槛1e−8；rcond约6.67e−11。主元范数较原对角增长约4.04e9倍；真实主元与消元账已保存，全部S输入块已独立复算。')
nativeFile=VAL/'pivotedschurreal02/pivoted_schur_real_summary.json';native=read(nativeFile);process=read(VAL/'pivotedschurreal02_process.json');assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
assert native['sourceBinding']==f['sourceBinding'] and native['candidateRejected'] and not native['factorReady'] and native['errorIdentifier']=='MATLAB:ilu:ZeroPivot'
out['studies'].append(dict(label='CV相邻排列 / 行主元标量Schur',caseID='pivotedschurreal02',linearSystem='C=-Jgraph',arnoldiSteps=0,CApplications=0,ILUApplications=0,
 innerRelativeCResidual=None,basisOrthogonalityError=None,modelRelationRelativeDifference=None,orthogonalityThreshold=1e-8,modelRelationThreshold=1e-8,candidateUsable=False,
 newActualNSEvaluations=0,newProductionSolvedAreaM2=0,evidenceFile=str(nativeFile),evidenceSHA256=sha(nativeFile),displayState='零主元拒绝',
 note='ILUTP阈值0.01；禁止补对角/替换零主元。未生成因子，未开始右端或GMRES；库未返回零主元位置，不沿用原块的面元位置。'))
scaleFile=VAL/'schur_scalar_scales01.json';scale=read(scaleFile);assert scale['zeroScalarRows']==scale['zeroScalarColumns']==0
out['rawSchurScalarDiagnostics']=dict(sourceFile=str(scaleFile),sourceSHA256=sha(scaleFile),zeroRows=0,zeroColumns=0,
 structuralRank=scale['provenStructuralRank'],numericalFullRankClaimed=False,rowNormSpan=scale['rowInfinityNormSpan']['maximumToMinimum'],columnNormSpan=scale['columnInfinityNormSpan']['maximumToMinimum'],
 initialColumnThresholdPotentialSmallFraction=scale['potentialFraction'],actualZeroPivotCauseProved=False)
out['status']='source-bound-actual-pivot-checkpoint-and-zero-pivot-candidate-research-display';out['diagnosticBuilderProcessID']=os.getpid();out['builderSHA256']=sha(__file__)
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],studies=len(out['studies']),knownOriginalBlockInverseDefect=details['localInverseDefectValue']),ensure_ascii=False))
