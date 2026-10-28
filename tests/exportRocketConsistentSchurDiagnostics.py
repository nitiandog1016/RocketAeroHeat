# -*- coding: utf-8 -*-
"""将已独立核对的动量回代溢出接回工作台，保持原物理场与旧失败位置独立。"""
from pathlib import Path
import os,json,hashlib
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'block_linear_research_summary_v5.json';assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
out=read(VAL/'block_linear_research_summary_v4.json');summary=VAL/'consistentdiagrefine01/consistent_schur_real_summary.json';r=read(summary)
audit_file=VAL/'consistentdiagrefine01_independent_consistent_schur.json';audit=read(audit_file);failure=audit['fixedMomentumFailure'];f=r['failedMomentumObservation']['firstFailure']
assert audit['status'].startswith('passed-') and audit['diagnosticOnlyRefinementKeepsAllOtherPayloadValues']
for case in ('consistentschurreal02','consistentdiagrefine01'):
 p=read(VAL/(case+'_process.json'));assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
assert r['sourceBinding']['sourceFileSHA256']==out['sourceEvidenceSHA256'] and r['momentumMapRejected'] and r['consistentSchurActions']==0
assert failure['firstAlgebraicRow']==f['algebraicRow']==167894 and failure['originalCV']==r['failedOriginalCV']==4570
assert failure['independentFirstNonfiniteStageAndRowExactlyMatch'] and not audit['completeGlobalNumericalRankOrUsableInverseProven']
note=('原CV4570 / 原面4570 / 第1层：前代第167894行首次非有限，独立重放一致。首次凝聚动量逆增益约'+
      format(audit['sameQCondensedMomentumInverseGainSeconds'],'.3g')+'秒；未完成任何一致Schur作用。该位置属于代数回代诊断，不表示物理面元网格故障；不能沿用为旧标量零主元位置。')
out['studies'].append(dict(label='三动量因子 / 一致Schur前置检查',caseID='consistentdiagrefine01',linearSystem='C=-Jgraph',
 arnoldiSteps=0,CApplications=0,ILUApplications=0,innerRelativeCResidual=None,basisOrthogonalityError=None,modelRelationRelativeDifference=None,
 orthogonalityThreshold=1e-8,modelRelationThreshold=1e-8,candidateUsable=False,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,
 evidenceFile=str(summary),evidenceSHA256=sha(summary),displayState='动量回代溢出（已定位）',note=note))
out['fixedMomentumInverseDiagnostic']=dict(auditFile=str(audit_file),auditSHA256=sha(audit_file),
 firstAlgebraicRow=failure['firstAlgebraicRow'],originalCV=failure['originalCV'],originalFace=failure['originalFace'],layer=failure['layer'],
 stage=failure['stage'],sameQMomentumInverseGainSeconds=audit['sameQCondensedMomentumInverseGainSeconds'],
 rawMomentumRelativeResidual=audit['sameQRawMomentumResidualOfFactorSolve'],localScaledDiagonalSingularRatioRange=audit['localScaledDiagonalSingularRatioRange'],
 completeGlobalRankProven=False,originalPhysicalGeometryFaultProven=False,oldNativeScalarZeroPivotLocationStillUnknown=True,
 consistentSchurActionsCompleted=0,fixedMomentumSolveAttempts=2,physicalStateChanged=False)
out['status']='source-bound-fixed-momentum-overflow-and-consistent-Schur-prerequisite-rejection';out['diagnosticBuilderProcessID']=os.getpid();out['builderSHA256']=sha(__file__)
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],researchRows=len(out['studies']),newActualNSEvaluations=0)))
