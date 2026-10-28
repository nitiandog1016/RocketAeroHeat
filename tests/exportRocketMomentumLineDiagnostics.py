# -*- coding: utf-8 -*-
"""把原A的有限但未收敛研究接回工作台，原C指标保持未知，避免错列变量。"""
from pathlib import Path
import os,json,hashlib
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'block_linear_research_summary_v6.json';assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
out=read(VAL/'block_linear_research_summary_v5.json');summary_file=VAL/'momentumlinereport02/momentum_lines_real_summary.json';assert summary_file.stat().st_size<128*1024;r=read(summary_file)
audit_file=VAL/'momentumlinereal03_independent_momentum_lines.json';audit=read(audit_file);inner=audit['innerMomentum'];factor=r['factorizationReport']
assert audit['status'].startswith('passed-') and inner['completed'] and not inner['internalTargetMet']
assert not r['internalMomentumTargetMet'] and not r['wholeFiveFieldPreconditionerQualified'] and r['sourceBinding']['sourceFileSHA256']==out['sourceEvidenceSHA256']
for case,code in [('momentumlinereal03',4294967295),('momentumlinereport02',0)]:
 p=read(VAL/(case+'_process.json'));assert p['status']=='terminal' and p['exitCode']==code and all(p['sourceStableAfter'].values())
assert r['summaryRecoveryNewAApplications']==r['summaryRecoveryNewPreconditionerApplications']==0
note=('19752列、每列33量；2列显式采用原CV对角预条件，原A全部耦合仍参与。32方向后原A相对剩余'+
 format(inner['originalAResidual'],'.9g')+'，未达内部1e-3；不是完整C或原NS收敛。数值MAT已保存，报告写入中止后仅恢复精简摘要，没有重算。')
out['studies'].append(dict(label='原A三动量 / 法向列与CV对角预条件',caseID='momentumlinereal03',linearSystem='original-three-momentum-A-subblock-of-C',
 arnoldiSteps=32,CApplications=0,ILUApplications=0,innerRelativeCResidual=None,basisOrthogonalityError=r['innerReport']['krylov']['maximumBasisOrthogonalityError'],
 modelRelationRelativeDifference=inner['originalACombinationRelationRelativeDifference'],orthogonalityThreshold=1e-8,modelRelationThreshold=1e-8,candidateUsable=False,
 newActualNSEvaluations=0,newProductionSolvedAreaM2=0,evidenceFile=str(summary_file),evidenceSHA256=sha(summary_file),displayState='有限但内部未收敛',note=note))
out['normalMomentumLineDiagnostics']=dict(auditFile=str(audit_file),auditSHA256=sha(audit_file),rawAResidual=inner['originalAResidual'],rawAInnerTarget=1e-3,
 internalTargetMet=False,normalLineCount=19752,unknownsPerLine=33,originalCVDiagonalFallbackColumns=2,allOriginalACouplingsKept=True,
 maximumCertifiedOriginalInputInverseDefect=factor['originalUnitMaximumInverseDefect'],maximumOrdinaryBLASProductDefect=factor['maximumOrdinaryBLASDefect'],
 decimal80CheckedLines=audit['decimal80CheckedLines'],rejectedFullLineInversesClaimedValid=False,
 rawARowNormSpan=inner['originalARowNormSpan'],rowScaledResidualDiagnosticOnly=inner['rowScaledResidualDiagnosticOnly'],rowScaledNormUsedForAcceptance=False,
 numericalRunAApplications=37,numericalRunLineInverseApplications=35,summaryRecoveryNewAApplications=0,newActualNSEvaluations=0,
 fullFiveFieldPreconditionerQualified=False,physicalStateChanged=False,browserVisualOrWebGLRuntimeVerified=False)
out['status']='source-bound-finite-normal-momentum-preconditioner-with-unresolved-original-A-residual';out['diagnosticBuilderProcessID']=os.getpid();out['builderSHA256']=sha(__file__)
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],researchRows=len(out['studies']))))
