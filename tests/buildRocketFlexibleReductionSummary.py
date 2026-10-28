# -*- coding: utf-8 -*-
"""将完整父C可变内层研究接回工作台，原C内部资格与原NS/区域资格分开。"""
from pathlib import Path
import json,os
from screenRocketFullFLeastSquares import ROOT,VAL,sha


def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))


target=VAL/'block_linear_research_summary_v14.json';assert not target.exists()
d=read(VAL/'block_linear_research_summary_v13.json')
file=VAL/'flexiblereductionreal01/flexible_reduction_summary.json';audit_file=VAL/'flexiblereductionreal01_independent_flexible_reduction.json'
r,a=read(file),read(audit_file)
assert r['innerRunComplete'] and not r['candidateRejected'] and a['status'].startswith('passed-')
assert len(d['studies'])==19 and r['sourceBinding']['sourceFileSHA256']==d['sourceEvidenceSHA256']
assert not r['physicalSolution'] and r['newActualNSEvaluations']==r['newProductionSolvedAreaM2']==0
f=r['innerReport'];k=f['krylov'];passed=bool(f['internalCTargetMet'])
assert abs(f['relativeCResidual']-a['result']['finalOriginalCRelativeResidual'])<1e-8
row=dict(d['studies'][-1])
row.update(label='完整父 C / 可变 F 内层耦合',caseID='flexiblereductionreal01',linearSystem='complete-original-C-fixed-operator-flexible-right-preconditioning',
    arnoldiSteps=k['arnoldiSteps'],CApplications=f['CApplications']+r['finalCDiagnosticApplications'],ILUApplications=0,
    innerRelativeCResidual=f['relativeCResidual'],basisOrthogonalityError=k['maximumBasisOrthogonalityError'],
    modelRelationRelativeDifference=k['cycles'][0]['arnoldiRelationRelativeDifference'],candidateUsable=passed,
    evidenceFile=str(file),evidenceSHA256=sha(file),displayState='原 C 内部达标，原 NS 待验' if passed else '完整父 C 仍未达标',
    note=f"全部原C未知量重新参与外层，{k['arnoldiSteps']}个实际Z方向、每次两条最多32步归一F内解；原C残差{f['relativeCResidual']:.6g}，目标0.001。{f['FInnerTargetsNotMet']}次F内解未达其目标，仅作近似预条件；旧S不声称恒等。全部实际外层列/块回代独立核对，未核对真实全部内层LSQR列。无新NS/物理解面积。")
d['studies'].append(row);d['status']='source-bound-complete-parent-C-flexible-right-preconditioned-study'
d['flexibleReductionDiagnostics']=dict(caseID=row['caseID'],algorithm='FGMRES',relativeOriginalCResidual=f['relativeCResidual'],
    internalCTargetMet=passed,allRealOuterArnoldiColumnsIndependentlyVerified=True,actualOuterArnoldiSteps=k['arnoldiSteps'],
    FInnerTargetsNotMet=f['FInnerTargetsNotMet'],fixedLinearPreconditionerClaimed=False,completeParentOperatorFixed=True,
    oldSchurIdentityClaimed=False,independentAllRealInnerLSQRColumnsClaimed=False,auditFile=str(audit_file),auditSHA256=sha(audit_file),
    originalNSOrRegionalQualificationProven=False,newActualNSEvaluations=0,newProductionSolvedAreaM2=0)
d['flexibleSummaryBuilder']=dict(processID=os.getpid(),scriptSHA256=sha(__file__))
target.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status='built-twenty-source-bound-research-rows',internalCTargetMet=passed,path=str(target),sha256=sha(target))))
