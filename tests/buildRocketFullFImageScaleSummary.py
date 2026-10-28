# -*- coding: utf-8 -*-
"""把完整像列归一的真实研究结果接回工作台，区分F目标与原父矩阵/NS资格。"""
from pathlib import Path
import json
import os
from screenRocketFullFLeastSquares import ROOT, VAL, sha


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


target=VAL/'block_linear_research_summary_v13.json'
assert not target.exists()
d=read(VAL/'block_linear_research_summary_v12.json')
summary_file=VAL/'fullFimagescalereal01/full_F_image_scaled_summary.json'
audit_file=VAL/'fullFimagescalereal01_independent_full_F_image_scale.json'
r,a=read(summary_file),read(audit_file)
assert len(d['studies'])==18 and a['status'].startswith('passed-')
assert r['sourceBinding']['sourceFileSHA256']==d['sourceEvidenceSHA256']
assert not r['physicalSolution'] and r['newActualNSEvaluations']==r['newProductionSolvedAreaM2']==0
f=r['fullFReport'];parent=r['parentDiagnostic']['relativeParentOriginalCResidual']
assert abs(f['relativeFullFResidual']-a['finalFullFResidual'])<1e-12
assert not f['originalCOrNSConvergenceClaimed']
state=('F/父 C 达内部目标，原 NS 待验' if parent<=1e-3 else 'F 内部达标，父 C 仍未达标') if f['fullFInternalTargetMet'] else 'F 显著改善，仍未达目标'
row=dict(d['studies'][-1])
row.update(label='完整 A_FF / 完整像列归一 LSQR',caseID='fullFimagescalereal01',arnoldiSteps=0,CApplications=1,ILUApplications=0,
    innerRelativeCResidual=None,basisOrthogonalityError=None,modelRelationRelativeDifference=None,candidateUsable=False,
    evidenceFile=str(summary_file),evidenceSHA256=sha(summary_file),displayState=state,
    note=f"从原保存F端点以新右侧尺度作{f['nativeIterations']}步LSQR（预算128），原F残差{f['initialRelativeFullFResidual']:.6g}→{f['relativeFullFResidual']:.6g}；原F目标0.001。固定C父矩阵残差{parent:.6g}。全部完整像列范数独立核对，原方程和范数不变；旧S未替换，不代表原NS收敛或新增区域热解。")
d['studies'].append(row)
d['status']='source-bound-complete-image-scaled-F-LSQR-with-separate-parent-and-NS-qualification'
d['fullFImageScaleDiagnostics']=dict(caseID=row['caseID'],algorithm='LSQR',nativeIterations=f['nativeIterations'],
    auditFile=str(audit_file),auditSHA256=sha(audit_file),initialFullFResidual=f['initialRelativeFullFResidual'],
    finalFullFResidual=f['relativeFullFResidual'],fullFInternalTargetMet=f['fullFInternalTargetMet'],
    parentOriginalCResidual=parent,parentOriginalCTargetMet=parent<=1e-3,
    parentCRowsRemainingSquaredNormFraction=(r['parentDiagnostic']['CRowsRelativeToOriginalRHS']/parent)**2,
    columnNormSpan=r['construction']['columnNormSpan'],allRealColumnsIndependentlyCompared=a['allPreconditionedColumns'],
    adjointIdentityRelativeDifference=f['adjointIdentityRelativeDifference'],
    originalFResidualNormAndRightHandSideUnchanged=True,independentAllRealBidiagonalizationColumnsClaimed=False,
    originalNSOrRegionalQualificationProven=False,newActualNSEvaluations=0,newProductionSolvedAreaM2=0)
d['imageScaleSummaryBuilder']=dict(processID=os.getpid(),scriptSHA256=sha(__file__))
target.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status='built-nineteen-source-bound-research-rows',path=str(target),sha256=sha(target))))
