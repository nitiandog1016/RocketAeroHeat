# -*- coding: utf-8 -*-
"""把真实MATLAB LSQR与独立核对接回科研工作台，保留原十七行与场量身份。"""
from pathlib import Path
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


target = VAL/'block_linear_research_summary_v12.json'
assert not target.exists()
data = read(VAL/'block_linear_research_summary_v11.json')
summary_file = VAL/'fullFlsqrreal02/full_F_least_squares_summary.json'
audit_file = VAL/'fullFlsqrreal02_independent_full_F_least_squares.json'
r, audit = read(summary_file), read(audit_file)
assert len(data['studies']) == 17 and audit['status'].startswith('passed-')
assert not r['physicalSolution'] and r['newActualNSEvaluations'] == r['newProductionSolvedAreaM2'] == 0
assert r['sourceBinding']['sourceFileSHA256'] == data['sourceEvidenceSHA256']
f = r['fullFReport']
assert not f['fullFInternalTargetMet'] and abs(f['relativeFullFResidual']-audit['finalRelativeFullFResidual']) < 1e-12
row = dict(data['studies'][-1])
row.update(label='完整 A_FF / LSQR 全局校正', caseID='fullFlsqrreal02',
    arnoldiSteps=0, CApplications=1, ILUApplications=0, innerRelativeCResidual=None,
    basisOrthogonalityError=None, modelRelationRelativeDifference=None, candidateUsable=False,
    evidenceFile=str(summary_file), evidenceSHA256=sha(summary_file),
    displayState='F 有所改善，父 C 未改善',
    note=f"从同一最佳F端点作{f['nativeIterations']}步LSQR，原F残差{f['initialRelativeFullFResidual']:.6g}→{f['relativeFullFResidual']:.6g}，仍未达0.001；固定C父残差9.65848→{r['parentDiagnostic']['relativeParentOriginalCResidual']:.6g}，反而变差。原F/父C分母分开。伴随和最终像已核对；LSQR不适用Arnoldi指标，未逐列独立核对真实双对角化。没有新NS或新增合格面积。")
data['studies'].append(row)
data['status'] = 'source-bound-real-complete-F-LSQR-improved-but-unresolved'
data['fullFLeastSquaresDiagnostics'] = dict(caseID=row['caseID'], algorithm='LSQR',
    nativeIterations=f['nativeIterations'], auditFile=str(audit_file), auditSHA256=sha(audit_file),
    initialFullFResidual=f['initialRelativeFullFResidual'], finalFullFResidual=f['relativeFullFResidual'],
    parentOriginalCResidual=r['parentDiagnostic']['relativeParentOriginalCResidual'],
    parentOriginalCResidualImproved=False,
    adjointIdentityRelativeDifference=f['adjointIdentityRelativeDifference'],
    fullFInternalTargetMet=False, independentAllRealBidiagonalizationColumnsClaimed=False,
    originalCAndFResidualDenominatorsKeptDistinct=True, originalNSOrRegionalQualificationProven=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
data['leastSquaresSummaryBuilder'] = dict(processID=os.getpid(), scriptSHA256=sha(__file__))
target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status='built-eighteen-source-bound-research-rows', path=str(target), sha256=sha(target))))
