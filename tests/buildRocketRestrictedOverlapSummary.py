# -*- coding: utf-8 -*-
"""将已测受限重叠研究追加到工作台，保留前十六行和原物理结果身份。"""
from pathlib import Path
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


target = VAL / 'block_linear_research_summary_v11.json'
assert not target.exists()
data = read(VAL / 'block_linear_research_summary_v10.json')
summary_file = VAL / 'overlapFreal01/restricted_overlap_summary.json'
audit_file = VAL / 'overlapFreal01_independent_restricted_overlap.json'
r, audit = read(summary_file), read(audit_file)
assert len(data['studies']) == 16 and audit['status'].startswith('passed-')
assert not r['physicalSolution'] and not r['fullFInternalTargetMet']
assert r['newActualNSEvaluations'] == r['newProductionSolvedAreaM2'] == 0
assert r['sourceBinding']['sourceFileSHA256'] == data['sourceEvidenceSHA256']
assert abs(r['innerFullFReport']['relativeFullFResidual'] - audit['finalFullFResidual']) < 1e-12
k = r['innerFullFReport']['krylov']
row = dict(data['studies'][-1])
row.update(label='完整 A_FF / 受限重叠核心输出', caseID='overlapFreal01',
    arnoldiSteps=k['arnoldiSteps'], CApplications=1, ILUApplications=0,
    innerRelativeCResidual=None, basisOrthogonalityError=k['maximumBasisOrthogonalityError'],
    modelRelationRelativeDifference=k['cycles'][0]['arnoldiRelationRelativeDifference'],
    candidateUsable=False, evidenceFile=str(summary_file), evidenceSHA256=sha(summary_file),
    displayState='受限重叠未改善，F 内部未收敛',
    note='同完整F/右端32步残差0.308766，较原0.278621更差；27461核心组、至多2外邻、1组明确CV对角替代。固定C回代原C为14.5520。全部逆/最终像独立通过，但真实完整基未保存；没有新NS或新增合格面积，旧S未替换。')
data['studies'].append(row)
data['status'] = 'source-bound-restricted-overlap-tested-with-unimproved-complete-F-residual'
data['restrictedOverlapDiagnostics'] = dict(auditFile=str(audit_file), auditSHA256=sha(audit_file),
    fullFResidual=audit['finalFullFResidual'], parentOriginalCResidual=audit['fixedCParentOriginalCResidual'],
    fullFInternalTargetMet=False, originalCAndFResidualDenominatorsKeptDistinct=True,
    numericalOverlapInversesConstructed=True, allCoreUnknownsReturnedExactlyOnce=True,
    allExtendedInverseCoefficientsIndependentlyChecked=True, fullFArnoldiBasisSaved=False,
    independentAllRealFullFArnoldiColumnsClaimed=False, actualFactorArrayBytes=r['actualFactorArrayBytes'],
    actualGraphAndFactorCacheBytes=r['actualCombinedCacheBytes'],
    constructionArrayEstimateBytes=r['factorizationReport']['constructionArrayEstimateBytes'],
    processPeakMemoryClaimed=False, originalNSOrRegionalQualificationProven=False,
    fullFResidualImproved=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    previousPreflightRetainedAsHistorical=True)
data['restrictedOverlapSummaryBuilder'] = dict(processID=os.getpid(), scriptSHA256=sha(__file__))
target.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(status='built-seventeen-source-bound-research-rows', path=str(target), sha256=sha(target))))
