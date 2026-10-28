# -*- coding: utf-8 -*-
"""核对十一条研究身份、A/C分列、原面数据不变及实际导出拒绝。"""
from pathlib import Path
import os, json, hashlib
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
folder = VAL / 'workbench23'
target = folder / 'independent_coupled_patch_study_view.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


data = read(folder / 'wall-velocity-study-data.json')
old = read(VAL / 'workbench22/wall-velocity-study-data.json')
research = data.pop('linearResearch')
parent = old.pop('linearResearch')
assert data == old
assert research == read(VAL / 'block_linear_research_summary_v7.json')
assert research['studies'][:9] == parent['studies'] and len(research['studies']) == 11
assert research['sourceEvidenceSHA256'] == sha(VAL / 'krylov17/wall_velocity_krylov_evidence.mat')
assert not research['physicalStateChanged'] and not research['productionQualified']
assert research['newProductionSolvedAreaM2'] == 0
bindings = {}
for entry in research['studies']:
    source = Path(entry['evidenceFile'])
    source = source if source.is_absolute() else ROOT / source
    assert sha(source) == entry['evidenceSHA256']
    assert not entry['candidateUsable'] and entry['newActualNSEvaluations'] == entry['newProductionSolvedAreaM2'] == 0
    bindings[str(source.relative_to(ROOT))] = sha(source)
row, coupled = research['studies'][-2:]
assert row['caseID'] == 'rowmomentumreal01' and row['innerRelativeCResidual'] is None and row['CApplications'] == 0
assert row['displayState'] == '原单位残差未改善'
row_d = research['rowMomentumScaleDiagnostics']
assert row_d['originalAResidual'] > 1 and not row_d['originalUnitTargetMet'] and not row_d['rowScaledTargetMet']
assert row_d['originalTargetNotReplacedByScaledFlag']
assert coupled['caseID'] == 'coupledpatchreal01' and coupled['CApplications'] == 33
assert coupled['innerRelativeCResidual'] > 1e-3 and coupled['displayState'] == '完整C内部未收敛'
assert coupled['basisOrthogonalityError'] < 1e-8 and coupled['modelRelationRelativeDifference'] < 1e-8
coupled_d = research['coupledPatchDiagnostics']
assert not coupled_d['realFullArnoldiBasisSaved'] and not coupled_d['independentAllRealArnoldiColumnAuditClaimed']
assert not coupled_d['internalCTargetMet'] and not coupled_d['completeFiveFieldOrOriginalNSConvergenceProven']
for diagnostic in (row_d, coupled_d):
    audit = Path(diagnostic['auditFile'])
    assert sha(audit) == diagnostic['auditSHA256'] and read(audit)['status'].startswith('passed-')
    bindings[str(audit.relative_to(ROOT))] = sha(audit)
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent / name
        before = VAL / 'workbench22' / relative
        after = folder / relative
        assert sha(before) == sha(after)
        payload.append(dict(path=str(after.relative_to(ROOT)), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder / 'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding_test = loadmat(VAL / 'coupledpatchui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert binding_test['status'].startswith('passed-') and all(binding_test['rejected'])
process = read(VAL / 'coupledpatchui01_process.json')
assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
out = dict(status='passed-independent-eleven-source-bound-rows-original-A-C-scope-and-unchanged-facet-data',
    processID=os.getpid(), workbench='outputs/v83_validation/workbench23/index.html', researchRows=11,
    originalAAndScaledResidualsNotMislabelledAsC=True, internalTargetFailuresKeptExplicit=True,
    realFullArnoldiColumnAuditNotClaimed=True, allThreeFalseSourceOrQualificationExportsRejected=True,
    productionQualificationNotClaimed=True, sourceFacetPayload=payload,
    previousFullMATFieldAuditReusedByByteIdentity=True, originalNSAndAcceptanceHistoryExactlyUnchanged=True,
    browserVisualOrWebGLRuntimeVerified=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
