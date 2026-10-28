# -*- coding: utf-8 -*-
"""独立核对新增重叠研究身份、前十六行与六份原面数据字节不变。"""
from pathlib import Path
import json
import os
import hashlib
from scipy.io import loadmat

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

folder = VAL / 'workbench27'
target = folder / 'independent_restricted_overlap_study_view.json'
assert not target.exists()
data = read(folder / 'wall-velocity-study-data.json')
old = read(VAL / 'workbench26/wall-velocity-study-data.json')
research, previous = data.pop('linearResearch'), old.pop('linearResearch')
assert data == old and len(research['studies']) == 17
assert research['studies'][:16] == previous['studies']
assert research == read(VAL / 'block_linear_research_summary_v11.json')
assert not research['physicalStateChanged'] and not research['productionQualified']
assert research['newProductionSolvedAreaM2'] == research['newActualNSEvaluations'] == 0
bindings = {}
for row in research['studies']:
    file = Path(row['evidenceFile'])
    file = file if file.is_absolute() else ROOT / file
    assert sha(file) == row['evidenceSHA256'] and not row['candidateUsable']
    assert row['newActualNSEvaluations'] == row['newProductionSolvedAreaM2'] == 0
    bindings[str(file.relative_to(ROOT))] = sha(file)
row = research['studies'][-1]
assert row['innerRelativeCResidual'] is None and row['CApplications'] == 1
assert row['linearSystem'] == 'exact-complete-spatial-A_FF-with-separate-fixed-C-parent-diagnostic'
assert row['basisOrthogonalityError'] < 1e-8 and row['modelRelationRelativeDifference'] < 1e-8
d = research['restrictedOverlapDiagnostics']
assert sha(d['auditFile']) == d['auditSHA256'] and read(d['auditFile'])['status'].startswith('passed-')
assert d['numericalOverlapInversesConstructed'] and not d['fullFInternalTargetMet'] and not d['fullFResidualImproved']
assert d['originalCAndFResidualDenominatorsKeptDistinct'] and not d['independentAllRealFullFArnoldiColumnsClaimed']
assert d['previousPreflightRetainedAsHistorical'] and research['restrictedOverlapPreflight'] == previous['restrictedOverlapPreflight']
assert research['sourceEvidenceSHA256'] == sha(VAL / 'krylov17/wall_velocity_krylov_evidence.mat')
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent / name
        before, after = VAL / 'workbench26' / relative, folder / relative
        assert sha(before) == sha(after)
        payload.append(dict(path=str(after.relative_to(ROOT)), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder / 'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560
assert node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding = loadmat(VAL / 'overlapFui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert binding['status'].startswith('passed-') and all(binding['rejected'])
process = read(VAL / 'overlapFui01_process.json')
assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
out = dict(status='passed-independent-seventeen-source-bound-rows-with-restricted-overlap-and-unchanged-facet-data',
    processID=os.getpid(), researchRows=17, originalNSAndAcceptanceHistoryExactlyUnchanged=True,
    originalCAndFResidualDenominatorsAndSystemScopesKeptDistinct=True,
    allThreeFalseSourceOrQualificationExportsRejected=True,
    allExistingSixteenRowsAndHistoricalPreflightExactlyPreserved=True,
    newOverlapExplicitlyUnconverged=True, productionQualificationNotClaimed=True,
    sourceFacetPayload=payload, previousFullMATFieldAuditReusedByByteIdentity=True,
    browserVisualOrWebGLRuntimeVerified=False, newActualNSEvaluations=0,
    newProductionSolvedAreaM2=0, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps({k: v for k, v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
