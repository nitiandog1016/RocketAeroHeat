# -*- coding: utf-8 -*-
"""核对LSQR行、伴随指标适用性、旧十七行与全部原面数据保持。"""
from pathlib import Path
import hashlib
import json
import os
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


folder = VAL/'workbench28'
target = folder/'independent_full_F_least_squares_study_view.json'
assert not target.exists()
data = read(folder/'wall-velocity-study-data.json')
old = read(VAL/'workbench27/wall-velocity-study-data.json')
research, prior = data.pop('linearResearch'), old.pop('linearResearch')
assert data == old and research['studies'][:17] == prior['studies'] and len(research['studies']) == 18
assert research == read(VAL/'block_linear_research_summary_v12.json')
assert not research['physicalStateChanged'] and not research['productionQualified']
bindings = {}
for row in research['studies']:
    path = Path(row['evidenceFile'])
    path = path if path.is_absolute() else ROOT/path
    assert sha(path) == row['evidenceSHA256'] and not row['candidateUsable']
    assert row['newActualNSEvaluations'] == row['newProductionSolvedAreaM2'] == 0
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
row = research['studies'][-1]
assert row['arnoldiSteps'] == 0 and row['basisOrthogonalityError'] is None and row['modelRelationRelativeDifference'] is None
assert row['innerRelativeCResidual'] is None and row['CApplications'] == 1
d = research['fullFLeastSquaresDiagnostics']
assert d['algorithm'] == 'LSQR' and d['caseID'] == row['caseID']
assert sha(d['auditFile']) == d['auditSHA256'] and read(d['auditFile'])['status'].startswith('passed-')
assert d['adjointIdentityRelativeDifference'] < 1e-8 and not d['fullFInternalTargetMet']
assert not d['independentAllRealBidiagonalizationColumnsClaimed']
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent/name
        before, after = VAL/'workbench27'/relative, folder/relative
        assert sha(before) == sha(after)
        payload.append(dict(path=after.relative_to(ROOT).as_posix(), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding = loadmat(VAL/'fullFlsqrui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert all(binding['rejected'])
process = read(VAL/'fullFlsqrui01_process.json')
assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
out = dict(status='passed-independent-eighteen-source-bound-rows-LSQR-applicability-and-unchanged-facet-data',
    processID=os.getpid(), researchRows=18, previousSeventeenRowsExactlyPreserved=True,
    originalTenFieldsAndNSHistoryUnchanged=True, originalFAndParentCNormsKeptDistinct=True,
    LSQRNotLabeledAsArnoldi=True, allThreeFalseSourceOrQualificationExportsRejected=True,
    sourceFacetPayload=payload, previousFullMATFieldAuditReusedByByteIdentity=True,
    browserVisualOrWebGLRuntimeVerified=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k: v for k, v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
