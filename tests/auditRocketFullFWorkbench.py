# -*- coding: utf-8 -*-
"""核对十六条来源、F/C范围区分、失败身份保留及原面全部场量字节不变。"""
from pathlib import Path
import os, json, hashlib
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
folder = VAL/'workbench26'
target = folder/'independent_full_F_study_view.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


data = read(folder/'wall-velocity-study-data.json')
old = read(VAL/'workbench25/wall-velocity-study-data.json')
research, parent = data.pop('linearResearch'), old.pop('linearResearch')
assert data == old and len(research['studies']) == 16
assert research['studies'][:14] == parent['studies'] and research == read(VAL/'block_linear_research_summary_v10.json')
assert research['sourceEvidenceSHA256'] == sha(VAL/'krylov17/wall_velocity_krylov_evidence.mat')
assert not research['physicalStateChanged'] and not research['productionQualified'] and research['newProductionSolvedAreaM2'] == 0
bindings = {}
for row in research['studies']:
    file = Path(row['evidenceFile'])
    file = file if file.is_absolute() else ROOT/file
    assert sha(file) == row['evidenceSHA256'] and not row['candidateUsable']
    assert row['newActualNSEvaluations'] == row['newProductionSolvedAreaM2'] == 0
    bindings[str(file.relative_to(ROOT))] = sha(file)
for row in research['studies'][-2:]:
    assert row['innerRelativeCResidual'] is None and row['CApplications'] == 1
    assert row['linearSystem'] == 'exact-complete-spatial-A_FF-with-separate-fixed-C-parent-diagnostic'
    assert row['basisOrthogonalityError'] < 1e-8 and row['modelRelationRelativeDifference'] < 1e-8
for key in ('fullFDiagnostics', 'fullFStrengthDiagnostics'):
    d = research[key]
    file = Path(d['auditFile'])
    assert sha(file) == d['auditSHA256'] and read(file)['status'].startswith('passed-')
    assert d['originalCAndFResidualDenominatorsKeptDistinct'] and not d['fullFInternalTargetMet']
    assert not d['fullFArnoldiBasisSaved'] and not d['independentAllRealFullFArnoldiColumnsClaimed']
failure = research['fullFStrengthDiagnostics']['retainedShapeFailure']
assert sha(failure['evidenceFile']) == failure['evidenceSHA256'] and failure['driverExitCode'] == 1
assert not failure['iterationHadStarted'] and failure['repairedCoverageRowAndColumnTestPassed']
overlap = research['restrictedOverlapPreflight']
assert sha(overlap['evidenceFile']) == overlap['evidenceSHA256']
assert not overlap['completeConstructionOrProcessBudgetProven'] and not overlap['numericalOverlapInversesConstructed']
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent/name
        before, after = VAL/'workbench25'/relative, folder/relative
        assert sha(before) == sha(after)
        payload.append(dict(path=str(after.relative_to(ROOT)), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding_test = loadmat(VAL/'fullFui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert binding_test['status'].startswith('passed-') and all(binding_test['rejected'])
p = read(VAL/'fullFui01_process.json')
assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
out = dict(status='passed-independent-sixteen-source-bound-rows-separated-full-F-and-unchanged-facet-data', processID=os.getpid(),
    workbench='outputs/v83_validation/workbench26/index.html', researchRows=16,
    originalNSAndAcceptanceHistoryExactlyUnchanged=True, originalCAndFResidualDenominatorsAndSystemScopesKeptDistinct=True,
    diagnosticShapeFailureAndRepairKeptExplicit=True, structuralOverlapNotClaimedAsNumericalSolver=True,
    allThreeFalseSourceOrQualificationExportsRejected=True, productionQualificationNotClaimed=True,
    sourceFacetPayload=payload, previousFullMATFieldAuditReusedByByteIdentity=True, browserVisualOrWebGLRuntimeVerified=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
