# -*- coding: utf-8 -*-
"""独立核对十四行研究来源、原C/S指标分开和所有原面场量/历史字节不变。"""
from pathlib import Path
import os, json, hashlib
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
folder = VAL/'workbench25'
target = folder/'independent_spatial_reduction_study_view.json'
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
old = read(VAL/'workbench24/wall-velocity-study-data.json')
research, parent = data.pop('linearResearch'), old.pop('linearResearch')
assert data == old and len(research['studies']) == 14
assert research['studies'][:12] == parent['studies'] and research == read(VAL/'block_linear_research_summary_v9.json')
assert research['sourceEvidenceSHA256'] == sha(VAL/'krylov17/wall_velocity_krylov_evidence.mat')
assert not research['physicalStateChanged'] and not research['productionQualified'] and research['newProductionSolvedAreaM2'] == 0
bindings = {}
for entry in research['studies']:
    evidence = Path(entry['evidenceFile'])
    evidence = evidence if evidence.is_absolute() else ROOT/evidence
    assert sha(evidence) == entry['evidenceSHA256']
    assert not entry['candidateUsable'] and entry['newActualNSEvaluations'] == entry['newProductionSolvedAreaM2'] == 0
    bindings[str(evidence.relative_to(ROOT))] = sha(evidence)
main, schur = research['studies'][-2:]
assert main['caseID'] == 'spatialreductionreal01' and main['innerRelativeCResidual'] > 1e-3 and main['CApplications'] == 33
assert schur['caseID'] == 'spatialschurreal01' and schur['innerRelativeCResidual'] is None and schur['CApplications'] == 1
assert main['basisOrthogonalityError'] < 1e-8 and schur['basisOrthogonalityError'] < 1e-8
assert 'S残差' in schur['note'] and '单次代回原C残差' in schur['note']
for key in ('spatialReductionDiagnostics', 'spatialSchurDiagnostic'):
    diagnostic = research[key]
    audit = Path(diagnostic['auditFile'])
    assert sha(audit) == diagnostic['auditSHA256'] and read(audit)['status'].startswith('passed-')
    bindings[str(audit.relative_to(ROOT))] = sha(audit)
d = research['spatialSchurDiagnostic']
assert d['originalSAndCResidualDenominatorsKeptDistinct'] and d['all32ApproximateSPhysicalDirectionsIndependentlyVerified']
assert d['FRowsRemainingSquaredNormFraction'] > .999 and not d['approximateSTargetMet'] and not d['originalCInternalTargetMet']
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent/name
        before, after = VAL/'workbench24'/relative, folder/relative
        assert sha(before) == sha(after)
        payload.append(dict(path=str(after.relative_to(ROOT)), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding_test = loadmat(VAL/'spatialreductionui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert binding_test['status'].startswith('passed-') and all(binding_test['rejected'])
process = read(VAL/'spatialreductionui01_process.json')
assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
out = dict(status='passed-independent-fourteen-source-bound-rows-separated-C-S-and-unchanged-facet-data', processID=os.getpid(),
    workbench='outputs/v83_validation/workbench25/index.html', researchRows=14,
    originalNSAndAcceptanceHistoryExactlyUnchanged=True, originalCAndApproximateSResidualDenominatorsKeptDistinct=True,
    approximateSNativeAndIndependentArnoldiGatesNotClaimedForOriginalC=True,
    allThreeFalseSourceOrQualificationExportsRejected=True, productionQualificationNotClaimed=True,
    sourceFacetPayload=payload, previousFullMATFieldAuditReusedByByteIdentity=True,
    browserVisualOrWebGLRuntimeVerified=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
