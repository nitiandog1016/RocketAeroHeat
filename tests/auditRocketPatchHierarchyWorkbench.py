# -*- coding: utf-8 -*-
"""独立核对十二行来源、七层未收敛身份及原面全部场量/历史字节不变。"""
from pathlib import Path
import os, json, hashlib
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
folder = VAL/'workbench24'
target = folder/'independent_patch_hierarchy_study_view.json'
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
old = read(VAL/'workbench23/wall-velocity-study-data.json')
research, parent = data.pop('linearResearch'), old.pop('linearResearch')
assert data == old and len(research['studies']) == 12
assert research['studies'][:11] == parent['studies']
assert research == read(VAL/'block_linear_research_summary_v8.json')
assert research['sourceEvidenceSHA256'] == sha(VAL/'krylov17/wall_velocity_krylov_evidence.mat')
assert not research['physicalStateChanged'] and not research['productionQualified'] and research['newProductionSolvedAreaM2'] == 0
bindings = {}
for entry in research['studies']:
    evidence = Path(entry['evidenceFile'])
    evidence = evidence if evidence.is_absolute() else ROOT/evidence
    assert sha(evidence) == entry['evidenceSHA256']
    assert not entry['candidateUsable'] and entry['newActualNSEvaluations'] == entry['newProductionSolvedAreaM2'] == 0
    bindings[str(evidence.relative_to(ROOT))] = sha(evidence)
row = research['studies'][-1]
diagnostic = research['patchHierarchyDiagnostics']
assert row['caseID'] == 'patchhierarchyreal01' and row['displayState'] == '粗校正放大，内部未收敛'
assert row['innerRelativeCResidual'] > 1e-3 and row['basisOrthogonalityError'] < 1e-8 and row['modelRelationRelativeDifference'] < 1e-8
assert not diagnostic['internalCTargetMet'] and not diagnostic['fullCResidualImproved']
assert diagnostic['singleVcycleFullCResidual'] > 70000
assert diagnostic['levelCellCounts'] == [217272, 29523, 4335, 732, 135, 23, 3]
assert not diagnostic['realFullArnoldiBasisSaved'] and not diagnostic['independentAllRealArnoldiColumnAuditClaimed']
audit_file = Path(diagnostic['auditFile'])
assert sha(audit_file) == diagnostic['auditSHA256'] and read(audit_file)['status'].startswith('passed-')
preflight = research['spatialReductionPreflight']
assert not preflight['numericalFactorsConstructed'] and not preflight['completeMemoryBudgetClaimed'] and not preflight['convergenceOrPhysicalQualificationProven']
assert sha(preflight['evidenceFile']) == preflight['evidenceSHA256']
payload = []
for case in data['cases']:
    for name in ('results.json', 'case-data.js'):
        relative = Path(case['resultsJSONURL']).parent/name
        before, after = VAL/'workbench23'/relative, folder/relative
        assert sha(before) == sha(after)
        payload.append(dict(path=str(after.relative_to(ROOT)), sha256=sha(after), bytes=after.stat().st_size))
assert len(payload) == 6
node = read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding_test = loadmat(VAL/'patchhierarchyui01/linear_study_binding_evidence.mat', simplify_cells=True)['report']
assert binding_test['status'].startswith('passed-') and all(binding_test['rejected'])
process = read(VAL/'patchhierarchyui01_process.json')
assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
out = dict(status='passed-independent-twelve-source-bound-research-rows-and-unchanged-facet-data', processID=os.getpid(),
    workbench='outputs/v83_validation/workbench24/index.html', researchRows=12,
    originalNSAndAcceptanceHistoryExactlyUnchanged=True, nativeAlgebraicGateAndConvergenceKeptDistinct=True,
    realFullArnoldiColumnAuditNotClaimed=True, preflightNotClaimedAsNumericalSolve=True,
    allThreeFalseSourceOrQualificationExportsRejected=True, productionQualificationNotClaimed=True,
    sourceFacetPayload=payload, previousFullMATFieldAuditReusedByByteIdentity=True,
    browserVisualOrWebGLRuntimeVerified=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256', 'sourceFacetPayload')}))
