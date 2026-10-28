# -*- coding: utf-8 -*-
"""登记真实尺度研究及工作台证据；资源停止不算数值通过，不增加物理解面积。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
target = ROOT / 'docs/validation/v83/progress_equilibration.json'
assert not target.exists()

def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            result.update(block)
    return result.hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def identity(path):
    path = Path(path)
    return dict(path=str(path.relative_to(ROOT)), sha256=sha(path), bytes=path.stat().st_size)

names = ['schurscaleunit01', 'schurscaleunit02', 'schurscaleunit03', 'schurscaleunit04', 'schurscalereal01', 'schurscaleui01']
processes = {}
for name in names:
    p = read(VAL/(name+'_process.json'))
    assert p['status'] == 'terminal' and all(p['sourceStableAfter'].values())
    expected = 0 if name in ('schurscaleunit04','schurscaleui01') else 4294967295 if name=='schurscalereal01' else 1
    assert p['exitCode'] == expected
    processes[name] = dict(exitCode=p['exitCode'], sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=identity(VAL/(name+'_process.json')))

small = read(VAL/'schurscaleunit04_independent_pivoted_schur.json')
real = read(VAL/'schurscalereal01_independent_equilibration.json')
assert small['status'].startswith('passed-') and small['equilibrationIndependentlyVerified']
assert real['status'].startswith('passed-') and not real['coefficientsAddedOrDeleted']
stop = read(VAL/'schurscalereal01/agent_stop_record.json')
assert not stop['numericalFailureOrZeroPivotClaimed'] and not stop['usableFactorProduced']
for name in ('pivoted_schur_factor_evidence.mat','pivoted_schur_real_summary.json','pivoted_field_schur_actions.mat'):
    assert not (VAL/'schurscalereal01'/name).exists()

view = VAL/'workbench20'
data = read(view/'wall-velocity-study-data.json')
research = read(VAL/'block_linear_research_summary_v4.json')
assert data['linearResearch'] == research and len(research['studies']) == 7
new_row = research['studies'][-1]
assert new_row['displayState'] == '资源停止（数值未定）' and not new_row['candidateUsable']
for key in ('innerRelativeCResidual','basisOrthogonalityError','modelRelationRelativeDifference'):
    assert new_row[key] is None
node = read(view/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15
assert not node['browserVisualOrWebGLRuntimeVerified']
unchanged = []
for c in data['cases']:
    original = Path(c['resultsJSONURL'])
    for rel in (original, original.parent/'case-data.js'):
        a = view/rel
        b = VAL/'workbench19'/rel
        assert sha(a) == sha(b)
        unchanged.append(identity(a))
assert len(unchanged)==6
ui = dict(status='passed-independent-source-bindings-seven-research-rows-unknown-metrics-and-unchanged-facet-data',
          processID=os.getpid(), workbench='outputs/v83_validation/workbench20/index.html',
          sourceFacetPayload=unchanged, previousFullMATFieldAuditReusedByByteIdentity=True,
          numericalAcceptanceOfStoppedAttemptClaimed=False, browserVisualOrWebGLRuntimeVerified=False)
(view/'independent_equilibration_study_view.json').write_text(json.dumps(ui,indent=2)+'\n',encoding='utf8')

old = read(ROOT/'docs/validation/v83/progress_schur_pivot.json')
production = old['production']
for key in ('latestEndpoint','originalNSSource'):
    assert sha(ROOT/production[key]['path']) == production[key]['sha256']
assert subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip()=='42601c9'
files = [ROOT/'internalFlow/cutCell'/name for name in (
    'equilibrateRocketSchurScales.m','factorRocketPivotedSchur.m','solveRocketPivotedSchur.m')]
files += [ROOT/'tests'/name for name in (
    'testRocketSchurEquilibration.m','studyRocketPivotedSchurActions.m','auditRocketPivotedSchur.py',
    'auditRocketRealPivotedSchur.py','auditRocketSchurEquilibration.py','exportRocketSchurEquilibrationDiagnostics.py')]
out = dict(status='tested-reversible-equilibration-with-independent-all-real-coefficients-audit-native-factorization-resource-stopped',
    recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),version=83,processes=processes,
    smallDomainAudit=small,completeRealScaleAudit=real,realStop=stop,
    realNativeFactorizationCompleted=False,realNumericalFailureOrZeroPivotConcluded=False,
    fullFiveFieldRealPreconditionerAvailable=False,threeRealRhsOrGMRESStarted=False,
    repairHistory=dict(failedUnitAttemptsRetained=names[:3],
        changes=['Use actual raw CSR for bit-exact scales; keep independent S reference comparison.',
                 'Rename local full variable to avoid shadowing full().',
                 'Retain failed tridiagonal stress; use an explicitly well-conditioned balanced known-solution case for the tight accuracy contract. No physical gate changed.']),
    userInterface=ui,nodeDataVerification=node,production=production,delivery=old['delivery'],
    codeSourceIdentity=[identity(p) for p in files],
    draftVerificationEntry=dict(path='tests/auditRocketRealPivotedSchur.py',syntaxChecked=True,numericalRun=False,
        reason='No real pivoted factor exists. Source is drafted, not numerical passing evidence.'),
    workAccounting=dict(newSelectedJApplications=0,newActualNSEvaluations=0,newAcceptedPhysicalStates=0,newProductionSolvedAreaM2=0),
    nextAction='Do not restart the stopped native ILUTP or repeat raw S assembly/scale construction/threshold scans. Use the saved equilibrated matrix to inspect weak modes, approximate Schur structure, and sparse fill/resource feasibility before a new controlled solver choice. Keep all physical CVs, original C/NS, direction/temperature/descent/convergence gates and co-priority regional validation/UI.',
    githubReadiness=dict(conditionalReviewVersion=90,versionIsNotACommitment=True,regionMilestoneMet=False,commitOrPushPerformed=False),
    goalStatus='active',wholeGoalComplete=False,recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],processID=os.getpid(),manifest=identity(target),sourceFacetFilesUnchanged=len(unchanged))))
