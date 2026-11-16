"""模块职责：封存真实NK控制、容量失败和科研出口，核对保护资产并生成精确发布清单。"""
from pathlib import Path
import ast
import json
import shutil
import subprocess
import sys
import time
from runSU2RocketPhaseV222 import ROOT, read, save, sha

BASE = ROOT / 'outputs/v249_validation'
DOC = ROOT / 'docs/validation/v249'
CLIS = ['scripts/' + n + '.py' for n in (
    'prepareOuterStencilV249', 'reviewNonlinearEntryV249', 'runNewtonControlV249',
    'auditNonlinearControlV249', 'exportNonlinearControlV249')]
OWN = CLIS + [
    'scripts/terminalCapacityV249.py', 'scripts/finalizeNonlinearControlV249.py',
    'native/su2PhaseAwareAir/newtonControlV249.cpp',
    'tests/testTerminalCapacityV249.py', 'tests/testNonlinearControlV249.js',
    'docs/NATIVE_NEWTON_CONTROL_V249.md', 'docs/validation/v249/PLAN.md',
    'docs/validation/v249/NATIVE_CONTROL_PLAN.md', 'docs/validation/v249/AUDIT_RECOVERY_PLAN.md',
    'docs/validation/v249/NEXT_DECISION.md',
    'AGENTS.md', 'docs/PROJECT_GOALS.md', 'README.md', '.gitattributes']
SPLIT = ('AGENTS.md', 'docs/PROJECT_GOALS.md')


def source(name):
    return BASE / 'owned_source' / name if name in SPLIT else ROOT / name


def main():
    begin = time.monotonic()
    out = BASE / 'final01'
    evidence = DOC / 'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：核对原物理资产、工程候选和并发工作，不将本地补充误公开。
    protected = read(ROOT / 'outputs/v176_validation/protected_initial.json')
    for name, expected in protected.items():
        assert sha(ROOT / name) == expected, name
    engineering = read(ROOT / 'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for item in engineering['files']:
        assert sha(ROOT / item['path']) == item['sha256']
    candidates = read(ROOT / 'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for item in candidates['files']:
        assert sha(ROOT / item['file']) == item['baselineSHA256']
        assert sha(ROOT / 'harness_tasks/engineering_bl_20261108/candidates' / item['file']) == item['candidateSHA256']
    prev = read(ROOT / 'docs/validation/v248/evidence/final_report.json')
    assets = prev['assets']
    unrelated = prev['unrelatedTrackedSHA256']
    for name, expected in dict(assets, **unrelated).items():
        assert sha(ROOT / name) == expected, name
    private = read(BASE / 'concurrent_prefixes.json')
    split = {}
    for name in SPLIT:
        live = (ROOT / name).read_text(encoding='utf8')
        prefix = private[name]['prefix']
        owned = source(name).read_text(encoding='utf8')
        assert live.count(prefix) == 1 and live.replace(prefix, '', 1) == owned
        assert 'public_benchmarks_20261116' not in owned
        split[name] = dict(localSupplementPreserved=True,
                           publishedOwnVersionSHA256=sha(source(name)), combinedWorktreeSHA256=sha(ROOT / name))

    # 模块二：分别登记原监督退出、独立容量失败和只读恢复；变化的源码绑定原执行快照。
    records = {p.parent.name: read(p) for p in BASE.glob('process_*/report.json')}
    expected_processes = {'process_' + n for n in (
        'stencil', 'method', 'prepare', 'build', 'native', 'capacity_test', 'audit', 'audit2', 'audit3', 'export')}
    assert set(records) == expected_processes
    totals, cache, bindings, required, failures = {}, {}, [], set(), {}
    for name, receipt in records.items():
        assert receipt['status'] == 'terminal' and receipt['inputsStable'] and not receipt['liveGroupAfterCleanup'], name
        if not receipt['passed']:
            failures[name] = receipt['exitCode']
        kind = receipt['budgetKind']
        totals[kind] = totals.get(kind, 0) + receipt['elapsedSeconds']
        for path, expected in receipt['inputSHA256'].items():
            if path not in cache:
                cache[path] = sha(path)
            if cache[path] == expected:
                continue
            matches = list(BASE.glob('process_*/execution_sources/' + expected + Path(path).suffix))
            matches += [p for p in BASE.glob('failed_source*/*') if p.is_file() and sha(p) == expected]
            matches += [p for p in BASE.glob('stage*/' + Path(path).name) if p.is_file() and sha(p) == expected]
            assert matches, (name, path)
            snapshot = matches[0]
            assert sha(snapshot) == expected
            required.add(snapshot)
            bindings.append(dict(process=name, path=path, executedSHA256=expected,
                                 currentSHA256=cache[path], snapshot=snapshot.relative_to(BASE).as_posix()))
    assert failures == {'process_audit': 1, 'process_audit2': 1}
    assert totals['offline'] < 120 and totals['build'] < 120 and totals['cfd'] < 180 and totals['audit'] < 120
    capacity = read(BASE / 'process_native/terminal_capacity.json')
    assert records['process_native']['passed'] and records['process_native']['exitCode'] == 0
    assert not capacity['passed'] and capacity['originalSupervisorPassed'] and capacity['originalReportPreserved']
    assert capacity['limitBytes'] == 512 * 1024**2 and capacity['logicalBytes'] > capacity['limitBytes']
    for name in ('process_capacity_test', 'process_audit3', 'process_export'):
        assert read(BASE / name / 'terminal_capacity.json')['passed'], name
    assert 'Ran 2 tests' in (BASE / 'process_capacity_test/stdout.txt').read_text(encoding='utf8')
    audit = read(BASE / 'audit03/report.json')
    ui = read(BASE / 'workbench01/report.json')
    method = read(BASE / 'method01/report.json')
    native = read(BASE / 'native01/native_completion.json')
    assert audit['passed'] and ui['passed'] and method['passed']
    assert audit['nativeControlCompleted'] and audit['completeUpdates'] == 1 and audit['rejections'] == 0
    assert native['normalReturn'] and native['updates'] == 1 and native['actualNewtonFactory'] and native['forcedNativeFinal']
    assert audit['initialElevenBitEqual'] and audit['finalNativeSixCoordinateBitEqual']
    assert not audit['originalCapacityPassed'] and not audit['physicalQualification']
    assert not audit['nonlinearConvergenceProved'] and not audit['fullJacobianAccuracyProved']
    assert ui['faces'] == 19752 and ui['updates'] == [0, 1] and ui['fields'] == 10
    assert ui['exactResearchValuesChecked'] == 158016
    assert abs(ui['originalAreaM2'] - 1.9646292565715553) < 1e-12
    frontend = read(BASE / 'workbench01/frontend_test.txt')
    assert frontend['passed'] and frontend['values'] == 158016
    assert method['PETSc']['allSourcesBlobVerified'] and not method['PETSc']['newInstallation']
    assert not method['fullResidualMeritLineSearch'] and not method['productionAdopted']
    for item in method['PETSc']['sources']:
        assert sha(item['localPath']) == item['sha256']
    closed = {}
    for name in OWN:
        if name.endswith('.py'):
            ast.parse(source(name).read_text(encoding='utf8'))
    for name in CLIS:
        proc = subprocess.run([sys.executable, str(ROOT / name), 'not-admitted'],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=10)
        message = proc.stdout.decode('utf8', errors='replace')
        assert proc.returncode != 0 and 'closed' in message, name
        closed[name] = message.strip()

    # 模块三：只发布小证据和精确最终源码；完整场/网格/外部源留E盘并登记SHA，不删除。
    out.mkdir()
    evidence.mkdir()
    local = {}

    def copy(path, rel=None):
        dest = evidence / (rel if rel is not None else path.relative_to(BASE))
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)

    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or out in path.parents or '__pycache__' in path.parts:
            continue
        rel = path.relative_to(BASE)
        top = rel.parts[0]
        keep = path in required or path.suffix in ('.json', '.txt', '.csv', '.png') or top.startswith(('process_', 'failed_source'))
        if path.suffix in ('.bin', '.exe', '.d', '.su2', '.npz', '.pyc') or top == 'owned_source':
            keep = False
        if top == 'external_sources':
            keep = path.name in ('manifest.json', 'LICENSE')
        if path.name in ('concurrent_prefixes.json', 'rocket_fas_control.json') or path.name.startswith('original_faces_'):
            keep = False
        if path.stat().st_size >= 2 * 1024**2:
            keep = False
        if path in required:
            assert keep, path
        if keep:
            copy(path)
        else:
            local[rel.as_posix()] = dict(bytes=path.stat().st_size, sha256=sha(path))
    for name in OWN:
        copy(source(name), Path('final_source') / (name if name != '.gitattributes' else '.gitattributes.snapshot'))
    save(evidence / 'local_asset_index.json', local)
    report = dict(
        passed=True, goalStatus='active', wholeSurfaceGoalComplete=False,
        originalFaces=19752, originalAreaM2=1.9646292565715553, finePoints=335726, volumeCells=705490,
        completedUpdates=1, completeSavedUpdates=[0, 1], originalFaceUpdates=ui['updates'],
        actualFactory=audit['actualFactory'], flowSpaceEvaluations=audit['flowSpaceEvaluations'],
        inferredMatrixFreeProducts=audit['inferredMatrixFreeProducts'], rejections=0,
        initialElevenBitEqual=True, finalNativeSixCoordinateBitEqual=True,
        independentSavedStates=audit['independentFineStates'],
        originalMeshIndependentStarPassed=audit['originalMeshIndependentStarPassed'],
        stencilPoints=audit['patchNodes'], stencilCells=audit['patchCells'], allSeedBoundaryCounts=0,
        originalNativeResourceContractPassed=False, originalCapacityFailure=capacity,
        originalCapacityReceiptPreserved=True, terminalCapacityTestsPassed=2,
        nativeExitedNormally=True, nativeRunRerun=False,
        workingRHSIsFinalResidual=False, nonlinearConvergenceProved=False, fullJacobianAccuracyProved=False,
        scientificImagesVisuallyInspected=True, sourceAndNumericTestsPassed=True,
        physicalQualification=False, newQualifiedAreaM2=0,
        exportedFields=ui['fields'], originalFaceValuesVerified=ui['exactResearchValuesChecked'],
        numericalObjectsRecompiled=0, sharedWebGLChanged=False, MATLABExecuted=False, browserInteractionTested=False,
        protectedCount=len(protected), engineeringFilesUnchanged=len(engineering['files']), candidateCount=len(candidates['files']),
        assets=assets, unrelatedTrackedSHA256=unrelated, splitDocumentationOwnership=split,
        concurrentLocalOnlyDataPublished=False, allReceiptsTerminalStableGroupsEmpty=True,
        executionFailuresPreserved=failures, preSupervisorFailuresPreserved=['failed_source01', 'failed_source03'],
        processSeconds=totals, maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        historicalSourceBindings=bindings, closedCLIs=closed, uniqueDataDeleted=False, newPackageInstalled=False,
        PETScVersion=method['PETSc']['version'], PETScInstalled=False,
        fullResidualMeritLineSearch=False, unmeteredDevelopmentReadsAndPublicationDisclosed=True,
        ownedFiles=OWN)
    logical_bytes = sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    report['readOnlyRecoveryOutputLimitBytes'] = 1024 * 1024**2
    report['outputLogicalBytesBeforeFinalReport'] = logical_bytes
    report['readOnlyRecoveryCapacityPassed'] = logical_bytes <= report['readOnlyRecoveryOutputLimitBytes']
    assert report['readOnlyRecoveryCapacityPassed']
    report['finalizationSeconds'] = time.monotonic() - begin
    save(out / 'report.json', report)
    copy(out / 'report.json', Path('final_report.json'))
    manifest = {p.relative_to(DOC).as_posix(): sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()}
    save(DOC / 'evidence_manifest.json', manifest)
    files = OWN + ['docs/validation/v249/evidence_manifest.json'] + [
        p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files) == len(set(files))
    save(out / 'delivery_files.json', files)
    # 最后的清单也计入只读封存容量；原512MiB报告始终保持false。
    final_bytes = sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert final_bytes <= 1024 * 1024**2
    print(json.dumps(dict(passed=True, files=len(files), evidenceFiles=len(manifest),
                         evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),
                         finalLocalBytes=final_bytes, originalNativeResourceContractPassed=False,
                         processSeconds=totals)))


if __name__ == '__main__':
    main()
