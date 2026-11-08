"""模块职责：绑定已结束离线任务、关闭入口、保护历史并归档小证据，不重做数值/网格任务。"""
from pathlib import Path
import ast
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'outputs/v216_validation'
DOC = ROOT / 'docs/validation/v216'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf8')


def main():
    started = time.monotonic()
    out = BASE / 'archive01'
    out.mkdir(exist_ok=False)
    names = ['inspectCompleteRegionInputsV216.py', 'designCompleteRegionsV216.py',
             'assessFrozenDifferenceResolutionV216.py', 'closeCompleteRegionTopologyV216.py',
             'finalizeCompleteRegionDecisionV216.py', 'plotCompleteRegionsV216.py',
             'runCompleteRegionDecisionV216.py']
    # 模块一：只审核已结束的监督回执，记录早期未保存文本的明确缺口。
    reports = {}
    missing = []
    for path in sorted(BASE.glob('process_*/report.json')):
        r = json.loads(path.read_text(encoding='utf8'))
        assert r['status'] != 'running'
        assert r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None
        expected = path.parent.name not in ('process_inventory', 'process_inventory2', 'process_design', 'process_design2')
        assert r['passed'] == expected and r['nativeExitCode'] == (0 if expected else 1)
        for p, h in r['inputSHA256'].items():
            name = p.replace('\\', '/').split('/')[-1]
            if name not in names:
                continue
            local = path.parent / ('executed_entry.py' if name == 'runCompleteRegionDecisionV216.py' else 'executed_task.py')
            if local.exists():
                assert sha(local) == h, (path.parent.name, name)
            else:
                missing.append(dict(stage=path.parent.name, sourceName=name, recordedSHA256=h))
        reports[path.parent.name] = {k: r[k] for k in ('passed', 'nativeExitCode', 'elapsedSeconds',
                                                       'peakObservedGroupRSSBytes', 'inputsStable')}
    assert len(reports) == 11 and len(missing) == 3
    elapsed = sum(r['elapsedSeconds'] for r in reports.values())
    assert elapsed < 120
    # 绘图实际耗时符合15秒，但历史监督器配置仍为通用90秒，不追改回执。
    assert reports['process_plot']['elapsedSeconds'] < 15
    final = json.loads((BASE / 'final01/report.json').read_text())
    assert final['passed'] and final['allRegionsConnected'] and not final['numericalWholeOperatorAccepted']
    for p, h in final['sourceSHA256'].items():
        assert sha(ROOT / p.replace('\\', '/')) == h
    assert sha(BASE / 'final01/decision.json') == final['decisionSHA256']
    plot = json.loads((BASE / 'plot01/report.json').read_text())
    assert plot['passed'] and sha(BASE / 'plot01/regions.png') == plot['figureSHA256']
    for p, h in plot['sourceSHA256'].items():
        assert sha(ROOT / p.replace('\\', '/')) == h
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected) == 29 and all(sha(ROOT / p.replace('\\', '/')) == h for p, h in protected.items())

    # 模块二：存准确执行版，再关闭CLI；导入的函数仍可供后续明确预算的新任务复用。
    closed = BASE / 'closed_entries'
    closed.mkdir(exist_ok=False)
    for name in names:
        source = ROOT / 'scripts' / name
        shutil.copyfile(source, closed / name)
        text = source.read_text(encoding='utf8')
        guard = re.search(r"^if __name__\s*==\s*['\"]__main__['\"]\s*:", text, re.M)
        assert guard
        source.write_text(text[:guard.start()] +
                          "if __name__ == '__main__':\n    raise ValueError('v216离线设计已经结束；请按NEXT_DECISION登记新任务，当前无构建/查询/网格/CFD准入。')\n",
                          encoding='utf8')
        ast.parse(source.read_text(encoding='utf8'))
    result = subprocess.run([sys.executable, str(ROOT / 'scripts/runCompleteRegionDecisionV216.py'), 'final'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5)
    assert result.returncode != 0 and b'ValueError' in result.stderr
    (out / 'closed_entry_stderr.txt').write_bytes(result.stderr)
    report = dict(passed=True, processes=reports, offlineSupervisedSeconds=elapsed,
                  offlinePeakObservedRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in reports.values()),
                  missingImmediateSourceSnapshots=missing, allLaterExecutedSnapshotsBound=True,
                  initialWSLGenericCPUField=2, initialWSLActualTasksetCPUs=[0],
                  plotPlannedSeconds=15, plotSupervisorConfiguredSeconds=90, plotActualWithin15Seconds=True,
                  CLIClosed=names, CLIRefusalPassed=True, protectedFilesUnchanged=29,
                  geometryPartitionAccepted=True, numericalWholeOperatorAccepted=False, CFDAdmitted=False,
                  CFDAdvanceRuns=0, newQualifiedAreaM2=0, goalComplete=False, uniqueDeleted=0,
                  bookkeepingSeconds=time.monotonic()-started,
                  bookkeepingResourceScope='archive/hash/closure only; not solver peak or numerical recomputation')
    save(out / 'report.json', report)

    # 模块三：小证据逐字公开，完整索引不复制到Git；清单登记本地大小/摘要与恢复限制。
    dest = DOC / 'evidence'
    dest.mkdir(exist_ok=False)
    files, local = [], []
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(BASE)
        item = dict(source=p.relative_to(ROOT).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
        if p.suffix == '.npz':
            local.append(item)
            continue
        assert p.stat().st_size < 1024**2
        q = dest / rel
        q.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, q)
        assert sha(q) == item['sha256']
        item['path'] = q.relative_to(ROOT).as_posix()
        files.append(item)
    source_files = ['AGENTS.md', '.gitattributes', 'docs/PROJECT_GOALS.md',
                    'docs/HYBRID_COMPLETE_REGION_DECISION_V216.md']
    source_files += ['scripts/' + name for name in names]
    source_files.append(Path(__file__).relative_to(ROOT).as_posix())
    source_files += [p.relative_to(ROOT).as_posix() for p in sorted(DOC.glob('*.md'))]
    assert sum(r['bytes'] for r in files) < 3*1024**2
    manifest = dict(schema='RAH.V216Evidence.v1', files=files, localAssets=local,
                    sourceFiles=source_files, sourceSHA256={p:sha(ROOT/p) for p in source_files},
                    count=len(files), bytes=sum(r['bytes'] for r in files),
                    localAssetsRestore='No new remote backup. Retain local NPZ and original mother mesh; do not rerun closed CLI without a new budget.',
                    protectedFilesUnchanged=29, geometryPartitionAccepted=True,
                    numericalWholeOperatorAccepted=False, CFDAdmitted=False, goalComplete=False)
    save(DOC / 'evidence_manifest.json', manifest)
    print(json.dumps(dict(passed=True, files=len(files), bytes=manifest['bytes'],
                         offlineSeconds=elapsed, missingEarlySnapshots=len(missing))))


if __name__ == '__main__':
    main()
