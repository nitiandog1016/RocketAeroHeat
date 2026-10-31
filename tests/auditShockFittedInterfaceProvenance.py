"""绑定v103实际终态、显式输入和保守MATLAB调用闭包；不归档无关全项目副本。"""

from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v103_validation'
DEST = ROOT / 'docs/validation/v103'
JOBS = {
    'shockTrace103a': ['testShockFittedTraceReconstruction'],
    'shockNodal103a': ['testShockFittedNodalMotion', 'testShockFittedNodalDynamics'],
    'shockNodalCap103a': ['studyShockFittedCap'],
    'shockBoundary103a': ['testShockFittedBoundaryFlux'],
    'shockRiemannCap103a': ['studyShockFittedCap'],
    'shockNodalRiemannCap103a': ['studyShockFittedCap'],
    'shockRadialSensitivity103a': ['studyShockFittedCap'],
    'shockInterfaceView103a': ['exportShockFittedInterfaceStudy'],
}
UI = ['web/ale-study.js', 'web/ale-study.html', 'web/ale-study.css', 'web/ale-section-core.js',
      'web/ale-study-geometry.js', 'web/geometry.js', 'web/geometry.html', 'web/geometry.css']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def preserve(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == raw, str(path)
    else:
        with path.open('xb') as stream:
            stream.write(raw)


def main():
    """模块一：只接受已经结束且运行期间源稳定的实际进程，不根据旧状态重启任务。"""
    assert not (DEST / 'provenance.json').exists()
    candidates = {}
    for file in (ROOT / 'docs/validation').rglob('*.raw'):
        if '.m.' in file.name or file.name.endswith('.m.raw'):
            candidates.setdefault(sha(file), file)
    versions = {}

    def exact(relative, expected):
        key = (relative, expected)
        if key in versions:
            return (ROOT / versions[key]['exactBytesAt']).read_bytes()
        current = ROOT / relative
        file = current if current.is_file() and sha(current) == expected else candidates.get(expected)
        if file is None:
            result = subprocess.run(['git', 'show', 'HEAD:' + relative], cwd=ROOT,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == expected, relative
            raw = result.stdout
        else:
            raw = file.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == expected
        # 保存闭包中确切源码的原字节，避免Git对正常.m文件换行规范化后破坏历史SHA。
        archived = DEST / 'raw/source_versions' / (Path(relative).name + '.' + expected[:16] + '.raw')
        preserve(archived, raw)
        versions[key] = dict(restorePath=relative, sha256=expected,
                             exactBytesAt=archived.relative_to(ROOT).as_posix())
        return raw

    def closure(table, roots):
        """模块二：按工程函数名保守追踪，包括字符串形式的明确identity函数名单。"""
        names = {}
        for p, h in table.items():
            if p.endswith('.m'):
                names.setdefault(Path(p).stem, []).append((p, h))
        seen, pending, bindings = set(), list(roots), []
        while pending:
            name = pending.pop()
            if name in seen:
                continue
            seen.add(name)
            assert name in names and len(names[name]) == 1, name
            relative, expected = names[name][0]
            text = exact(relative, expected).decode('utf-8-sig')
            assert not re.search(r'\b(?:eval|evalin|feval|str2func)\s*\(', text), '动态调用须另行登记：' + relative
            references = set(re.findall(r'\b[A-Za-z]\w*\b', text)).intersection(names)
            pending.extend(sorted(references - seen))
            bindings.append(dict(path=relative, sha256=expected))
        return dict(method='conservative-project-function-name-token-closure-no-dynamic-evaluation',
                    roots=roots, bindings=bindings)

    processes, inputs = [], {}
    for job, roots in JOBS.items():
        file = VAL / (job + '_process.json')
        process = json.loads(file.read_text(encoding='utf-8-sig'))
        assert process['status'] == 'terminal' and process['exitCode'] == 0
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        preserve(DEST / 'raw/processes' / (file.name + '.raw'), file.read_bytes())
        for ending in ('.log', '_stdout.log'):
            log = VAL / (job + ending)
            preserve(DEST / 'raw/logs' / (log.name + '.raw'), log.read_bytes())
        table = process['sourceSHA256']
        graph = closure(table, roots)
        if job == 'shockInterfaceView103a':
            for relative in UI:
                exact(relative, table[relative])
        for relative in process['explicitInputs']:
            expected = table[relative]
            assert sha(ROOT / relative) == expected, relative
            inputs[relative] = expected
        processes.append(dict(job=job, status='terminal', exitCode=0,
                              engineProcessIDs=process['engineProcessIDs'],
                              processSHA256=sha(file), wholeSnapshotStable=True,
                              actualMATLABClosure=graph,
                              unrelatedWholeSnapshotDependenciesClaimed=False))
    # 模块三：独立审核与实际部署资源另外绑定，完整MAT和对象缓存不上传。
    audits = ['tests/auditShockFittedInterface.py', 'tests/auditShockFittedInterfaceExport.py',
              'tests/testShockFittedInterfaceStudy.js', 'tests/auditShockFittedInterfaceProvenance.py',
              'tests/auditShockFittedALEKernel.py', 'tests/auditShockFittedCap.py',
              'tests/auditShockFittedDynamics.py', 'tests/auditShockFittedDynamicsExport.py']
    report = dict(status='passed-terminal-processes-inputs-and-exact-conservative-source-closures',
                  processes=processes, explicitInputs=inputs, sourceVersions=list(versions.values()),
                  independentAuditorSHA256={p: sha(ROOT / p) for p in audits},
                  wholeUnrelatedProjectSnapshotsClaimedArchived=False,
                  productionQualified=False, newProductionSolvedAreaM2=0, goalRemainsActive=True)
    preserve(DEST / 'provenance.json', (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(dict(status='passed', processes=len(processes), exactSourceVersions=len(versions),
                         explicitInputs=len(inputs), productionQualified=False)))


if __name__ == '__main__':
    main()
