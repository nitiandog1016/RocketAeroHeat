"""绑定v104实际终态、输入与准确源码版本；复用已有结果，不重启MATLAB。"""

from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v104_validation'
DEST = ROOT / 'docs/validation/v104'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def sha(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def preserve(file, raw):
    file.parent.mkdir(parents=True, exist_ok=True)
    if file.exists():
        assert file.read_bytes() == raw, str(file)
    else:
        with file.open('xb') as stream:
            stream.write(raw)


def main():
    # 模块一：仅按记录SHA寻找确切原字节；历史失败也保留，不用当前源码替代。
    candidates = {}
    for file in (ROOT / 'docs/validation').rglob('*.raw'):
        if '.m.' in file.name or file.name.endswith('.m.raw'):
            candidates.setdefault(sha(file), file)
    versions = {}
    git_cache = {}

    def exact(relative, expected):
        key = (relative, expected)
        if key in versions:
            return (ROOT / versions[key]['exactBytesAt']).read_bytes()
        file = ROOT / relative
        raw = file.read_bytes() if file.is_file() else b''
        if digest(raw) != expected:
            raw = candidates[expected].read_bytes() if expected in candidates else b''
        if digest(raw) != expected:
            if relative not in git_cache:
                trial = subprocess.run(['git', 'show', '0b5e285:' + relative], cwd=ROOT,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                git_cache[relative] = trial.stdout if trial.returncode == 0 else b''
            old = git_cache[relative]
            for variant in (old, old.replace(b'\r\n', b'\n'),
                            old.replace(b'\r\n', b'\n').replace(b'\n', b'\r\n')):
                if digest(variant) == expected:
                    raw = variant
                    break
        if digest(raw) != expected and relative == 'tests/studyEngineeringWholeSurface.m':
            # 首次失败只涉及测试入口的整型角度；重建候选必须与当时完整SHA相同。
            text = file.read_text(encoding='utf-8-sig')
            text = text.replace('% 小型输入为节省空间可存整数；GUI/计算入口的角度统一转为double。\n', '')
            text = re.sub(r'double\((loaded\.fs\.\w+)\)', r'\1', text)
            for variant in (text.encode('utf-8'), text.replace('\n', '\r\n').encode('utf-8')):
                if digest(variant) == expected:
                    raw = variant
                    break
        assert digest(raw) == expected, '缺少准确源码版本: ' + relative + ' ' + expected
        target = DEST / 'source_versions' / (Path(relative).name + '.' + expected[:16] + '.raw')
        preserve(target, raw)
        versions[key] = dict(restorePath=relative, sha256=expected,
                             exactBytesAt=target.relative_to(ROOT).as_posix())
        return raw

    def closure(table, expression):
        # 模块二：保守的工程函数名词元闭包包含注释和evalc字符串，不冒充实际执行轨迹。
        names = {}
        for relative, expected in table.items():
            if relative.endswith('.m'):
                names.setdefault(Path(relative).stem, []).append((relative, expected))
        roots = sorted(set(re.findall(r'\b[A-Za-z]\w*\b', expression)).intersection(names))
        seen, pending, bindings, dynamic = set(), list(roots), [], []
        while pending:
            name = pending.pop()
            if name in seen:
                continue
            seen.add(name)
            assert len(names[name]) == 1, name
            relative, expected = names[name][0]
            text = exact(relative, expected).decode('utf-8-sig')
            calls = sorted(set(re.findall(r'\b(eval|evalc|evalin|feval|str2func)\s*\(', text)))
            if calls:
                dynamic.append(dict(path=relative, functions=calls,
                                    scope='literal-project-name-tokens-included; not-a-runtime-call-trace'))
            references = set(re.findall(r'\b[A-Za-z]\w*\b', text)).intersection(names)
            pending.extend(sorted(references - seen))
            bindings.append(dict(path=relative, sha256=expected))
        return dict(method='conservative-project-function-name-token-closure-including-strings',
                    roots=roots, bindings=bindings, dynamicCallSites=dynamic,
                    actualExecutedCallTraceClaimed=False)

    processes, inputs = [], {}
    for file in sorted(VAL.glob('*_process.json')):
        job = file.name[:-len('_process.json')]
        process = json.loads(file.read_text(encoding='utf-8-sig'))
        assert process['status'] == 'terminal'
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        preserve(DEST / 'processes' / (file.name + '.raw'), file.read_bytes())
        for ending in ('.log', '_stdout.log'):
            log = VAL / (job + ending)
            preserve(DEST / 'logs' / (log.name + '.raw'), log.read_bytes())
        graph = closure(process['sourceSHA256'], process['expression'])
        for relative in process['explicitInputs']:
            expected = process['sourceSHA256'][relative]
            assert sha(ROOT / relative) == expected, relative
            inputs[relative] = dict(sha256=expected, bytes=(ROOT / relative).stat().st_size,
                                   localOnly=relative.startswith('outputs/'))
        processes.append(dict(job=job, exitCode=process['exitCode'], status='terminal',
                              expression=process['expression'], processSHA256=sha(file),
                              engineProcessIDs=process['engineProcessIDs'],
                              wholeSnapshotStable=True, conservativeSourceClosure=graph))
    assert len(processes) == 13
    assert sum(p['exitCode'] == 0 for p in processes) == 9
    assert sum(p['exitCode'] == 1 for p in processes) == 3
    assert [p['job'] for p in processes if p['exitCode'] == 4294967295] == ['kinney_regression01']
    interruption = VAL / 'kinney_regression01_interruption.json'
    preserve(DEST / 'processes' / (interruption.name + '.raw'), interruption.read_bytes())

    # 模块三：复制小型摘要、修复来源；大型MAT和文献只记录位置/SHA及未上传身份。
    summaries = {}
    for folder in ('baseline02', 'tests03', 'actual03', 'workflow02', 'provenance01'):
        file = VAL / folder / 'summary.json'
        preserve(DEST / 'summaries' / (folder + '.json'), file.read_bytes())
        summaries[folder] = sha(file)
    small = VAL / 'provenance01/provenance.mat'
    preserve(DEST / 'baseline_provenance.mat', small.read_bytes())
    pdf = VAL / 'research/kinney_cbaero_2007.pdf'
    inputs[pdf.relative_to(ROOT).as_posix()] = dict(sha256=sha(pdf), bytes=pdf.stat().st_size,
        localOnly=True, url='https://ntrs.nasa.gov/citations/20070017947', inspectedPrintedPages=[4, 5])
    report = dict(status='passed-source-binding-not-all-processes-passed',
                  processes=processes, explicitInputs=inputs, sourceVersions=list(versions.values()),
                  summarySHA256=summaries, baselineProvenanceSHA256=sha(small),
                  independentAuditorSHA256={p: sha(ROOT / p) for p in
                      ('tests/auditKinneyAcreage.py', 'tests/testKinneyViewer.js', 'tests/auditKinneyProvenance.py')},
                  wholeUnrelatedProjectSnapshotsClaimedArchived=False,
                  productionQualified=False, newProductionSolvedAreaM2=0, goalRemainsActive=True)
    preserve(DEST / 'provenance.json', (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(dict(status=report['status'], processes=len(processes), exit0=9, exit1=3,
                         interrupted=1, sourceVersions=len(versions), inputs=len(inputs))))


if __name__ == '__main__':
    main()
