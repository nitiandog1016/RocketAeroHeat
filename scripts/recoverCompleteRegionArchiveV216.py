"""模块职责：只读核对已归档小证据并完成清单，保留首归档相对路径错误，不重新计算或关闭入口。"""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import time
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'outputs/v216_validation'
DOC = ROOT / 'docs/validation/v216'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(p, data):
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')


def main():
    start = time.monotonic()
    out = BASE/'archive02'
    out.mkdir(exist_ok=False)
    failure = dict(passed=False, entryExitCode=1, errorType='ValueError',
                   source='scripts/archiveCompleteRegionDecisionV216.py',
                   executedSHA256=sha(BASE/'archive01/executed_archive.py'),
                   observationSource='exec_command tool returned exit1 and traceback; not a captured process stderr file',
                   message="'scripts\\archiveCompleteRegionDecisionV216.py' does not start with 'E:\\RocketAeroHeat'",
                   scope='after evidence copies and CLI checks; manifest source list failed on relative __file__',
                   numericTasksRerun=0)
    save(BASE/'archive01/entry_failure.json', failure)
    old = json.loads((BASE/'archive01/report.json').read_text(encoding='utf8'))
    assert old['passed'] and old['CLIRefusalPassed'] and not old['CFDAdmitted']
    assert not (DOC/'evidence_manifest.json').exists()
    names = ['scripts/'+p.name for p in sorted((ROOT/'scripts').glob('*V216.py'))]
    for name in names:
        ast.parse((ROOT/name).read_text(encoding='utf8'))
    protected = json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items())
    report = dict(passed=True, originalArchiveEntryRemainsFailed=True, repeatedClosureOrNumerics=False,
                  protectedFilesUnchanged=29, seconds=time.monotonic()-start,
                  resourceScope='read-only manifest recovery; no solver peak claim')
    save(out/'report.json', report)
    files, local = [], []
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():
            continue
        item = dict(source=p.relative_to(ROOT).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
        if p.suffix=='.npz':
            local.append(item)
            continue
        q = DOC/'evidence'/p.relative_to(BASE)
        if q.exists():
            assert sha(q)==item['sha256']
        else:
            assert p.stat().st_size<1024**2
            q.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p,q)
        item['path']=q.relative_to(ROOT).as_posix()
        files.append(item)
    source_files=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md',
                  'docs/HYBRID_COMPLETE_REGION_DECISION_V216.md']+names
    source_files += [p.relative_to(ROOT).as_posix() for p in sorted(DOC.glob('*.md'))]
    assert sum(r['bytes'] for r in files)<3*1024**2
    manifest=dict(schema='RAH.V216Evidence.v1', files=files, localAssets=local,
                  sourceFiles=source_files,sourceSHA256={p:sha(ROOT/p) for p in source_files},
                  count=len(files),bytes=sum(r['bytes'] for r in files),
                  localAssetsRestore='No new remote backup. Retain local NPZ and original mother mesh; closed CLI requires a new task budget.',
                  protectedFilesUnchanged=29,geometryPartitionAccepted=True,
                  numericalWholeOperatorAccepted=False,CFDAdmitted=False,goalComplete=False)
    save(DOC/'evidence_manifest.json',manifest)
    print(json.dumps(dict(passed=True,files=len(files),bytes=manifest['bytes'],localAssets=len(local))))


if __name__ == '__main__':
    main()
