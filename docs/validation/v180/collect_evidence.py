"""收集v180原始小证据；大场仅登记路径/SHA，保护并发文件，不删除任何源数据。"""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'outputs/v180_validation'
DEST = ROOT / 'docs/validation/v180/evidence'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def collect():
    assert not DEST.exists(), '原归档不得覆盖'
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks = {p: sha(ROOT / p) == expected for p, expected in protected.items()}
    assert all(checks.values()), '并发保护文件变化，不能自动接管'
    entries, local = [], []
    for source in sorted(BASE.rglob('*')):
        if not source.is_file() or '__pycache__' in source.parts or source.suffix == '.pyc':
            continue
        if source.stat().st_size > 1024**2:
            local.append(dict(path=source.relative_to(ROOT).as_posix(), bytes=source.stat().st_size, sha256=sha(source)))
            continue
        target = DEST / source.relative_to(BASE)
        if target.suffix == '.log':
            target = target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source.relative_to(ROOT).as_posix(),
                            bytes=source.stat().st_size, sha256=sha(target)))
    sources = ['runEilmerCase.m', 'scripts/prepareEilmerRocketCase.py', 'scripts/prepareEilmerDevelopmentCase.py',
               'scripts/readEilmerNative.py', 'scripts/runEilmerProcess.py', 'scripts/runPreparedEilmerCase.py',
               'scripts/runEilmerRun14Control.py', 'scripts/auditEilmerCompletion.py', 'scripts/inspectEilmerRun.py',
               'tests/testEilmerCompletion.py', 'tests/testEilmerDevelopmentEntry.m', 'tests/testEilmerDevelopmentProfile.py',
               'tests/testEilmerGroupMemory.py', 'tests/verifyEilmerDevelopmentRun.py', 'tests/testEilmerNativeInterface.py']
    result = dict(files=entries, filesCount=len(entries), totalSmallBytes=sum(x['bytes'] for x in entries),
                  largeLocalOnly=local, finalSourcesSHA256={p: sha(ROOT / p) for p in sources},
                  protectedFiles=checks, uniqueFilesDeleted=0, newQualifiedAreaM2=0,
                  note='Final sources belong to this commit; execution-time variants remain in receipts/raw snapshots. '
                       'Large files remain local, not GitHub backups. Logical sizes include hard-linked aliases.')
    (DEST.parent / 'evidence_manifest.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps(dict(files=len(entries), bytes=result['totalSmallBytes'], localOnly=len(local), protected=len(checks))))


if __name__ == '__main__':
    collect()
