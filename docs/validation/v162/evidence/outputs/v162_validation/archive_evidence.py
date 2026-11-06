"""归档实际监督、文本夹具和准确源码，网格/初场/运行时按散列保留本地。"""
from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import sha

BASE = ROOT/'outputs/v162_validation'
DEST = ROOT/'docs/validation/v162/evidence'
assert not DEST.exists()
selected = [p for p in BASE.iterdir() if p.is_file()]
for name in ('process_test01', 'failed_precision_source', 'final_source'):
    selected += [p for p in (BASE/name).rglob('*') if p.is_file()]
for name in ('debug01', 'debug02'):
    folder = BASE/name
    selected += [p for p in folder.iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.cfg', '.gdb')]
    if (folder/'history.csv').exists():
        selected.append(folder/'history.csv')
    selected += [p for p in (folder/'executed_sources').rglob('*') if p.is_file()]
entries = []
for source in sorted(set(selected)):
    target = DEST/source.relative_to(ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(source), str(target))
    assert source.read_bytes() == target.read_bytes()
    entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source.relative_to(ROOT).as_posix(),
        bytes=source.stat().st_size, sha256=sha(source)))
report = dict(files=entries, count=len(entries), bytes=sum(row['bytes'] for row in entries),
    localOnly=['mesh', 'initial fields', 'native executables'], nativeBuilds=0,
    newCompleteFlowStates=0, uniqueResultsDeleted=0, physicsQualified=False)
(DEST.parent/'evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n', encoding='utf8')
print(dict(count=report['count'], bytes=report['bytes']))
