"""归档本轮小证据和准确执行源，大场、二进制及其他协作者的WebGL副本留在本地。"""
from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import sha

BASE = ROOT/'outputs/v163_validation'
DEST = ROOT/'docs/validation/v163/evidence'
assert not DEST.exists()
selected = [p for p in BASE.iterdir() if p.is_file()]
for name in ('fixture01', 'initial_source', 'small01', 'small02', 'real01', 'main01', 'matlab'):
    selected += [p for p in (BASE/name).rglob('*') if p.is_file() and p.suffix in
                 ('.json', '.py', '.cpp', '.hpp', '.md', '.log', '.cfg', '.csv', '.vtk')]
for name in ('build01', 'build02', 'build03'):
    folder = BASE/name
    selected += [p for p in folder.iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.py', '.diff')]
    if name == 'build03':
        selected += [p for p in (folder/'executed_sources').rglob('*') if p.is_file()]
    else:
        selected.append(folder/'executed_sources/SU2_CFD/src/SU2_CFD.cpp')
folder = BASE/'real02'
selected += [p for p in folder.iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.cfg')]
selected.append(folder/'history.csv')
selected += [p for p in (folder/'executed_sources').rglob('*') if p.is_file()]
selected += [p for p in (folder/'audit0').iterdir() if p.suffix in ('.json', '.log')]
selected += [p for p in (BASE/'matlab_sources').rglob('*.m') if p.is_file()]
entries = []
for source in sorted(set(selected)):
    assert source.stat().st_size < 2*1024**2, str(source)
    target = DEST/source.relative_to(ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(source), str(target))
    assert sha(source) == sha(target)
    entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source.relative_to(ROOT).as_posix(),
        bytes=source.stat().st_size, sha256=sha(source)))
report = dict(count=len(entries), bytes=sum(row['bytes'] for row in entries), files=entries,
    localOnly=['real full fields', 'large oracle CSV', 'native runtimes', 'existing collaborator WebGL snapshot'],
    uniqueResultsDeleted=0, cleanupBytes=0, physicsQualified=False, newQualifiedArea=0)
(DEST.parent/'evidence_manifest.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
print(dict(count=report['count'], bytes=report['bytes']))
