"""归档本轮完整域的小证据、导出和准确源码，独有体场仅本地保存。"""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'outputs/v160_validation'
DEST = ROOT/'docs/validation/v160/evidence'
assert not DEST.exists()
selected = [p for p in SOURCE.iterdir() if p.is_file()]
for folder in ('field_review01', 'field_review02', 'matlab', 'pre_matlab_source', 'failed_main01_source', 'final_source'):
    selected += [p for p in (SOURCE/folder).rglob('*') if p.is_file()]
selected += [p for p in (SOURCE/'main02').iterdir() if p.suffix in ('.json', '.csv')]
case = SOURCE/'real01'
selected += [p for p in case.iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.cfg')]
selected.append(case/'history.csv')
selected += [p for p in (case/'executed_sources').rglob('*') if p.is_file()]
selected += [p for p in (case/'ideal_air_audit').iterdir() if p.is_file() and p.suffix in ('.json', '.log', '.py')]
entries = []
for source in sorted(set(selected)):
    target = DEST/source.relative_to(ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(source), str(target))
    raw = source.read_bytes()
    assert raw == target.read_bytes()
    entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source.relative_to(ROOT).as_posix(),
                        bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
manifest = dict(files=entries, count=len(entries), bytes=sum(row['bytes'] for row in entries),
    localOnly=['volume and surface native fields', 'full mesh and initial state', 'full independent thermodynamic tables', 'result MAT'],
    physicalSolutionQualified=False, uniqueResultsDeleted=0)
(DEST.parent/'evidence_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
print(dict(count=manifest['count'], bytes=manifest['bytes']))
