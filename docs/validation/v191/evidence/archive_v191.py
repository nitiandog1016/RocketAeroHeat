"""本轮小证据归档；只遍历v191，不删除/清理任何原始文件。"""
from pathlib import Path
import hashlib
import json
import shutil

root=Path(__file__).resolve().parents[2]
base=root/'outputs/v191_validation';target=root/'docs/validation/v191/evidence'
assert not target.exists()
target.mkdir(parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
copied=[];local=[]
for p in sorted(base.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rel=p.relative_to(base)
    if p.name in ('fluid-0000','grid-0000'):
        local.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),backedUpOnGitHub=False))
        continue
    q=target/rel
    if q.suffix=='.log':q=q.with_suffix('.stdout.txt')
    q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    assert sha(p)==sha(q)
    copied.append(dict(source=str(p),path=q.relative_to(root).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text())
checks={p:sha(root/p)==digest for p,digest in protected.items()}
assert all(checks.values())
manifest=dict(schema='RAH.EvidenceManifest.v1',version=191,files=copied,localOnly=local,
              protectedChecks=checks,totalEvidenceBytes=sum(r['bytes'] for r in copied),
              localOutputBytes=sum(p.stat().st_size for p in base.rglob('*') if p.is_file()),
              independentDeletionCount=0,physicalQualification=False)
(target.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(files=len(copied),evidenceBytes=manifest['totalEvidenceBytes'],
                     localBytes=manifest['localOutputBytes'],protected=len(checks),localOnlyFiles=len(local))))
