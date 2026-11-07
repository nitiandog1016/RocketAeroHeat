"""只归档本轮小输入/源/证据；原生体场留本地，不打工程包或删除独有文件。"""
from pathlib import Path
import hashlib
import json
import shutil

root=Path(__file__).resolve().parents[2];base=root/'outputs/v192_validation'
target=root/'docs/validation/v192/evidence';assert not target.exists();target.mkdir(parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
copied=[];local=[];reused=[]
for p in sorted(base.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rel=p.relative_to(base)
    if p.name in ('fluid-0000','grid-0000'):
        local.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),backedUpOnGitHub=False));continue
    if rel.as_posix() in ('coarse01/mesh.su2','coarse01/surface.json'):
        previous=root/'docs/validation/v191/evidence/control01'/p.name
        assert sha(previous)==sha(p)
        reused.append(dict(source=str(p),publicPath=previous.relative_to(root).as_posix(),sha256=sha(p)));continue
    q=target/rel
    if q.suffix=='.log':q=q.with_suffix('.stdout.txt')
    q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
    copied.append(dict(source=str(p),path=q.relative_to(root).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text())
checks={p:sha(root/p)==digest for p,digest in protected.items()};assert all(checks.values())
m=dict(schema='RAH.EvidenceManifest.v1',version=192,files=copied,localOnly=local,reusedPublic=reused,
    totalEvidenceBytes=sum(r['bytes'] for r in copied),localOutputBytes=sum(p.stat().st_size for p in base.rglob('*') if p.is_file()),
    protectedChecks=checks,independentDeletionCount=0,physicalQualification=False)
(target.parent/'evidence_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print(json.dumps(dict(files=len(copied),evidenceBytes=m['totalEvidenceBytes'],localBytes=m['localOutputBytes'],
                     rawLocalPaths=len(local),reusedPublic=len(reused),protected=len(checks))))
