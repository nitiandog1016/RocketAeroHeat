"""归档小证据并保持大网格/原生数组仅本地；不删除任何独有结果。"""
from pathlib import Path
import hashlib
import json
import shutil

root=Path(__file__).resolve().parents[2];base=root/'outputs/v193_validation'
dest=root/'docs/validation/v193/evidence';assert not dest.exists();dest.mkdir(parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files=[];local=[]
for p in sorted(base.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    if p.suffix in ('.npz','.su2','.dat'):
        local.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),backedUpOnGitHub=False));continue
    rel=p.relative_to(base);q=dest/rel
    if q.suffix=='.log':q=q.with_suffix('.stdout.txt')
    q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
    files.append(dict(source=str(p),path=q.relative_to(root).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text())
checks={p:sha(root/p)==d for p,d in protected.items()};assert all(checks.values())
processes=[json.loads(p.read_text()) for p in base.glob('*_process.json')]
assert len(processes)==8 and all(p['status']=='terminal' for p in processes)
manifest=dict(schema='RAH.EvidenceManifest.v1',version=193,files=files,localOnly=local,
    totalEvidenceBytes=sum(r['bytes'] for r in files),localOutputBytes=sum(p.stat().st_size for p in base.rglob('*') if p.is_file()),
    parentSupervisedSeconds=sum(p['elapsedSeconds'] for p in processes),protectedChecks=checks,
    independentDeletionCount=0,physicalQualification=False)
(dest.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(evidenceFiles=len(files),evidenceBytes=manifest['totalEvidenceBytes'],localOnlyPaths=len(local),
                     localOutputBytes=manifest['localOutputBytes'],protectedFiles=len(checks))))
