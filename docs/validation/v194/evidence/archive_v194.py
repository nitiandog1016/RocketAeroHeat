"""保存小证据与大文件位置清单；不删除/移动独有数据，不复制共享历史。"""
from pathlib import Path
import hashlib
import json
import shutil

root=Path(__file__).resolve().parents[2];base=root/'outputs/v194_validation'
dest=root/'docs/validation/v194/evidence';assert not dest.exists();dest.mkdir(parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files=[];local=[]
for p in sorted(base.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:continue
    if p.suffix in ('.npz','.su2') or p.name in ('grid-0000','fluid-0000'):
        local.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),backedUpOnGitHub=False));continue
    q=dest/p.relative_to(base)
    if q.suffix=='.log':q=q.with_suffix('.stdout.txt')
    q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
    files.append(dict(source=str(p),path=q.relative_to(root).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
for name in ['scripts/optimizeHybridRocketEulerMesh.py','scripts/recoverHybridRocketNetgen.py',
             'scripts/prepareHybridRocketEulerNative.py','tests/testHybridRocketEulerPreparation.py',
             'tests/verifyHybridNetgenRecovery.py']:
    p=root/name;q=dest/'final_source'/name;q.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(p,q);assert sha(p)==sha(q)
    files.append(dict(source=str(p),path=q.relative_to(root).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text())
checks={p:sha(root/p)==d for p,d in protected.items()};assert len(checks)==29 and all(checks.values())
processes=[json.loads(p.read_text()) for p in base.glob('*_process.json')]
assert len(processes)==11 and all(p['status']=='terminal' and p['stoppedByBudget'] is None for p in processes)
manifest=dict(schema='RAH.EvidenceManifest.v1',version=194,files=files,localOnly=local,
    totalEvidenceBytes=sum(r['bytes'] for r in files),localOutputBytes=sum(p.stat().st_size for p in base.rglob('*') if p.is_file()),
    parentSupervisedSeconds=sum(p['elapsedSeconds'] for p in processes),parentSupervisedProcesses=len(processes),
    protectedChecks=checks,independentDeletionCount=0,physicalQualification=False,
    CFDEvolutionRuns=0,nativeOptimizerCalls=1,nativePreparationSequences=1)
assert manifest['localOutputBytes']<80*1024**2 and manifest['totalEvidenceBytes']<20*1024**2
(dest.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(files=len(files),evidenceBytes=manifest['totalEvidenceBytes'],
    localOnlyFiles=len(local),localBytes=manifest['localOutputBytes'],parents=len(processes),seconds=manifest['parentSupervisedSeconds'])))
