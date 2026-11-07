"""仅归档本轮小证据并核对保护项，不复制场文件或删除任何数据。"""
from pathlib import Path
import hashlib
import json
import shutil
import time

ROOT=Path('/mnt/e/RocketAeroHeat');base=ROOT/'outputs/v188_validation'
target=ROOT/'docs/validation/v188';start=time.monotonic()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,o):p.write_text(json.dumps(o,indent=2)+'\n')
protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
for name,h in protected.items():assert sha(ROOT/name.replace('\\','/'))==h,name
resource=json.loads((base/'resource01/report.json').read_text())
test=json.loads((base/'test02/report.json').read_text())
assert resource['passed'] and test['passed'] and len(resource['rows'])==41
assert all(sha(p)==h for p,h in resource['inputSHA256'].items())
assert all(sha(ROOT/p)==h for p,h in test['sourceSHA256'].items())
files=['resource_source01.py','reader_source01.py','executed_test01.py','executed_test02.py',
       'resource_parent.json','resource_parent.stdout.txt','test_parent.json','test_parent.stdout.txt',
       'test02_parent.json','test02_parent.stdout.txt','resource01/report.json','resource01/source_excerpts.json',
       'test01/report.json','test02/report.json','deliver.py']
entries=[]
for name in files:
    src=base/name;dst=target/'evidence'/name;dst.parent.mkdir(parents=True,exist_ok=True)
    assert not dst.exists();shutil.copyfile(src,dst)
    assert sha(src)==sha(dst)
    entries.append(dict(path=dst.relative_to(ROOT).as_posix(),source=src.relative_to(ROOT).as_posix(),
                        bytes=dst.stat().st_size,sha256=sha(dst)))
parents=[json.loads((base/n).read_text()) for n in ['resource_parent.json','test_parent.json','test02_parent.json']]
assert [p['exitCode'] for p in parents]==[0,1,0]
write(target/'budget_receipt.json',dict(budgetSeconds=120,parentTotalSeconds=sum(p['elapsedSeconds'] for p in parents),
    phases=[dict(name=n,**p) for n,p in zip(['resource','failed-test','recovered-test'],parents)],
    sourceReadingEditingCoordinationSeparate=True,nativeRuns=0,builds=0,
    preliminaryInheritedArithmeticSeconds=None,preliminaryArithmeticNotFormalEvidence=True))
source_names=['.gitattributes','AGENTS.md','docs/PROJECT_GOALS.md','scripts/eilmerResourceBudget.py',
              'scripts/readEilmerNative.py','tests/testEilmerResourceBudget.py','docs/EILMER_RESOURCE_FEASIBILITY_V188.md',
              'docs/validation/v188/PLAN.md','docs/validation/v188/RESOURCE_ACCOUNTING.md',
              'docs/validation/v188/NEXT_DECISION.md','docs/validation/v188/README.md']
write(target/'evidence_manifest.json',dict(files=entries,finalSourceSHA256={n:sha(ROOT/n) for n in source_names},
    previousCommit='1488a3abc714967c50a3b5b80f6398c1998ebd33',physicalQualification=False,
    protectedFilesUnchanged=len(protected),resourceInputsStable=True,testSourcesStable=True,
    uniqueFilesDeleted=0,systemResourceConfigurationChanged=False,
    knownLargeData='Existing fields, meshes, runtimes and objects remain local; no new CFD or binary.',
    elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(evidenceFiles=len(entries),evidenceBytes=sum(e['bytes'] for e in entries),
    protectedFilesUnchanged=len(protected),boundBytes=resource['minimumAllocationLowerBoundBytes'],
    sourceInputCount=len(resource['inputSHA256']),elapsedSeconds=time.monotonic()-start)))
