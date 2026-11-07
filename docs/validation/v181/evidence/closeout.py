"""封存实际失败、执行时源码和分项预算；不运行CFD或改动原回执。"""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[2];base=ROOT/'outputs/v181_validation';start=time.monotonic()
def read(name):return json.loads((base/name).read_text(encoding='utf8'))
def sha(data):return hashlib.sha256(data).hexdigest()
# 准入测试后只增加导出身份；按测试时SHA恢复该中间读取器，原测试不重跑。
reader=(ROOT/'scripts/readEilmerNative.py').read_bytes()
block=b"    method='Eilmer5 thermally-perfect air / Wilcox k-omega'\n    numerics={}\n    if m['condition'].get('development',{}).get('profile')=='rah-nk-real-v1':\n        method += ' / rah-nk-real-v1 real finite-difference Jv+ILU (comparison failed)'\n        numerics=dict(numericsProfile='rah-nk-real-v1',nativeNumberType='real',\n                      compatibilityStatus='small-case-comparison-failed')\n"
assert block in reader
old=reader.replace(block,b'').replace(b"int(a['outsign'][j]),method,'research',False]",b"int(a['outsign'][j]),'Eilmer5 thermally-perfect air / Wilcox k-omega','research',False]").replace(b',**numerics)',b')')
expected=read('admission01/report.json')['sourceSHA256']['scripts/readEilmerNative.py'];assert sha(old)==expected
target=base/'readEilmerNative-at-admission.py.raw';assert not target.exists();target.write_bytes(old)
run=read('small01_run/report.json');prep=read('small01_prepare/report.json')
for name in ['scripts/readEilmerNative.py','scripts/runPreparedEilmerCase.py','scripts/auditEilmerCompletion.py']:
    locked=run['commands'][0]['inputSHA256'][str(ROOT/name)]
    assert sha((base/'source_at_small_run'/(Path(name).name+'.raw')).read_bytes())==locked
protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
checks={p:sha((ROOT/Path(p.replace('\\','/'))).read_bytes())==h for p,h in protected.items()};assert all(checks.values())
rows=[]
for name in ['source01/report.json','profile01/report.json','small01_prepare/report.json','small01_run/report.json',
             'small01_check/report.json','difference01/report.json','admission01/report.json','identity01/report.json']:
    r=read(name);rows.append(dict(receipt='evidence/'+name,passed=r['passed'],seconds=r['elapsedSeconds']))
budget=dict(entries=rows,measuredEntrySumSeconds=sum(x['seconds'] for x in rows),
            nestedSmallProcessSeconds=run['commands'][0]['elapsedSeconds'],
            smallOrchestratorToolWallSeconds=5.2612614,scope='Entries are measured program time; nested process/tool and reading/development time are not added as a total.',
            firstFieldReadFailure="KeyError: tke; actual fields tq-tke/tq-omega; corrected read-only, no CFD rerun",
            goalDocumentWriteFailures=['Windows pipeline Python SyntaxError Non-UTF-8; no file written',
                                       'Windows Python 3.7.7 also rejected UTF-8 source helper; fixed-runtime WSL Python completed same helper'],
            nativeSmallRuns=1,nativeFullDomainRuns=0,builds=0,matlabRuns=0,newMainUIRuns=0,uniqueFilesDeleted=0,
            limits=dict(initialSourceAndAnalyticSeconds=120,smallPrepareSeconds=60,smallRunSeconds=30,smallReadSeconds=20,smallParentSeconds=120),
            unknownFullDomainMemoryFit=True,physicalQualification=False)
(ROOT/'docs/validation/v181/budget_receipt.json').write_text(json.dumps(budget,indent=2)+'\n')
result=dict(passed=True,actualSmallComparisonPassed=False,realFullDomainEnabled=False,protectedFiles=checks,
            nativeExecutionSourcesBound=True,admissionSourceRestoredSHA256=expected,allOwnedProcessesTerminal=not run['commands'][0]['liveGroupAfterCleanup'],
            elapsedSeconds=time.monotonic()-start,newQualifiedAreaM2=0)
(base/'closeout01.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['passed','actualSmallComparisonPassed','realFullDomainEnabled','allOwnedProcessesTerminal','elapsedSeconds']}))
