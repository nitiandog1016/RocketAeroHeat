"""收集v185真实构建、小例、独立比较与原失败；大构建物只登记，不删除独有数据。"""
from pathlib import Path
import hashlib
import json
import shutil
import time

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'outputs/v185_validation'
DEST = ROOT / 'docs/validation/v185/evidence'
PACKAGE = ROOT / 'tools/runtime/eilmer-explicit-dense-v1'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(name):
    return json.loads((BASE/name).read_text('utf-8-sig'))


def collect():
    """模块一：以准确结果与执行源为准核对终态，不重跑原生或旧控制。"""
    start = time.monotonic(); assert not DEST.exists()
    protected = json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks = {p:sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items()}
    assert len(checks)==29 and all(checks.values())
    scope = read('scope02/report.json'); build = read('build_result.json')
    prep = read('small01_prepare/report.json'); run = read('small01_run/report.json')
    verify = read('verify02/report.json'); first = read('verify01/report.json')
    assert scope['passed'] and build['passed'] and prep['passed'] and run['passed'] and verify['passed']
    assert first['passed'] is False and read('scope01/failure.json')['passed'] is False
    assert all(sha(p)==h for p,h in scope['originalSHA256'].items()), '共享原源/对象/安装发生变化'
    for name,digest in verify['sourceSHA256'].items(): assert sha(ROOT/name)==digest
    for name,digest in run['outputSHA256'].items(): assert sha(name)==digest
    assert all(verify['fieldAndWallBytesIdentical'].values()) and verify['maxScaledError']==0
    assert len(verify['rejected'])==12 and verify['comparedScalars']==658
    assert not build['process']['liveGroupAfterCleanup'] and not run['commands'][0]['liveGroupAfterCleanup']
    parents = {name:read(name+'.json') for name in ['build_parent','prepare_parent','run_parent','verify_parent','verify02_parent']}
    assert all(p['status']=='terminal' for p in parents.values())
    assert all(p['exitCode']==(1 if name=='verify_parent' else 0) for name,p in parents.items())
    runtime = json.loads((PACKAGE/'runtime.json').read_text())
    assert sha(PACKAGE/'runtime.json')==build['runtimePackageSHA256']
    assert all(sha(p)==h for p,h in runtime['filesSHA256'].items())
    meta = BASE/'runtime_metadata';meta.mkdir(exist_ok=False)
    for name in ['runtime.json','buildinfo.json']:
        shutil.copyfile(PACKAGE/name,meta/name)
    licenses=BASE/'upstream_license';licenses.mkdir(exist_ok=False)
    for name in ['COPYING','gpl.txt']:
        shutil.copyfile(Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source')/name,licenses/name)
    compatibility = dict(schema='RAH.EilmerDenseCompatibility.v1',passed=True,
        runtimePackage=str(PACKAGE/'runtime.json'),runtimePackageSHA256=sha(PACKAGE/'runtime.json'),
        binarySHA256=runtime['filesSHA256'][runtime['binary']],
        runReceipt=str(BASE/'small01_run/report.json'),runReceiptSHA256=sha(BASE/'small01_run/report.json'),
        verification=str(BASE/'verify02/report.json'),verificationSHA256=sha(BASE/'verify02/report.json'),
        caseManifestSHA256=run['caseManifestSHA256'],storedFieldAndWallBytesIdentical=True,
        comparedScalars=658,maxScaledError=0,physicalQualification=False,newQualifiedAreaM2=0,
        scope='One four-cell, seven-wall-face, six-update case only; no full-domain memory or convergence qualification.')
    (BASE/'compatibility.json').write_text(json.dumps(compatibility,indent=2)+'\n')

    # 模块二：仅本轮小证据逐字复制；运行包大对象/二进制仅记录路径和SHA。
    entries=[];large=[]
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix=='.pyc':continue
        record=dict(source=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path))
        if record['bytes']>1024**2:
            large.append(record);continue
        target=DEST/path.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        assert sha(target)==record['sha256']
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),**record))
    runtime_files=[dict(source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
                   for p in sorted(PACKAGE.rglob('*')) if p.is_file()]
    sources=['scripts/buildEilmerDenseRuntime.py','scripts/eilmerDenseRuntime.py','scripts/eilmerExplicitProfile.py',
             'scripts/prepareEilmerDevelopmentCase.py','scripts/runPreparedEilmerCase.py','scripts/readEilmerNative.py',
             'scripts/auditEilmerExplicitCompletion.py','scripts/inspectEilmerRun.py',
             'tests/verifyEilmerDenseRuntime.py','patches/eilmer/explicit-dense-checked-v1.patch']
    manifest=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(x['bytes'] for x in entries),
        largeLocalOnly=large,runtimeLocalOnly=runtime_files,finalSourcesSHA256={p:sha(ROOT/p) for p in sources},
        protectedFiles=checks,sharedEngineInputsStable=True,originalFailuresRetained=True,
        uniqueFilesDeleted=0,newQualifiedAreaM2=0)
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

    # 模块三：各阶段父与内部时间分开，不把重叠计时相加为总计算时间。
    budget=dict(sourceAuditSeconds=scope['elapsedSeconds'],failedSourceAuditToolSeconds=read('scope01/failure.json')['toolWallSecondsIncludingStartup'],
        sourceAuditBudgetSeconds=120,buildEntrySeconds=build['elapsedSeconds'],buildBudgetSeconds=300,
        buildParentSeconds=parents['build_parent']['elapsedSeconds'],buildParentLimitSeconds=360,
        compiledDSourceCount=347,productionFlowTargetsRecompiled=1,compilerBuilds=0,
        buildPeakGroupRSSBytes=build['process']['peakObservedGroupRSSBytes'],
        preparationSeconds=prep['elapsedSeconds'],preparationBudgetSeconds=30,preparationParentSeconds=parents['prepare_parent']['elapsedSeconds'],
        nativeSupervisionSeconds=run['commands'][0]['elapsedSeconds'],runEntrySeconds=run['elapsedSeconds'],runEntryBudgetSeconds=30,
        runParentSeconds=parents['run_parent']['elapsedSeconds'],runParentLimitSeconds=60,
        failedVerificationSeconds=first['elapsedSeconds'],verificationSeconds=verify['elapsedSeconds'],
        verificationParentsSeconds={n:p['elapsedSeconds'] for n,p in parents.items() if n.startswith('verify')},
        postprocessingBudgetSeconds=180,postprocessingHasSingleEnclosingParent=False,
        nativeRuns=1,nativeUpdates=6,fullDomainRuns=0,matlabRuns=0,guiRuns=0,
        nativePeakGroupRSSBytes=run['commands'][0]['peakObservedGroupRSSBytes'],
        outputDirectoryLogicalBytes=sum(x['bytes'] for x in entries+large),runtimeDirectoryLogicalBytes=sum(x['bytes'] for x in runtime_files),
        collectionSeconds=time.monotonic()-start,uniqueFilesDeleted=0,physicalQualification=False,newQualifiedAreaM2=0,
        timingScope='Build/preparation/run/each verifier have separate saved parents. Collection/coordination are not inside a single postprocessing parent. Nested times are not additive; editing/reading time excluded. Logical bytes are not newly allocated disk space.')
    (DEST.parent/'budget_receipt.json').write_text(json.dumps(budget,indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),smallBytes=manifest['totalSmallBytes'],runtimeBytes=budget['runtimeDirectoryLogicalBytes'],protected=29,collectionSeconds=budget['collectionSeconds'])))


if __name__=='__main__':collect()
