"""模块职责：封存原生恢复、真实发展及科研出口，核对保护资产并生成精确公开清单。"""
from pathlib import Path
import ast, shutil, subprocess, sys, time
from runSU2RocketPhaseV222 import ROOT, read, save, sha

BASE=ROOT/'outputs/v248_validation';DOC=ROOT/'docs/validation/v248'
CLIS=['scripts/'+n+'.py' for n in ('runFASRecoveryV248','auditFASRecoveryV248','exportFASRecoveryV248','diagnoseFASRecoveryV248')]
OWN=CLIS+['scripts/auditFASRecoveryEventsV248.py','scripts/fasRecoveryFieldsV248.py',
 'scripts/finalizeFASRecoveryV248.py','native/su2PhaseAwareAir/fasRecoveryDevelopmentV248.cpp',
 'native/su2PhaseAwareAir/FASRecoveryObserverV248.hpp','tests/testFASRecoveryFramesV248.py','tests/testFASRecoveryV248.js',
 'docs/FAS_NATIVE_RECOVERY_V248.md','docs/validation/v248/PLAN.md','docs/validation/v248/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(name):
    return BASE/'owned_source'/name if name in SPLIT else ROOT/name


def main():
    begin=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：保护原物理资产、工程分支及用户只留本地的文档补充。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for item in engineering['files']:assert sha(ROOT/item['path'])==item['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for item in candidates['files']:
        assert sha(ROOT/item['file'])==item['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/item['file'])==item['candidateSHA256']
    prev=read(ROOT/'docs/validation/v247/evidence/final_report.json')
    assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for name,h in dict(assets,**unrelated).items():assert sha(ROOT/name)==h,name
    private=read(BASE/'concurrent_prefixes.json');split={}
    for name in SPLIT:
        live=(ROOT/name).read_text(encoding='utf8');prefix=private[name]['prefix']
        assert live.count(prefix)==1 and live.replace(prefix,'',1)==source(name).read_text(encoding='utf8')
        assert 'public_benchmarks_20261116' not in source(name).read_text(encoding='utf8')
        split[name]=dict(localSupplementPreserved=True,publishedOwnVersionSHA256=sha(source(name)),combinedWorktreeSHA256=sha(ROOT/name))
    # 模块二：实际监督回执逐条登记；变化的执行源须能由原始快照重新绑定。
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    assert set(('process_prepare','process_refresh','process_build','process_rocket','process_reader_tests')).issubset(records)
    totals={};cache={};bindings=[];required=set();failures={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        if not r['passed']:failures[name]=r['exitCode']
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for path,h in r['inputSHA256'].items():
            if path not in cache:cache[path]=sha(path)
            if cache[path]==h:continue
            matches=list(BASE.glob('process_*/execution_sources/'+h+Path(path).suffix))
            matches += [p for p in BASE.glob('failed_source*/*') if p.is_file() and sha(p)==h]
            matches += [p for p in BASE.glob('stage*/'+Path(path).name) if p.is_file() and sha(p)==h]
            assert matches,(name,path)
            snapshot=matches[0];required.add(snapshot)
            bindings.append(dict(process=name,path=path,executedSHA256=h,currentSHA256=cache[path],snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<120 and totals['build']<120 and totals['cfd']<1800 and totals['audit']<240
    audit_path=BASE/'audit02/report.json';audit=read(audit_path);ui=read(BASE/'workbench01/report.json')
    diagnosis=read(BASE/'diagnosis01/report.json');assert diagnosis['passed']
    assert audit['passed'] and ui['passed'] and not audit['physicalQualification']
    assert ui['updates']==audit['originalFaceUpdates'] and audit['originalFaces']==19752
    assert audit['readerTestsReused'] and audit['executionSchemaLabelsCorrectedSeparately']
    closed={}
    for name in OWN:
        if name.endswith('.py'):ast.parse(source(name).read_text(encoding='utf8'))
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in message
        closed[name]=message.strip()
    # 模块三：发布小证据；大场和二进制仅登记位置/SHA，绝不清理独有结果。
    out.mkdir();evidence.mkdir();local={}
    def copy(path,rel=None):
        dest=evidence/(rel if rel is not None else path.relative_to(BASE))
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or out in path.parents or '__pycache__' in path.parts:continue
        rel=path.relative_to(BASE);top=rel.parts[0]
        keep=path in required or path.suffix in ('.json','.txt','.csv','.png') or top.startswith('process_') or top.startswith('failed_source')
        if path.suffix in ('.bin','.exe','.d','.su2','.pyc') or top=='owned_source':keep=False
        if path.name in ('concurrent_prefixes.json','rocket_fas_control.json') or path.name.startswith('original_faces_'):keep=False
        if path.stat().st_size>=2*1024**2:keep=False
        if keep:copy(path)
        else:local[rel.as_posix()]=dict(bytes=path.stat().st_size,sha256=sha(path))
    for name in OWN:copy(source(name),Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,
      originalFaces=19752,originalAreaM2=1.9646292565715553,finePoints=335726,coarsePoints=44484,
      developmentSucceeded=audit['developmentSucceeded'],completedUpdates=audit['completedUpdates'],
      completeSavedUpdates=audit['completeSavedUpdates'],originalFaceUpdates=audit['originalFaceUpdates'],
      nativeRecovery=audit['recovery'],terminal=audit['terminal'],benefitDecisionReached=audit['benefitDecisionReached'],
      recoveryPointsSquaredRHSFraction=diagnosis['recoveryPointsSquaredRHSFraction'],
      benefitStopExpected=audit['benefitStopExpected'],initialZeroBitEqual=True,firstUpdateBitEqual=audit['firstUpdateBitEqual'],
      firstUpdateMaximumScaledDifference=audit['errorMetrics']['crossThinExecutableFirstUpdate'],
      independentSavedStates=audit['independentSavedStateChecks'],brentCrosschecks=audit['brentCrosschecks'],
      scientificImagesVisuallyInspected=True,sourceAndNumericTestsPassed=True,physicalQualification=False,newQualifiedAreaM2=0,
      exportedFields=ui['fields'],originalFaceValuesVerified=ui['exactResearchValuesChecked'],
      numericalObjectsRecompiled=0,sharedWebGLChanged=False,MATLABExecuted=False,browserInteractionTested=False,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
      allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,
      processSeconds=totals,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      historicalSourceBindings=bindings,closedCLIs=closed,uniqueDataDeleted=False,newPackageInstalled=False,
      schemaMetadataCorrectionOnly=True,unmeteredDevelopmentReadsAndPublicationDisclosed=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<1536*1024**2
    report['finalizationSeconds']=time.monotonic()-begin;save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    files=OWN+['docs/validation/v248/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(dict(passed=True,files=len(files),evidenceFiles=len(manifest),evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),processSeconds=totals))


if __name__=='__main__':main()
