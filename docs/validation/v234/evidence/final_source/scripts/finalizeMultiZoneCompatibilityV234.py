"""模块职责：保存整体失败/单区审核/准确源和资源账，核对保护文件，仅归档小证据。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v234_validation';DOC=ROOT/'docs/validation/v234'
OWN=['native/su2PhaseAwareAir/multiZoneCompatibilityV234.cpp',
 'scripts/runMultiZoneCompatibilityV234.py','scripts/observeMultiZoneFailureV234.py',
 'scripts/recoverMultiZoneControlV234.py','scripts/stagePhaseMultiZoneV234.py',
 'scripts/auditMultiZoneCompatibilityV234.py','scripts/finalizeMultiZoneCompatibilityV234.py',
 'docs/HYBRID_MULTIZONE_COMPATIBILITY_V234.md','docs/validation/v234/PLAN.md',
 'docs/validation/v234/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    # 模块一：原始终态和已修改执行源绑定；原single缺回执不能补造。
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==18
    failures={'process_prepare':1,'process_capture':86,'process_multi2':1,'process_multi3':1,
      'process_build_mz':15,'process_multi4':3221225477}
    historical=[]
    old_cpp=[BASE/'observer01/originalHelper.cpp',BASE/'observer01/controlObserved.cpp',BASE/'recovery_prepare01/control.cpp']
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert r['exitCode']==failures.get(name,0) and r['passed']==(name not in failures),name
        assert r['stopReason']==('output' if name=='process_build_mz' else None),name
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            options=[BASE/name/'execution_sources'/(h+Path(p).suffix)]+old_cpp
            match=next((q for q in options if q.is_file() and sha(q)==h),None)
            assert match is not None,('EXECUTED_SOURCE_NOT_PRESERVED',name,p,h)
            historical.append(dict(process=name,path=p,executedSHA256=h,snapshotPath=str(match),currentSHA256=sha(p)))
    assert not (BASE/'process_single/report.json').exists()
    a=read(BASE/'audit01/report.json');assert a['passed'] and a['singleControlVerified'] and not a['dualCompatibilityPassed']
    assert a['productionDualRouteClosed'] and a['newQualifiedAreaM2']==0
    for r in read(BASE/'audit01/source_binding.json'):assert sha(r['path'])==r['sha256']
    for group in read(BASE/'mz_build01/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256']
    # 完整母域/气体资产SHA保持，单区旧场不再次计算或覆盖。
    assets=read(ROOT/'docs/validation/v222/evidence/final_report.json')['assets']
    for p,h in assets.items():assert sha(ROOT/p)==h,p
    closed={}
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        if n.startswith('scripts/') and n!=OWN[6] and 'finalize' not in n:
            proc=subprocess.run([sys.executable,str(p),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
            msg=proc.stdout.decode('utf8',errors='replace');assert proc.returncode!=0 and 'closed' in msg,(n,msg)
            closed[n]=msg.strip()
    assert len(closed)==5
    # 模块二：小控制/失败数据完整公开；大归档和EXE只索引，不复制运行库。
    out.mkdir();evidence.mkdir()
    for p in BASE.glob('process_*'):shutil.copytree(p,evidence/p.name)
    for n in ('prepare01','observer01','recovery_prepare01','mz_prepare01','single01','single02',
              'multi01','multi02','multi03','multi04','capture01','capture04','audit01'):
        shutil.copytree(BASE/n,evidence/n)
    for n in ('build01','build02','build03','mz_build01','mz_link02'):
        target=evidence/n;target.mkdir()
        for p in (BASE/n).iterdir():
            if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    # 记录原生崩溃签名，不公开Windows报告中的无关机器标识。
    wer_paths=[Path('C:/ProgramData/Microsoft/Windows/WER/ReportArchive')/n/'Report.wer' for n in (
      'AppCrash_control.exe_39ce3f48425e95172535d7469df7f3185d9c83_5560adf6_84a152c2-55d1-45ef-bd30-c892ffa1708f',
      'AppCrash_control.exe_fb5dbba92b8fe6a61e60551bf580ae7ecfc05b_5560adf6_70b8ab8d-e712-4acf-a4b3-c662dc5f7273')]
    wer=[]
    for p in wer_paths:
        if not p.is_file():wer.append(dict(path=str(p),present=False));continue
        txt=p.read_text(encoding='utf16');selected=[s for s in txt.splitlines() if s.startswith(('EventType=','Sig[','AppName=','AppPath='))]
        wer.append(dict(path=str(p),present=True,sha256=sha(p),signature=selected))
    save(evidence/'failure_signatures.json',wer)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p))
      for n in ('build01','build02','build03','mz_build01','mz_link02') for p in (BASE/n).iterdir() if p.suffix in ('.exe','.obj','.a')}
    save(evidence/'local_asset_index.json',dict(base='outputs/v234_validation',files=local,
      largeAssetsBackedUpToGitHub=False,closedExperimentalIdentity='RAH_PHASE_AIR_MZ=17',
      existingSourcesAndArchivesRequired='v220 build02 / v167 exact source tree; no automatic rerun'))
    for n in OWN:
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,target)
    totals={}
    for r in records.values():totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
    elapsed=time.monotonic()-start;assert totals['offline']+elapsed<120
    assert totals['build']<120 and totals['identity-build']<180 and totals['recovery-query']<60
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<256*1024**2
    report=dict(passed=True,goalStatus='active',scientificAcceptance=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
      flowUpdates=0,singleControlVerified=True,dualCompatibilityPassed=False,dualProductionRouteClosed=True,
      missingViscousSecondaryConfirmed=True,missingNewIdentityThermalDispatchConfirmed=True,
      GDBProcessExitZeroIsNotNativeSuccess=True,originalSingleReceiptMissing=True,originalUnknownQueryReservedSeconds=60,
      completedReceipts=len(records),allRecordedGroupsEmpty=True,initialMissingReceiptGroupNotRetroactivelyClaimed=True,
      historicalSourceBindings=historical,protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),
      candidateCount=len(candidates['files']),assets=assets,closedCLIs=closed,processSeconds=totals,finalizationSeconds=elapsed,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      outputBytesIncludingExternalJunction=size,externalJunctionExistingBytes=53848128,
      originalWallFaces=19752,originalAreaM2=1.9646292565715553,
      MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,uniqueAssetsDeleted=False,
      newRocketMeshGenerated=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v234/evidence_manifest.json']+['docs/validation/v234/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,outputBytes=size,closedCLIs=len(closed))))


if __name__=='__main__':main()
