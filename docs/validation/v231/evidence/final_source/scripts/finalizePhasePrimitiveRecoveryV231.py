"""模块职责：保护旧资产，归档同轨迹捕获/原生拒绝修复及失败，不授整箭或热资格。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from capturePhaseEulerFailureV231 import ROOT,BASE,DOC,PREV,read,save,sha,CLOSED
from preparePhasePrimitiveRecoveryV231 import OLD_BUILD,PREP,BUILD,CLOSED as RECOVERY_CLOSED

OWN=['native/su2PhaseAwareAir/PhaseAirPrimitiveRecovery.hpp','native/su2PhaseAwareAir/primitiveRecoveryControl.cpp',
 'scripts/capturePhaseEulerFailureV231.py','scripts/preparePhasePrimitiveRecoveryV231.py',
 'scripts/finalizePhasePrimitiveRecoveryV231.py','docs/HYBRID_PHASE_PRIMITIVE_RECOVERY_V231.md',
 'docs/validation/v231/PLAN.md','docs/validation/v231/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：旧物理资产及并发候选只读；新链接不替换旧对象或包。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'],r['file']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(OLD_BUILD/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    assets=read(PREV/'final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    for rows in read(ROOT/'outputs/v229_validation/build02/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'],r['path']
    index=read(ROOT/'docs/validation/v230/evidence/local_asset_index.json')['files']
    for n,r in index.items():assert sha(PREV/n)==r['sha256'],n
    records={p.parent.name[8:]:read(p) for p in BASE.glob('process_*/report.json')}
    failures={'capture','object','recovery_prepare','recovery_prepare2','recovery_prepare4'}
    assert len(records)==18,len(records)
    for stage,r in records.items():
        assert r['passed']==(stage not in failures),(stage,r['passed'])
        assert r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None,stage
    assert records['capture']['exitCode']==86
    audit=read(BASE/'audit01/report.json');recovery=read(BASE/'recovery_audit01/report.json')
    assert audit['passed'] and audit['originalFirst30TrajectoryByteIdentical']
    assert audit['termination']['exception_what']=='RAH_PHASE_ENERGY_ENTROPY_DOMAIN'
    assert recovery['passed'] and recovery['rejectedTrialsRecoveredByNativeOldState']==5
    assert recovery['invalidOldStillRejected'] and not recovery['scientificAcceptance']
    assert sha(BASE/'capture01/termination_raw_state.csv')==audit['rawStateSHA256']
    for n,r in audit['comparison'].items():assert sha(BASE/'capture01'/n)==sha(PREV/'rocket01'/n)==r['sha256']
    assert (BASE/'control_original/primitive.csv').read_bytes().splitlines()[:3]==(BASE/'control_recovered/primitive.csv').read_bytes().splitlines()[:3]
    for name in ('object','control'):
        for n,h in read(PREP/(name+'_dependencies.json')).items():assert sha(n)==h,n
    sources=read(PREP/'sources.json')
    for n,h in sources.items():assert sha(n)==h,n
    for n in ('PhaseAirPrimitiveRecovery.hpp','primitiveRecoveryControl.cpp'):
        assert sha(PREP/n)==sha(ROOT/'native/su2PhaseAwareAir'/n)
    # 编译后的依赖准确对应预登记文件；最终驱动只在原命令上加新对象/输出路径。
    actual_dependencies={}
    expected={**read(PREV/'build02/dependencies.json'),**read(PREP/'object_dependencies.json'),**read(PREP/'control_dependencies.json')}
    for dep in [BUILD/'primitive.d',BUILD/'original.d',BUILD/'recovered.d',BASE/'recovery_build03/rocket.d']:
        text=dep.read_text(encoding='utf8').replace('\\\n',' ');tokens=text[text.index(': ')+2:].split()
        for n in tokens:
            p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
            assert expected[str(p)]==h,str(p);actual_dependencies[str(p)]=h
    cmd=read(BASE/'recovery_build03/command.json');original=read(PREV/'build02/command.json')
    strip=lambda c:[x for i,x in enumerate(c) if x not in ('-o','-MF') and (i==0 or c[i-1] not in ('-o','-MF'))]
    normalized=[x for x in strip(cmd) if x not in (str(BUILD/'primitive.obj'),'-I'+str(ROOT/'native/su2PhaseAwareAir'))]
    assert normalized==strip(original)
    assert records['driver2']['inputSHA256'][str(BUILD/'primitive.obj')]==sha(BUILD/'primitive.obj')
    assert CLOSED and RECOVERY_CLOSED
    for name in ('capturePhaseEulerFailureV231.py','preparePhasePrimitiveRecoveryV231.py'):
        r=subprocess.run([sys.executable,str(ROOT/'scripts'/name),'capture'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        assert r.returncode!=0 and b'closed' in r.stdout
    # 模块二：仅紧凑观测/源码/配置公开，完整未接受原始状态及二进制留E盘。
    out.mkdir();evidence.mkdir()
    for folder in BASE.glob('process_*'):shutil.copytree(folder,evidence/folder.name)
    for stage in ('prepare01','audit01','recovery_prepare01','recovery_prepare02','recovery_prepare03',
                  'recovery_prepare04','recovery_audit01','control_original','control_recovered'):
        shutil.copytree(BASE/stage,evidence/stage)
    target=evidence/'capture01';target.mkdir()
    for n in ('case.cfg','contract.json','termination.txt','history.csv','physical_history.csv'):
        shutil.copyfile(BASE/'capture01'/n,target/n)
    for stage in ('recovery_build01','recovery_build02','recovery_build03'):
        target=evidence/stage;target.mkdir()
        for p in (BASE/stage).iterdir():
            if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    shutil.copyfile(BASE/'raw_scan01.json',evidence/'unmetered_raw_scan01.json')
    save(evidence/'actual_dependencies.json',actual_dependencies)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in
      [BASE/'capture01/termination_raw_state.csv',BASE/'capture01/state_000000.csv',BASE/'capture01/state_000001.csv',
       BUILD/'primitive.obj',BUILD/'original.exe',BUILD/'recovered.exe',BUILD/'rocket.exe',BASE/'recovery_build03/rocket.exe']}
    save(evidence/'local_asset_index.json',dict(files=local,base='outputs/v231_validation',binariesPublished=False,
      rawFailureStatePublished=False,preferredUnexecutedDriver='recovery_build03/rocket.exe'))
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    totals={kind:sum(r['elapsedSeconds'] for r in records.values() if r['budgetKind']==kind) for kind in ('offline','build','query','cfd')}
    elapsed=time.monotonic()-started
    assert totals['offline']+elapsed<120 and totals['build']<120 and totals['query']<30 and totals['cfd']<60
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<256*1024**2
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,scientificAcceptance=False,
      originalFailureCaptured=True,originalFirst30TrajectoryByteIdentical=True,actualWhat='RAH_PHASE_ENERGY_ENTROPY_DOMAIN',
      nativePrimitiveDomainRejectionFixed=True,actualFullExceptionStackKnown=False,uniqueUnderlyingCauseClaimed=False,
      actualRawEnergyInvalidPoint=17337,comparisonCases=8,recoveredTrialClasses=5,oldInvalidRejected=True,
      maxIndependentPrimitiveRelativeError=recovery['maxIndependentPrimitiveRelativeError'],
      recoveryFullRocketDriverCompiled=True,recoveryFullRocketDriverExecuted=False,
      preferredDriver='outputs/v231_validation/recovery_build03/rocket.exe',preferredDriverSHA256=sha(BASE/'recovery_build03/rocket.exe'),
      primitiveObjectSHA256=sha(BUILD/'primitive.obj'),actualDependencies=len(actual_dependencies),newHeatAreaM2=0,newQualifiedAreaM2=0,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),engineeringCandidateCount=len(candidates['files']),
      oldObjectsUnchanged=len(old),oldDependenciesUnchanged=len(deps),assets=assets,processSeconds=totals,
      finalizationSeconds=elapsed,outputBytesBeforeReport=size,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      allGroupsEmpty=True,unmeteredTargetedReadOnlyInspectionAndRawScan=True,allFailuresPreserved=True,
      oldPhysicsDefaultUnchanged=True,phaseDomainUnchanged=True,MATLABExecuted=False,sharedWebModified=False,
      browserInteractionTested=False,newMeshGenerated=False,uniqueAssetsDeleted=False,executionCLIsClosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v231/evidence_manifest.json']+['docs/validation/v231/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,processSeconds=totals,
      fullRocketRecoveryExecuted=False,scientificAcceptance=False,preferredDriverSHA256=report['preferredDriverSHA256'])))

if __name__=='__main__':main()
