"""模块职责：封存本轮真实发展/失败、独立审核与科研出口，保护并发工作并生成精确发布清单。"""
from pathlib import Path
import ast
import json
import shutil
import subprocess
import sys
import time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v251_validation';DOC=ROOT/'docs/validation/v251'
CLIS=['scripts/'+n+'.py' for n in ('runNativeILUV251','auditNativeILUV251','exportNativeILUV251','recoverILUSurfaceV251')]
OWN=CLIS+['scripts/finalizeNativeILUV251.py','native/su2PhaseAwareAir/nativeILUV251.cpp','tests/testNativeILUV251.js',
 'docs/NATIVE_ILU_ACCEPTANCE_V251.md','docs/validation/v251/PLAN.md','docs/validation/v251/SURFACE_RECOVERY.md','docs/validation/v251/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(name):
    return BASE/'owned_source'/name if name in SPLIT else ROOT/name


def main():
    begin=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：原资产及并发本地补充只核对，不将无关修改加入本轮发布。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for item in engineering['files']:assert sha(ROOT/item['path'])==item['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for item in candidates['files']:
        assert sha(ROOT/item['file'])==item['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/item['file'])==item['candidateSHA256']
    prev=read(ROOT/'docs/validation/v250/evidence/final_report.json');assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for name,h in dict(assets,**unrelated).items():assert sha(ROOT/name)==h,name
    private=read(BASE/'concurrent_prefixes.json');split={}
    for name in SPLIT:
        live=(ROOT/name).read_text(encoding='utf8');prefix=private[name]['prefix'];own=source(name).read_text(encoding='utf8')
        assert live.count(prefix)==1 and live.replace(prefix,'',1)==own
        assert 'public_benchmarks_20261116' not in own
        split[name]=dict(localSupplementPreserved=True,publishedOwnVersionSHA256=sha(source(name)),combinedWorktreeSHA256=sha(ROOT/name))
    # 模块二：所有执行回执与准确运行源码绑定；科学拒绝、资源失败和正常退出分别保存。
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    assert {'process_prepare','process_prepare2','process_build','process_rocket','process_audit2','process_recovery','process_export','process_export2'}.issubset(records)
    totals={};failures={};bindings=[];required=set();cache={};capacities={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        if not r['passed']:failures[name]=r['exitCode']
        capacity=read(BASE/name/'terminal_capacity.json');capacities[name]=capacity
        for path,h in r['inputSHA256'].items():
            if path not in cache:cache[path]=sha(path)
            if cache[path]==h:continue
            matches=list(BASE.glob('process_*/execution_sources/'+h+Path(path).suffix))
            matches += [p for p in BASE.glob('failed_source*/*') if p.is_file() and sha(p)==h]
            matches += [p for p in BASE.glob('stage*/'+Path(path).name) if p.is_file() and sha(p)==h]
            assert matches,(name,path)
            snapshot=matches[0];assert sha(snapshot)==h;required.add(snapshot)
            bindings.append(dict(process=name,path=path,executedSHA256=h,currentSHA256=cache[path],snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<=120 and totals['build']<=120 and totals['cfd']<=1801 and totals['audit']<=180
    audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench02/report.json');recovery=read(BASE/'recovery01/report.json')
    assert audit['passed'] and ui['passed'] and not audit['physicalQualification']
    assert ui['updates']==audit['completeSavedUpdates']+[recovery['update']] and ui['faces']==19752 and audit['initialBitEqual']
    assert recovery['passed'] and not recovery['nativeFullFrame'] and ui['nativeCompleteSavedUpdates']==audit['completeSavedUpdates']
    assert audit['nativeNumericalStatementsUnchanged'] and read(BASE/'build01/report.json')['matureObjectsRecompiled']==0
    assert all(c['passed'] for c in capacities.values())
    assert read(BASE/'workbench02/frontend_test.txt')['passed']
    closed={}
    for name in OWN:
        if name.endswith('.py'):ast.parse(source(name).read_text(encoding='utf8'))
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        msg=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in msg;closed[name]=msg.strip()
    # 模块三：只发布小证据，完整流场和成熟二进制留本地；原执行源与最终关闭版分别封存。
    out.mkdir();evidence.mkdir();local={}
    def copy(path,rel=None):
        dest=evidence/(rel if rel is not None else path.relative_to(BASE));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or out in path.parents or '__pycache__' in path.parts:continue
        rel=path.relative_to(BASE);top=rel.parts[0]
        keep=path in required or path.suffix in ('.json','.txt','.csv','.png','.patch') or top.startswith(('process_','failed_source'))
        if path.suffix in ('.bin','.exe','.d','.a','.obj','.su2','.npz','.pyc') or top=='owned_source':keep=False
        if top=='publication_attempt01':keep=False
        if path.name in ('concurrent_prefixes.json','rocket_newton_development.json') or path.name.startswith('original_faces_'):keep=False
        if path.stat().st_size>=2*1024**2:keep=False
        if path in required:assert keep,path
        if keep:copy(path)
        else:local[rel.as_posix()]=dict(bytes=path.stat().st_size,sha256=sha(path))
    for name in OWN:copy(source(name),Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,originalFaces=19752,originalAreaM2=1.9646292565715553,
      finePoints=335726,volumeCells=705490,completedUpdates=audit['completedUpdates'],completeSavedUpdates=audit['completeSavedUpdates'],
      nativeExitedNormally=audit['nativeExitedNormally'],stopReason=audit['stopReason'],observedMatrixFreeProducts=audit['observedMatrixFreeProducts'],
      nativeLinearHistory=audit['nativeLinearHistory'],linearRefusal=audit['linearRefusal'],nativeMemoryGraph=audit['nativeMemoryGraph'],jointPressureTemperatureRejections=audit['jointPressureTemperatureRejections'],
      reconstructedOldSurface=recovery,availableResearchSurfaceUpdates=ui['updates'],
      recoveryCount=audit['recoveryCount'],recoveryRoles=audit['recoveryRoles'],terminal=audit['terminal'],benefitDecisionReached=audit['benefitDecisionReached'],
      firstUpdateBitEqual=audit['firstUpdateBitEqual'],firstUpdateMaximumScaledDifference=audit['firstUpdateMaximumScaledDifference'],
      independentSavedStates=audit['independentSavedStates'],nativeNumericalStatementsUnchanged=True,
      physicalQualification=False,newQualifiedAreaM2=0,nonlinearConvergenceProved=False,fullJacobianAccuracyProved=False,
      scientificImagesVisuallyInspected=True,sourceAndNumericTestsPassed=True,exportedFields=ui['fields'],originalFaceValuesVerified=ui['exactResearchValuesChecked'],
      readOnlyObserverObjectsCompiled=0,matureObjectsRecompiled=0,sharedWebGLChanged=False,MATLABExecuted=False,browserInteractionTested=False,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
      allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,sourceBindingPreparationFailurePreserved=True,
      auditAdmissionPreProcessFailurePreserved=True,firstFrontendSourceReconstructedWithExactAdmissionSHA=True,
      publicationWhitespaceCorrection=read(BASE/'publication_review.json'),
      terminalCapacities=capacities,processSeconds=totals,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),historicalSourceBindings=bindings,
      closedCLIs=closed,ownedFiles=OWN,uniqueDataDeleted=False,newPackageInstalled=False,unmeteredDevelopmentReadsAndPublicationDisclosed=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<1536*1024**2
    report['finalizationSeconds']=time.monotonic()-begin;save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    files=OWN+['docs/validation/v251/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    assert sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())<1536*1024**2
    print(json.dumps(dict(passed=True,files=len(files),evidenceFiles=len(manifest),processSeconds=totals)))


if __name__=='__main__':main()
