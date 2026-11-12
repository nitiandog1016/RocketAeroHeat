"""模块职责：保留失败/原执行，公开紧凑恢复数据，明确未重跑修复及并发保护。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from runPhaseNoseHeatingV227 import SOURCES
BASE=ROOT/'outputs/v227_validation';DOC=ROOT/'docs/validation/v227'
OWN=['solvePhaseAwareLaminarProfile.m','solvePhaseAwareNoseHeating.m','runPhaseAwareNoseHeatingCase.m',
 'buildPhaseAwareNoseHeatingViewerData.m','runAeroHeatSingleCase.m','buildAerothermalViewerData.m',
 'tests/testPhaseNoseHeatingV227.m','tests/testPhaseNoseHeatingViewerV227.js',
 'scripts/runPhaseNoseHeatingV227.py','scripts/auditPhaseNoseHeatingV227.py','scripts/finalizePhaseNoseHeatingV227.py',
 'docs/HYBRID_PHASE_NOSE_HEATING_V227.md','docs/validation/v227/PLAN.md','docs/validation/v227/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
def compact(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n',encoding='utf8')
    assert read(p)==obj
def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    other=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in other['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    assets=read(ROOT/'outputs/v226_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    changed=[]
    for n in SOURCES:
        if sha(ROOT/n)!=sha(BASE/'matlab02/executed_source'/n):changed.append(n)
    assert changed==['solvePhaseAwareLaminarProfile.m','buildPhaseAwareNoseHeatingViewerData.m'],changed
    recovery=read(BASE/'audit02/wall_recovery.json')
    assert recovery['correctedCodeSHA256']==sha(ROOT/'solvePhaseAwareLaminarProfile.m')
    old_view=(BASE/'matlab02/executed_source/buildPhaseAwareNoseHeatingViewerData.m').read_text(encoding='utf8')
    expected=old_view.replace("data.missingReason=h.reason;data.region(h.validFace)={'相态球鼻热工程候选'};",
      "data.missingReason=h.reason;data.region=repmat({'其他待求解区域'},size(mesh.F,1),1);\n"
      "data.region(h.validFace)={'相态球鼻热工程候选'};\n"
      "data.regions=struct('name',{'相态球鼻热工程候选','其他待求解区域'}, ...\n"
      "    'faces',{nnz(h.validFace),nnz(~h.validFace)}, ...\n"
      "    'area',{sum(data.faceArea(h.validFace)),sum(data.faceArea(~h.validFace))});")
    assert expected==(ROOT/'buildPhaseAwareNoseHeatingViewerData.m').read_text(encoding='utf8')
    records=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
    for name,wanted in [('matlab01',False),('matlab02',True),('audit01',False),('audit02',True)]:
        assert read(BASE/('process_'+name)/'report.json')['passed']==wanted
    assert read(BASE/'audit02/report.json')['passed']
    retained=[]
    for folder in ('matlab01','matlab02'):
        target=evidence/folder;target.mkdir()
        for p in (BASE/folder).iterdir():
            if p.name in ('heating.json','rocket_workbench.json'):
                retained.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size,retainedLocally=True,authoritative=False));continue
            if p.is_dir():shutil.copytree(p,target/p.name)
            else:shutil.copyfile(p,target/p.name)
    for src in BASE.glob('process_*'):
        if src.is_dir():shutil.copytree(src,evidence/src.name)
    for name in ('process_matlab01','process_matlab02'):
        raw=BASE/name/'stdout.txt';decoded=raw.read_bytes().decode('gb18030',errors='replace')
        (evidence/name/'stdout_utf8.txt').write_text(decoded,encoding='utf8')
        save(evidence/name/'decoding.json',dict(sourceSHA256=sha(raw),encodingAssumption='gb18030',replacementCharacters=decoded.count('\ufffd')))
    # 已审核数值只改JSON空白/Unicode表达；区域列表与当前导出器静态同规则恢复。
    recovered=evidence/'recovered';recovered.mkdir();h=read(BASE/'audit02/heating.json');m=read(BASE/'audit02/rocket_workbench.json')
    names=['相态球鼻热工程候选','其他待求解区域'];v=h['validFace'];areas=m['faceArea']
    m['region']=[names[0] if flag else names[1] for flag in v]
    m['regions']=[dict(name=name,faces=sum(r==name for r in m['region']),area=sum(a for a,r in zip(areas,m['region']) if r==name)) for name in names]
    assert sum(r['faces'] for r in m['regions'])==19752 and abs(sum(r['area'] for r in m['regions'])-m['totalArea'])<1e-12
    assert m['regions'][0]['faces']==1306 and abs(m['regions'][0]['area']-h['summary']['heatCandidateAreaM2'])<1e-15
    m['researchIdentity']['regionListRestoredOffline']=True
    compact(recovered/'heating.json',h);compact(recovered/'rocket_workbench.json',m)
    for name in ('report.json','wall_recovery.json','frontend.txt','rocket_nose_heat.png'):shutil.copyfile(BASE/'audit02'/name,recovered/name)
    for name in ('heating.json','rocket_workbench.json'):
        p=BASE/'audit02'/name;retained.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size,retainedLocally=True,
          note='compact equivalent published; workbench region list separately reconciled'))
    save(evidence/'local_raw_assets.json',retained)
    # 与审核当时的共享前端准确绑定；不把共享未提交代码复制进本交付。
    assert sha(ROOT/'web/viewer.js')==read(BASE/'audit02/report.json')['frontendSHA256']
    result=subprocess.run(['node',str(ROOT/'tests/testPhaseNoseHeatingViewerV227.js'),str(recovered)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout
    shutil.copyfile(out/'frontend.txt',evidence/'final_frontend.txt')
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        q=evidence/'final_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseNoseHeatingV227.py'),'matlab01'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    total={kind:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==kind) for kind in ('matlab','audit')}
    assert total['matlab']<360 and total['audit']<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    elapsed=time.monotonic()-start;assert elapsed<120
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,physicalQualification=False,newQualifiedAreaM2=0,
      heatCandidate=h['summary'],protectedCount=29,sharedEngineeringArtifactsUnchanged=len(other['files']),
      sharedCandidatesUnchanged=True,oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,
      originalGeometryUnchanged=True,rawMatlabScientificAcceptance=False,currentMatlabCorrectionExecuted=False,
      postExecutionChanges=changed,offlineWallRecoveryPassed=True,regionListStaticRecoveryPassed=True,
      firstFailuresRetained=True,processBudgetSeconds=total,finalizationSeconds=elapsed,outputBytesBeforeFinalReport=size,
      validationEntryClosed=True,CFDRuns=0,nativeQueries=0,engineBuilds=0,browserInteractionTested=False,
      sharedWebSourceCopied=False,finalScientificImageVisuallyInspected=True,uniqueAssetsDeleted=False,
      completeReplacementMatlabFiles=[n for n in OWN if n.endswith('.m')])
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v227/evidence_manifest.json']+['docs/validation/v227/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,files=len(manifest),seconds=elapsed,outputBytes=size,processBudgetSeconds=total,
      rawMatlabScientificAcceptance=False,currentMatlabCorrectionExecuted=False)))
if __name__=='__main__':main()
