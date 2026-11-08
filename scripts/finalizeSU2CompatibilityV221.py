"""模块职责：只读核对历史资产、构建对象、准确执行源和本轮小证据。"""
from pathlib import Path
import ast,hashlib,json,shutil,subprocess,sys,time
from runSU2CompatibilityV221 import ROOT,BASE,BUILD,read,save
EV=ROOT/'docs/validation/v221/evidence';DOC=EV.parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);EV.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    assets=read(ROOT/'outputs/v220_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    for rows in read(ROOT/'outputs/v220_validation/build02/members.json').values():
        for row in rows:assert sha(Path(row['path']))==row['sha256']
    inputs=read(BUILD/'inputs02.json')
    for n,h in inputs.items():assert sha(Path(n))==h,n
    assert read(BASE/'restart_audit01/report.json')['passed']
    scripts=['inspectSU2CompatibilityV221.py','runSU2CompatibilityV221.py','su2PhaseAirContract.py',
             'verifySU2CompatibilityV221.py','verifySU2RestartV221.py','finalizeSU2CompatibilityV221.py']
    for name in scripts:ast.parse((ROOT/'scripts'/name).read_text(encoding='utf8'))
    result=subprocess.run([sys.executable,str(ROOT/'scripts/runSU2CompatibilityV221.py'),'control'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert result.returncode!=0 and b'closed' in result.stdout
    processes=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];processes.append(r)
        shutil.copytree(p.parent,EV/p.parent.name)
    for name in ('control01','restart01','audit01','restart_audit01','source01','build01_source'):
        shutil.copytree(BASE/name,EV/name)
    for p in BASE.glob('entry_*.py'):shutil.copyfile(p,EV/p.name)
    (EV/'build').mkdir()
    for p in BUILD.iterdir():
        if p.suffix in ('.json','.d'):shutil.copyfile(p,EV/'build'/p.name)
    # 全部小状态都归档；引擎/对象及大箭体网格继续留E盘。
    for name in ('COPYING','AUTHORS.md'):shutil.copyfile(ROOT/'docs/validation/v220/evidence'/name,EV/name)
    candidates={}
    def add(data):candidates[hashlib.sha256(data).hexdigest()]=data
    def variants(text):
        text=text.replace('\r\n','\n');add(text.encode('utf8'));add(text.replace('\n','\r\n').encode('utf8'))
    for p in EV.rglob('*'):
        if p.is_file() and p.suffix in ('.py','.cpp','.md'):add(p.read_bytes())
    own=[ROOT/'scripts'/n for n in scripts]+[ROOT/'native/su2PhaseAwareAir/compatibility.cpp',DOC/'PLAN.md']
    for p in own:variants(p.read_text(encoding='utf8'))
    lines=(DOC/'PLAN.md').read_text(encoding='utf8').splitlines(True)
    for i in range(1,len(lines)+1):variants(''.join(lines[:i]).rstrip()+'\n')
    # 早期物性合同没有全局坐标检查；仅删除新增块恢复当时SHA，绝不追改其回执。
    identity=(ROOT/'scripts/su2PhaseAirContract.py').read_text(encoding='utf8')
    begin=identity.index('    # 原生会重排局部节点');end=identity.index('    packet=convert_flow',begin)
    variants(identity[:begin]+identity[end:])
    wanted={}
    for r in processes:
        for n,h in r['inputSHA256'].items():
            p=Path(n)
            if p in own or p==ROOT/'scripts/runSU2PhaseAirV220.py':wanted.setdefault(h,set()).add(n)
    snapshots=EV/'execution_sources';snapshots.mkdir();manifest=[];missing=[]
    for h,names in wanted.items():
        if h not in candidates:missing.append(dict(sha256=h,paths=sorted(names)));continue
        name=h+Path(next(iter(names))).suffix;(snapshots/name).write_bytes(candidates[h])
        manifest.append(dict(sha256=h,paths=sorted(names),snapshot=name))
    save(EV/'execution_source_manifest.json',manifest);save(EV/'missing_source_snapshots.json',missing);assert not missing,missing
    files=own+[ROOT/'AGENTS.md',ROOT/'README.md',ROOT/'docs/PROJECT_GOALS.md',ROOT/'docs/HYBRID_SU2_COMPATIBILITY_V221.md',DOC/'NEXT_DECISION.md',ROOT/'.gitattributes']
    for p in files:
        rel=p.relative_to(ROOT) if p.name!='.gitattributes' else Path('.gitattributes.snapshot')
        target=EV/'final_source'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    totals={k:sum(r['elapsedSeconds'] for r in processes if r['budgetKind']==k) for k in ('offline','build','query')}
    assert totals['offline']<=120 and totals['build']<=240 and totals['query']<=60
    report=dict(passed=True,protectedCount=29,assets=assets,oldObjectsStable=len(old),actualDependencies=len(inputs),
        sourceSnapshots=len(manifest),closedEntryOutput=result.stdout.decode('utf8',errors='replace'),processBudgetSecondsBeforeFinal=totals,
        peakObservedRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in processes),
        initialSourcePlanMemoryGuardMatched=False,initialSourcePlannedBytes=2*1024**3,initialSourceOnlineBytes=4*1024**3,
        nativeInterfaceControlAccepted=True,actualRestartBitwise=True,CFDIterations=0,newQualifiedAreaM2=0,
        physicalQualification=False,goalStatus='active',seconds=time.monotonic()-start)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',EV/'final_report.json')
    print(json.dumps(dict(passed=True,protectedCount=29,oldObjectsStable=len(old),sourceSnapshots=len(manifest))))
if __name__=='__main__':main()
