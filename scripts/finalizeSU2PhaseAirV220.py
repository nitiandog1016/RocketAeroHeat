"""模块职责：保留全部小证据/准确历史入口，核对保护资产并关闭本轮准入。"""
from pathlib import Path
import ast,hashlib,json,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v220_validation'
DOC=ROOT/'docs/validation/v220';EV=DOC/'evidence'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);EV.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    assets={
        'outputs/v144_validation/actual/hybrid_volume.mat':'22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c',
        'outputs/v167_validation/real01/mesh.su2':'00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181',
        'data/thermo/phase_aware_dry_air_v1.json':'7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'}
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    assert read(BASE/'audit01/report.json')['passed']
    scripts=['prepareSU2PhaseAirV220.py','buildSU2PhaseAirV220.py','runSU2PhaseAirV220.py',
             'recoverSU2PhaseAirV220.py','verifySU2PhaseAirV220.py','finalizeSU2PhaseAirV220.py']
    for n in scripts:ast.parse((ROOT/'scripts'/n).read_text(encoding='utf8'))
    # 两个历史执行入口均应在创建进程或目录之前拒绝。
    closed=[]
    for n in ('runSU2PhaseAirV220.py','recoverSU2PhaseAirV220.py'):
        result=subprocess.run([sys.executable,str(ROOT/'scripts'/n),'query'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        assert result.returncode!=0 and b'closed' in result.stdout
        closed.append(dict(script=n,exitCode=result.returncode,output=result.stdout.decode('utf8',errors='replace')))
    processes=[]
    for p in sorted(BASE.glob('process_*/report.json')):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup']
        shutil.copytree(p.parent,EV/p.parent.name);processes.append(r)
    # 只带必要证据；大源树、对象、引擎和历史资产留E盘。
    for n in ('prepare01','prepare02','prepare03','dependencies01','link01','audit01','control01','control02','control03',
              'helper01_source','helper02_source','observe01'):
        shutil.copytree(BASE/n,EV/n)
    (EV/'build02').mkdir()
    for n in ('dependencies.json','members.json','link_command.json','link_report.json'):
        shutil.copyfile(BASE/'build02'/n,EV/'build02'/n)
    (EV/'generated').mkdir()
    shutil.copyfile(BASE/'source/SU2_CFD/include/fluid/PhaseAirTable.hpp',EV/'generated/PhaseAirTable.hpp')
    for n in ('COPYING','AUTHORS.md'):shutil.copyfile(ROOT/'docs/validation/v219/evidence'/n,EV/n)
    # 按执行SHA恢复历次入口，而不是把最终修订版冒称旧执行源。
    candidates={}
    def candidate(data):candidates[hashlib.sha256(data).hexdigest()]=data
    def variants(text):
        text=text.replace('\r\n','\n');candidate(text.encode('utf8'));candidate(text.replace('\n','\r\n').encode('utf8'))
    for p in EV.rglob('*'):
        if p.is_file() and p.suffix in ('.py','.cpp','.hpp','.md'):candidate(p.read_bytes())
    own=[ROOT/'scripts'/n for n in scripts]+list((ROOT/'native/su2PhaseAwareAir').glob('*'))+[DOC/'PLAN.md']
    for p in own:variants(p.read_text(encoding='utf8'))
    plan=(DOC/'PLAN.md').read_text(encoding='utf8')
    lines=plan.splitlines(True)
    for i in range(1,len(lines)+1):variants(''.join(lines[:i]).rstrip()+'\n')
    recovery=(ROOT/'scripts/recoverSU2PhaseAirV220.py').read_text(encoding='utf8')
    recovery=recovery.replace('EXECUTION_CLOSED=True\n','').replace("    if EXECUTION_CLOSED:raise SystemExit('v220 recovery closed; new control requires a separate recorded plan')\n",'')
    variants(recovery)
    recovery=recovery.replace("'recover_prepare02'","'recover_prepare'").replace('str(Path(__file__).resolve())','__file__')
    variants(recovery)
    wanted={}
    for r in processes:
        for n,h in r['inputSHA256'].items():
            p=Path(n)
            if (p.parent==ROOT/'scripts' and p.name in scripts) or p.parent==ROOT/'native/su2PhaseAwareAir' or p==DOC/'PLAN.md':
                wanted.setdefault(h,[]).append(n)
    snapshots=EV/'execution_sources';snapshots.mkdir();snapshot_rows=[];missing=[]
    for h,names in wanted.items():
        data=candidates.get(h)
        if data is None:missing.append(dict(sha256=h,paths=sorted(set(names))));continue
        name=h+Path(names[0]).suffix;(snapshots/name).write_bytes(data)
        snapshot_rows.append(dict(sha256=h,paths=sorted(set(names)),snapshot=name))
    save(EV/'source_snapshot_manifest.json',snapshot_rows);save(EV/'missing_source_snapshots.json',missing)
    assert not missing,missing
    # 只快照本轮最终文件；.gitattributes改名以免其规则递归改变证据字节。
    final_files=own+[ROOT/'AGENTS.md',ROOT/'README.md',ROOT/'docs/PROJECT_GOALS.md',ROOT/'.gitattributes',
                    ROOT/'docs/HYBRID_SU2_PHASE_AIR_V220.md',DOC/'NEXT_DECISION.md']
    for p in final_files:
        relative=p.relative_to(ROOT)
        if p.name=='.gitattributes':relative=Path('.gitattributes.snapshot')
        dest=EV/'final_source'/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    totals={k:sum(r['elapsedSeconds'] for r in processes if r['budgetKind']==k) for k in ('offline','build','query')}
    assert totals['offline']<=120 and totals['build']<=480 and totals['query']<=60
    report=dict(passed=True,protectedFiles=len(protected),assets=assets,sourceSnapshots=len(snapshot_rows),
        allVersionedInputSourceBytesRecovered=True,closedEntries=closed,processCount=len(processes),
        processBudgetSecondsBeforeFinal=totals,peakObservedRSS=max(r['peakObservedGroupRSSBytes'] for r in processes),
        failedProcesses=[dict(stage=Path(r['argv'][0]).name,exitCode=r['exitCode'],seconds=r['elapsedSeconds'],stopReason=r['stopReason']) for r in processes if not r['passed']],
        completeNativeCompatibility=False,newEngineStarted=False,CFDRuns=0,newQualifiedAreaM2=0,goalStatus='active',seconds=time.monotonic()-start)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',EV/'final_report.json')
    manifest=[dict(path=p.relative_to(EV).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(EV.rglob('*')) if p.is_file()]
    save(EV/'manifest.json',manifest)
    print(json.dumps(dict(passed=True,protectedFiles=29,sourceSnapshots=len(snapshot_rows),evidenceBytes=sum(r['bytes'] for r in manifest))))

if __name__=='__main__':main()
