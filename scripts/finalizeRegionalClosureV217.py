"""模块职责：绑定真实观察与独立审核、保留监督限制差异、关闭入口并归档小证据。"""
from pathlib import Path
import ast,hashlib,json,re,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v217_validation';DOC=ROOT/'docs/validation/v217'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_bytes((json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf8'))

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    names=['inspectRegionalClosureV217.py','prepareRegionalObservationV217.py','auditRegionalStencilV217.py','runRegionalClosureV217.py']
    # 恢复初始调用入口的准确文本，只有匹配历史SHA才登记为已绑定。
    runner=(BASE/'entry_build/executed_entry.py').read_text(encoding='utf8')
    initial=runner.replace('import json,os,sys,shutil,time','import json,os,sys')
    initial=initial.replace("entries={'source':'inspectRegionalClosureV217.py','prepare':'prepareRegionalObservationV217.py'}", "entries={'source':'inspectRegionalClosureV217.py'}")
    a=initial.index("        if stage=='prepare':");b=initial.index('        used=',a);initial=initial[:a]+initial[b:]
    a=initial.index('    from runEilmerProcess import');b=initial.index("\nif __name__",a)
    initial=initial[:a]+"    raise ValueError('尚无原生阶段准入')\n"+initial[b:]
    source_receipt=json.loads((BASE/'offline_source/report.json').read_text())
    initial_expected=next(h for p,h in source_receipt['inputSHA256'].items() if p.replace('\\','/').endswith('/runRegionalClosureV217.py'))
    restored=hashlib.sha256(initial.encode('utf8')).hexdigest()==initial_expected
    if restored:(out/'source_entry_recovered.py').write_bytes(initial.encode('utf8'))
    entries={}
    for stage in ('source','prepare','audit'):
        record=BASE/('offline_'+stage);r=json.loads((record/'report.json').read_text())
        assert r['passed'] and r['inputsStable'] and r['nativeExitCode']==0 and not r['liveGroupAfterCleanup']
        task=record/'executed_task.py';assert sha(task) in r['inputSHA256'].values()
        # 这里副本属于被导入的Windows监督模块，准确角色记录而非冒充调用入口。
        assert sha(record/'executed_entry.py')==sha(ROOT/'scripts/runCompleteRegionDecisionV216.py')
        if stage=='source':matched=restored
        elif stage=='prepare':matched=sha(BASE/'entry_build/executed_entry.py') in r['inputSHA256'].values()
        else:matched=sha(ROOT/'scripts/runRegionalClosureV217.py') in r['inputSHA256'].values()
        entries[stage]=dict(elapsedSeconds=r['elapsedSeconds'],RSSBytes=r['peakObservedGroupRSSBytes'],callerSourceBound=matched)
    native={}
    for stage in ('build','query'):
        r=json.loads((BASE/(stage+'01')/'report.json').read_text());entry=BASE/('entry_'+stage)
        assert r['passed'] and r['nativeExitCode']==0 and r['inputsStable'] and r['stopReason'] is None and not r['liveGroupAfterCleanup']
        assert sha(entry/'executed_entry.py') in r['inputSHA256'].values()
        assert sha(entry/'executed_preparation.py') in r['inputSHA256'].values()
        native[stage]=dict(elapsedSeconds=r['elapsedSeconds'],RSSBytes=r['peakObservedGroupRSSBytes'],inputCount=len(r['inputSHA256']))
    b=json.loads((BASE/'build_result.json').read_text());assert sha(ROOT/'tools/runtime/eilmer-stencil-observation-v1/probe')==b['binarySHA256']
    audit=json.loads((BASE/'audit01/report.json').read_text())
    for p,h in audit['sourceSHA256'].items():assert sha(ROOT/p)==h
    assert audit['passed'] and audit['selectedOldOutputMaxScaledDifference']==0 and not audit['fullNumericalAcceptance']
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items())
    final_old=json.loads((ROOT/'outputs/v216_validation/final01/report.json').read_text())
    for p,h in final_old['sourceSHA256'].items():
        if p.endswith('.npz') or p.endswith('.mat'):assert sha(ROOT/p.replace('\\','/'))==h
    # 边界Lua/face动作另绑定准确固定源，尚未声称已实际分派新边界。
    old_manifest=json.loads((ROOT/'outputs/v212_validation/source01/report.json').read_text())
    lua=ROOT/'outputs/v212_validation/source01/source/src/lmr/lua-modules/bc.lua'
    assert sha(lua)==next(x['sha256'] for x in old_manifest['files'] if x['path']=='src/lmr/lua-modules/bc.lua')
    face=ROOT/'outputs/v209_validation/source02/source/src/lmr/bc/boundary_interface_effect.d'
    build=json.loads((BASE/'build01/report.json').read_text())
    assert sha(face)==build['inputSHA256']['/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/bc/boundary_interface_effect.d']
    source_rows=[]
    for p,ranges in [(lua,[(1260,1295),(1549,1572)]),(face,[(228,279),(281,340)])]:
        lines=p.read_text(encoding='utf8').splitlines()
        source_rows.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),excerpts=[dict(first=a,last=b,text='\n'.join('%d: %s'%(j+1,lines[j]) for j in range(a-1,b))) for a,b in ranges]))
    save(out/'boundary_source.json',dict(passed=True,sources=source_rows))
    # 先保存精确当前版本，再关闭CLI；只做启动前拒绝检查。
    closed=BASE/'closed_entries';closed.mkdir(exist_ok=False)
    for name in names:
        p=ROOT/'scripts'/name;shutil.copyfile(p,closed/name);text=p.read_text(encoding='utf8')
        match=re.search(r"^if __name__\s*==\s*['\"]__main__['\"]\s*:",text,re.M);assert match
        p.write_bytes((text[:match.start()]+"if __name__ == '__main__':\n    raise ValueError('v217观察与审核已结束；当前无构建/查询/网格/CFD准入。')\n").encode('utf8'))
        ast.parse(p.read_text(encoding='utf8'))
    refusal=subprocess.run([sys.executable,str(ROOT/'scripts/runRegionalClosureV217.py'),'query'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=5)
    assert refusal.returncode!=0 and b'ValueError' in refusal.stderr
    (out/'closed_entry_stderr.txt').write_bytes(refusal.stderr)
    save(out/'report.json',dict(passed=True,offline=entries,native=native,initialCallerRecoveredAndSHAMatched=restored,
        wrapperSnapshotRole='imported Windows supervisor, not v217 caller',
        windowsPreparationTighterStopNotEnforced=True,queryActualTasksetOneCPU=True,
        localPrimitiveTangentAccepted=True,fullNumericalAcceptance=False,physicalBoundaryAccepted=False,CFDAdmitted=False,
        protectedFilesUnchanged=29,originalMeshAndSevenPartitionUnchanged=True,CLIClosed=True,CLIRefusalPassed=True,
        CFDAdvanceRuns=0,newQualifiedAreaM2=0,goalComplete=False,bookkeepingSeconds=time.monotonic()-start))
    # 原生配置/网格留本地，只复制小源码、状态CSV、回执；不归档未结束任务。
    dest=DOC/'evidence';dest.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE);item=dict(source=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
        if rel.parts[0]=='control01' and p.name!='stencil.csv':local.append(item);continue
        assert p.stat().st_size<1024**2
        q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==item['sha256']
        item['path']=q.relative_to(ROOT).as_posix();rows.append(item)
    sources=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_LOCAL_STENCIL_EVIDENCE_V217.md']
    sources+=['scripts/'+name for name in names]+[Path(__file__).resolve().relative_to(ROOT).as_posix()]
    sources+=[p.relative_to(ROOT).as_posix() for p in sorted(DOC.glob('*.md'))]
    assert sum(r['bytes'] for r in rows)<4*1024**2
    save(DOC/'evidence_manifest.json',dict(schema='RAH.V217Evidence.v1',files=rows,localAssets=local,
         sourceFiles=sources,sourceSHA256={p:sha(ROOT/p) for p in sources},count=len(rows),bytes=sum(r['bytes'] for r in rows),
         binaryLocalPath='tools/runtime/eilmer-stencil-observation-v1/probe',binarySHA256=b['binarySHA256'],
         localBackupsCreated=False,protectedFilesUnchanged=29,fullNumericalAcceptance=False,CFDAdmitted=False,goalComplete=False))
    print(json.dumps(dict(passed=True,evidenceFiles=len(rows),evidenceBytes=sum(r['bytes'] for r in rows),initialCallerSourceBound=restored)))

if __name__=='__main__':main()
