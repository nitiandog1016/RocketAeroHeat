"""模块职责：复核实际文件往返、历史保护和源绑定，关闭本轮入口并归档小证据。"""
from pathlib import Path
import json,hashlib,shutil,sys,ast,subprocess,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v215_validation';DEST=ROOT/'docs/validation/v215/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf8')

def entry_versions():
    """模块一：依据精确编辑差异恢复未即时保存的入口文本，必须与执行回执SHA相同。"""
    runner=(ROOT/'scripts/runNativeFileBridgeV215.py').read_text(encoding='utf8')
    nofinal=runner.replace("elif stage in ('audit','final-audit'):","elif stage=='audit':")
    nofinal=nofinal.replace("        if stage=='final-audit':remaining-=json.loads((BASE/'audit01/report.json').read_text())['elapsedSeconds']\n",'')
    nofinal=nofinal.replace("'final_audit01' if stage=='final-audit' else 'audit01'","'audit01'")
    a=nofinal.index("    elif stage in ('export','export-recover','import'):")
    b=nofinal.index("    elif stage=='audit':",a)
    oldquery='''    elif stage in ('export','import'):
        assert json.loads((BASE/'build_result.json').read_text())['passed'];binary=PACKAGE/'probe'
        assert sha(binary)==json.loads((BASE/'build_result.json').read_text())['binarySHA256']
        case=BASE/'control01';locks=common+[binary]+[p for p in case.rglob('*') if p.is_file()]
        remaining=60
        if stage=='import':remaining-=json.loads((BASE/'export01/report.json').read_text())['elapsedSeconds']
        r=invoke(['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary),stage],case,stage+'01',remaining,2048,64,locks)
'''
    noquery=nofinal[:a]+oldquery+nofinal[b:]
    a=noquery.index("    elif stage=='recover':");b=noquery.index("    elif stage=='build':",a)
    initial=noquery[:a]+noquery[b:]
    prep=(ROOT/'scripts/prepareNativeFileBridgeV215.py').read_text(encoding='utf8')
    gas=prep.replace("'mesh.su2','phase-air.lua'","'mesh.su2','gas-model.lua'")
    original=gas.replace('auto receipt=File("export.json","w");receipt.writefln','auto out=File("export.json","w");out.writefln')
    original=original.replace('InterpolateOption.pt;','InterpolateOption.pT;').replace('auto expected=parseJSON(meta.toString());','auto expected=meta;')
    texts={'runner-initial.py':initial,'runner-compile-recovery.py':noquery,'runner-file-recovery.py':nofinal,
           'runner-final.py':runner,'prepare-initial.py':original,'prepare-syntax-recovery.py':gas,'prepare-final.py':prep}
    dest=BASE/'executed_entries';dest.mkdir(exist_ok=False)
    for name,text in texts.items():(dest/name).write_text(text,encoding='utf8')
    return {sha(p):p for p in dest.iterdir()}

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    versions=entry_versions();stages=['preparation01','build01','build02','export01','export02','audit01','import01','final_audit01']
    processes={};entries={}
    for name in stages:
        r=json.loads((BASE/name/'report.json').read_text());expected=name not in ('build01','export01')
        assert r['inputsStable'] and r['stopReason'] is None and not r['liveGroupAfterCleanup']
        assert r['passed']==expected and r['nativeExitCode']==(0 if expected else 1)
        for path,h in r['inputSHA256'].items():
            if Path(path).name in ('runNativeFileBridgeV215.py','prepareNativeFileBridgeV215.py'):
                assert h in versions,(name,path,h);entries[name+'/'+Path(path).name]=str(versions[h].relative_to(ROOT))
        processes[name]={k:r[k] for k in ['passed','nativeExitCode','elapsedSeconds','peakObservedGroupRSSBytes','inputsStable']}
    build=json.loads((BASE/'build02/report.json').read_text())
    for rel in ['native/eilmerFrozenMapped/lmr/bc/ghost_cell_effect/rah_frozen_file.d',
                'native/eilmerFrozenMapped/lmr/bc/ghost_cell_effect/rah_frozen_mapped.d',
                'outputs/v215_validation/prepared02/file_probe.d']:
        p=ROOT/rel;assert sha(p)==build['inputSHA256'][str(p)]
    binary=ROOT/'tools/runtime/eilmer-frozen-file-probe-v1/probe'
    assert sha(binary)==json.loads((BASE/'build_result.json').read_text())['binarySHA256']
    for name in ['audit01','final_audit01']:
        record=json.loads((BASE/name/'report.json').read_text())
        for rel in ['tests/testNativeFrozenFile.py','scripts/nativeFrozenSnapshot.py']:
            p=ROOT/rel;assert sha(p)==record['inputSHA256'][str(p)]
    result=json.loads((BASE/'control01/final-audit.json').read_text())
    assert result['passed'] and result['bitwiseNativeRoundtrip'] and result['nativeRejectedCount']==56
    assert not result['fullNumericalAcceptance'] and result['CFDAdvanceRuns']==0
    old=json.loads((ROOT/'docs/validation/v212/evidence/packet_work02/audit/report.json').read_text())
    assert old['passed'] and old['rejectedCount']==48
    for p,h in old['inputSHA256'].items():assert sha(ROOT/p)==h
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/k.replace('\\','/'))==v for k,v in protected.items())
    # 模块二：关闭已结束的原生任务入口；应用接口/独立读取器仍可复用。
    closed=BASE/'closed_entries';closed.mkdir(exist_ok=False)
    for name in ['runNativeFileBridgeV215.py','prepareNativeFileBridgeV215.py']:
        p=ROOT/'scripts'/name;shutil.copyfile(p,closed/name)
        text=p.read_text(encoding='utf8');index=text.rindex("if __name__=='__main__':")
        p.write_text(text[:index]+"if __name__=='__main__':raise ValueError('v215文件桥接验证已结束；当前无新构建/查询/CFD准入。')\n",encoding='utf8')
        r=subprocess.run([sys.executable,str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        assert r.returncode!=0 and b'ValueError' in r.stderr
    report=dict(passed=True,processes=processes,executedEntrySnapshots=entries,
        actualSourceBound=True,v1Unchanged=True,v1ExistingRejections=48,protectedFilesUnchanged=29,
        nativeFileBridgePassed=True,fullNumericalAcceptance=False,CFDAdvanceRuns=0,newQualifiedAreaM2=0,
        goalComplete=False,uniqueDeleted=0,seconds=time.monotonic()-start)
    save(out/'report.json',report)
    # 模块三：完整数据、网格和包留E盘；公开小回执、准确helper、源码与SHA清单。
    DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE).as_posix();item=dict(source=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
        keep=('/lmrsim/' not in rel and p.name not in ('native-snapshot.json','reordered-snapshot.json','roundtrip.json',
            'reordered-roundtrip.json','requests.json','negative-current.json','mapping.npz','mesh.su2','donors.json'))
        if not keep:local.append(item);continue
        assert p.stat().st_size<1024**2
        q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==sha(p)
        item['path']=q.relative_to(ROOT).as_posix();rows.append(item)
    sourceFiles=['AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_NATIVE_FILE_BRIDGE_V215.md',
        'native/eilmerFrozenMapped/lmr/bc/ghost_cell_effect/rah_frozen_file.d',
        'scripts/nativeFrozenSnapshot.py','scripts/prepareNativeFileBridgeV215.py','scripts/runNativeFileBridgeV215.py',
        'scripts/finalizeNativeFileBridgeV215.py','tests/testNativeFrozenFile.py']
    sourceFiles += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v215').glob('*.md')]
    for p in sourceFiles:
        if p.endswith('.py'):ast.parse((ROOT/p).read_text(encoding='utf8'))
    assert sum(r['bytes'] for r in rows)<2*1024**2
    save(ROOT/'docs/validation/v215/evidence_manifest.json',dict(schema='RAH.V215Evidence.v1',files=rows,localAssets=local,
        sourceFiles=sourceFiles,sourceSHA256={p:sha(ROOT/p) for p in sourceFiles},count=len(rows),bytes=sum(r['bytes'] for r in rows),
        nativeFileBridgePassed=True,fullNumericalAcceptance=False,protectedFilesUnchanged=29,goalComplete=False))
    print(json.dumps(dict(passed=True,evidenceFiles=len(rows),evidenceBytes=sum(r['bytes'] for r in rows))))

if __name__=='__main__':main()
