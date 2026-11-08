"""模块职责：核对已保存的正式NK与独立解析证据，归档交付；不重跑原生计算。"""
from pathlib import Path
import json,hashlib,shutil,ast,time,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v214_validation';DEST=ROOT/'docs/validation/v214/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(),out/'executed_entry.py')
    audit=json.loads((BASE/'audit01/report.json').read_text());assembly=json.loads((BASE/'assembly01/report.json').read_text())
    assert audit['auditPassed'] and not audit['nativeAcceptance'] and audit['directionRelativeError']<1e-11
    assert all(r['equivalenceError']==0 and r['zedError']==0 and not r['passed'] for r in audit['regimes'])
    assert all(r['groups']['integratedFlux']['failed']==0 for r in audit['regimes'])
    assert assembly['auditCompleted'] and not assembly['originalFiniteDifferenceAcceptance']
    assert all(r['derivativeGatePassed'] and r['valueGatePassed'] for r in assembly['records'])
    for p,h in audit['sourceSHA256'].items():assert sha(ROOT/p)==h
    assert assembly['inputSHA256']==sha(BASE/'control01/nk_samples.csv')
    processes={}
    for name in ['preparation01','build01','query01']:
        r=json.loads((BASE/name/'report.json').read_text())
        assert r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None
        assert r['passed']==(name!='query01') and r['nativeExitCode']==(1 if name=='query01' else 0)
        assert r['elapsedSeconds']<r['limitSeconds']
        processes[name]={k:r[k] for k in ['passed','nativeExitCode','elapsedSeconds','peakObservedGroupRSSBytes','inputsStable']}
    log=(BASE/'query01/stdout.txt').read_text();assert 'Single process running with 1 threads.' in log and 'RAH_NK_DERIVATIVE_GATE' in log
    build=json.loads((BASE/'build01/report.json').read_text());probe=BASE/'prepared01/nk_probe.d'
    assert sha(probe)==build['inputSHA256'][str(probe)]
    receipt=json.loads((BASE/'build_result.json').read_text());assert sha(ROOT/'tools/runtime/eilmer-frozen-nk-probe-v1/probe')==receipt['binarySHA256']
    for name in ['prepareNativeNKControlV214.py','runNativeNKControlV214.py']:
        assert sha(BASE/'closed_entries'/name)==build['inputSHA256'][str(ROOT/'scripts'/name)]
        p=ROOT/'scripts'/name;tree=ast.parse(p.read_text(encoding='utf8'));assert isinstance(tree.body[-1].body[0],ast.Raise)
        r=subprocess.run([sys.executable,str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE);assert r.returncode!=0 and b'ValueError' in r.stderr
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text());assert len(protected)==29
    assert all(sha(ROOT/k.replace('\\','/'))==h for k,h in protected.items())
    save(out/'report.json',dict(passed=True,processes=processes,protectedFilesUnchanged=29,nativeAcceptance=False,
        originalFailurePreserved=True,SAAssemblyVerified=True,independentFullJacobianVerified=False,
        offlineAssemblyRSSMeasured=False,offlineAssemblyResourceAcceptance=False,CFDAdvanceRuns=0,newQualifiedAreaM2=0,
        goalComplete=False,uniqueDeleted=0,seconds=time.monotonic()-start))
    selections=[
        (ROOT/'outputs/v213_validation/source02/source/src/lmr/newtonkrylovsolver.d',[(3304,3357),(3518,3650)]),
        (ROOT/'outputs/v211_validation/source01/source/lmr/turbulence.d',[(628,677),(780,810),(826,855),(1060,1103)])]
    excerpts=[]
    for p,ranges in selections:
        lines=p.read_text(encoding='utf8').splitlines();excerpts.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),
            ranges=[dict(start=a,end=b,text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
    save(BASE/'source_excerpts.json',excerpts)
    DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE).as_posix();record=dict(source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
        if '/lmrsim/' in rel or p.name=='nk_samples.csv':local.append(record);continue
        assert p.stat().st_size<1024**2
        q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==sha(p)
        record['path']=q.relative_to(ROOT).as_posix();rows.append(record)
    sourceFiles=['.gitattributes','AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_NATIVE_NK_OPERATOR_V214.md']
    sourceFiles+=['scripts/'+n for n in ['prepareNativeNKControlV214.py','runNativeNKControlV214.py',
        'auditNativeNKControlV214.py','auditNativeNKAssemblyV214.py','finalizeNativeNKControlV214.py']]
    sourceFiles += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'docs/validation/v214').glob('*.md'))]
    for name in sourceFiles:
        p=ROOT/name
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    size=sum(r['bytes'] for r in rows);assert size<3*1024**2
    save(ROOT/'docs/validation/v214/evidence_manifest.json',dict(schema='RAH.V214Evidence.v1',files=rows,localAssets=local,
        count=len(rows),bytes=size,sourceFiles=sourceFiles,protectedFilesUnchanged=29,nativeAcceptance=False,goalComplete=False))
    print(json.dumps(dict(passed=True,files=len(rows),bytes=size,seconds=time.monotonic()-start)))

if __name__=='__main__':main()
