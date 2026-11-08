"""模块职责：核对既存原生/离线证据、关闭入口和共享保护，归档小证据；不重算流动。"""
from pathlib import Path
import json,hashlib,shutil,ast,time,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'
DEST=ROOT/'docs/validation/v213/evidence'
NAMES=['inspectFrozenMappedEntryV213.py','prepareFrozenMappedControlV213.py','runFrozenMappedControlV213.py',
       'recoverFrozenMappedGeometryV213.py','discriminateFrozenOperatorV213.py']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(),out/'executed_entry.py')
    audit=json.loads((BASE/'operator_audit01/report.json').read_text());assert audit['auditPassed'] and not audit['nativeAcceptance']
    assert len(audit['regimes'])==3 and all(r['baseMax']==r['derivativeMax']==0 for r in audit['regimes'])
    for p,h in audit['sourceSHA256'].items():assert sha(ROOT/p)==h
    corrected=json.loads((BASE/'derivative_mismatch_corrected.json').read_text())
    original=json.loads((BASE/'derivative_mismatch.json').read_text());assert original['records']==corrected['records']
    assert corrected['sourceSHA256']==sha(BASE/'control03/samples.csv')
    processes={}
    for p in sorted(BASE.glob('*/report.json')):
        r=json.loads(p.read_text())
        if r.get('schema')!='RAH.EilmerProcess.v1':continue
        assert r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None
        processes[p.parent.name]={k:r[k] for k in ['nativeExitCode','passed','elapsedSeconds','peakObservedGroupRSSBytes','inputsStable']}
    assert len(processes)==11
    assert [processes['query%02d'%i]['nativeExitCode'] for i in (1,2,3)]==[1,1,1]
    assert [processes['build%02d'%i]['passed'] for i in (1,2)]==[True,True]
    assert sum(processes[n]['elapsedSeconds'] for n in ['build01','build02'])<240
    assert sum(processes[n]['elapsedSeconds'] for n in ['query01','query02','query03'])<60
    for i in (1,2,3):
        log=(BASE/('query%02d'%i)/'stdout.txt').read_text()
        assert ('RAH_OPERATOR_GATE' if i>1 else 'failed to initialise geometry') in log
    # 原始构建与各准确源码快照绑定；旧脚本修订不当作原始运行版本。
    b1=json.loads((BASE/'build01/report.json').read_text());b2=json.loads((BASE/'build02/report.json').read_text())
    key='/mnt/e/RocketAeroHeat/native/eilmerFrozenMapped/probe.d'
    assert sha(BASE/'query02/executed_source/probe.d')==b1['inputSHA256'][key]
    assert sha(ROOT/'native/eilmerFrozenMapped/probe.d')==b2['inputSHA256'][key]==sha(BASE/'query03/executed_probe.d')
    for package,receipt in [('v1','build_result.json'),('v2','discriminator_build.json')]:
        h=json.loads((BASE/receipt).read_text())['binarySHA256'];folder=ROOT/('tools/runtime/eilmer-frozen-mapped-probe-'+package)
        assert sha(folder/'probe')==h==sha(folder/'bin/probe')
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text());assert len(protected)==29
    assert all(sha(ROOT/k.replace('\\','/'))==h for k,h in protected.items())
    closed=[]
    for name in NAMES:
        p=ROOT/'scripts'/name;tree=ast.parse(p.read_text(encoding='utf8'));assert isinstance(tree.body[-1].body[0],ast.Raise)
        r=subprocess.run([sys.executable,str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        assert r.returncode!=0 and b'ValueError' in r.stderr
        closed.append(dict(file=name,exitCode=r.returncode))
    save(out/'report.json',dict(passed=True,nativeAcceptance=False,physicalQualification=False,CFDAdvanceRuns=0,
        newQualifiedAreaM2=0,protectedFilesUnchanged=29,closedEntries=closed,processes=processes,
        conversionResourceAcceptance=False,independentAuditSHA256=sha(BASE/'operator_audit01/report.json'),
        firstNativeFailurePreserved=True,originalNumericalFailurePreserved=True,uniqueDeleted=0,goalComplete=False,
        seconds=time.monotonic()-start))
    # 保存定向源摘录；固定上游的完整副本不重复上传。
    selections=[
      (ROOT/'outputs/v209_validation/source04/source/src/geom/elements/properties.d',[(290,305),(538,579)]),
      (BASE/'source02/source/src/geom/grid/usgrid.d',[(1440,1456)]),
      (BASE/'source02/source/src/lmr/newtonkrylovsolver.d',[(3359,3381),(3422,3477),(3517,3588)]),
      (ROOT/'outputs/v212_validation/source01/source/src/lmr/simcore_exchange.d',[(107,135),(152,180)]),
      (ROOT/'outputs/v211_validation/source01/source/lmr/lsqinterp.d',[(601,652)])]
    excerpts=[]
    for p,ranges in selections:
        lines=p.read_text(encoding='utf8').splitlines()
        excerpts.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),
            ranges=[dict(start=a,end=b,text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
    save(BASE/'source_excerpts.json',excerpts)
    # 每个执行结果和小输入逐字归档；完整逐量CSV及原生场留E盘并登记SHA。
    DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(BASE).as_posix()
        excluded=('/source/' in rel or '/lmrsim/' in rel or p.name=='samples.csv')
        record=dict(source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
        if excluded:local.append(record);continue
        assert p.stat().st_size<1024**2
        q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
        record['path']=q.relative_to(ROOT).as_posix();rows.append(record)
    sourceFiles=['.gitattributes','AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_FROZEN_NATIVE_OPERATOR_V213.md']
    sourceFiles += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'native/eilmerFrozenMapped').rglob('*')) if p.is_file()]
    sourceFiles += ['scripts/'+n for n in NAMES+['auditFrozenMappedGeometryV213.py','auditFrozenOperatorV213.py',
                    'analyzeFrozenDerivativeMismatchV213.py','finalizeFrozenNativeV213.py']]
    sourceFiles += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT/'docs/validation/v213').glob('*.md'))]
    for name in sourceFiles:
        p=ROOT/name
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    size=sum(r['bytes'] for r in rows);assert size<6*1024**2
    save(ROOT/'docs/validation/v213/evidence_manifest.json',dict(schema='RAH.V213Evidence.v1',files=rows,localAssets=local,
        count=len(rows),bytes=size,sourceFiles=sourceFiles,nativeAcceptance=False,goalComplete=False,protectedFilesUnchanged=29))
    print(json.dumps(dict(passed=True,count=len(rows),bytes=size,localFiles=len(local),seconds=time.monotonic()-start)))

if __name__=='__main__':main()
