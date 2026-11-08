"""监督隔离实数/复数气体helper；不链接流动主程序、不启动CFD。"""
from pathlib import Path
import json,sys,shutil,time
from runEilmerProcess import run,sha,PREFIX,COMMIT

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v210_validation'
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def build():
    start=time.monotonic();out=BASE/'build04';out.mkdir(exist_ok=False);prep=BASE/'prepared02'
    shutil.copyfile(__file__,out/'executed_entry.py')
    r=json.loads((prep/'report.json').read_text());assert r['passed']
    assert sha(prep/'overlay/gas/rah_phase_aware_air.d')==r['nativeSourceSHA256']
    assert sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json')==r['dataSHA256']
    src=PREFIX/'source';compiler=PREFIX/'ldc2-1.41.0-linux-x86_64/bin/ldc2'
    lua=src/'extern/lua-5.4.3/install/lib/liblua.a';eqc=src/'src/extern/eqc/source/libeqc.a'
    files=[p for p in prep.rglob('*') if p.is_file()]+[Path(__file__).resolve(),ROOT/'docs/validation/v210/PLAN.md',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'docs/validation/v210/BUILD_RECOVERY.md',
        ROOT/'docs/validation/v210/HELPER_RECOVERY.md',compiler,lua,eqc]
    # 自动导入编译仍绑定全部可能的上游D源，避免只锁入口却丢依赖身份。
    files += sorted(p for p in (src/'src').rglob('*.d') if p.is_file())
    before={str(p):sha(p) for p in files}
    save(out/'admission.json',dict(sourceCommit=COMMIT,eachBuildSeconds=60,totalBuildSeconds=120,
        parentSeconds=150,memoryMiB=4096,outputMiB=96,cpuMaximum=2,automaticRetry=0,inputSHA256=before))
    shutil.copyfile(__file__,out/'executed_entry.py');results=[]
    for variant in ['real','complex']:
        command=[str(compiler),'-i','-O0','-g','-boundscheck=on','-dip1008','-preview=in',
            '-I'+str(prep/'overlay'),'-I'+str(src/'src'),'-I'+str(src/'src/extern/eqc/source'),'-J'+str(ROOT/'data/thermo'),
            '-of='+str(out/('probe_'+variant)),'-deps='+str(out/(variant+'_dependencies.txt')),
            str(prep/'probe.d'),str(eqc),str(lua),'-L-ldl']
        if variant=='complex':command.insert(1,'-d-version=complex_numbers')
        remaining=120-sum(row['elapsedSeconds'] for row in results)
        result=run(command,out,out/variant,min(60,remaining),96,files+[out/'admission.json'],4096,4096)
        results.append(result);save(out/'progress.json',results)
        if not result['passed']:break
    stable=all(sha(p)==h for p,h in before.items())
    report=dict(passed=len(results)==2 and all(x['passed'] for x in results) and stable,results=results,
        inputStable=stable,elapsedSeconds=time.monotonic()-start,engineBuilds=0,CFDRuns=0,physicalQualification=False)
    if report['passed']:report['binarySHA256']={v:sha(out/('probe_'+v)) for v in ['real','complex']}
    save(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='results'}))
    return report['passed']

def query():
    start=time.monotonic();out=BASE/'query02';out.mkdir(exist_ok=False);prep=BASE/'prepared02';built=BASE/'build04'
    r=json.loads((built/'report.json').read_text());assert r['passed']
    inputs=BASE/'query_inputs02'
    locks=[p for d in [prep,inputs] for p in d.rglob('*') if p.is_file()]+[Path(__file__).resolve(),ROOT/'docs/validation/v210/PLAN.md',ROOT/'docs/validation/v210/QUERY_RECOVERY.md',built/'report.json']
    shutil.copyfile(__file__,out/'executed_entry.py')
    save(out/'admission.json',dict(eachQuerySeconds=15,totalSeconds=30,memoryMiB=1024,outputMiB=16,automaticRetry=0))
    results=[]
    for variant in ['real','complex']:
        binary=built/('probe_'+variant);assert sha(binary)==r['binarySHA256'][variant]
        result=run([str(binary),str(prep/'gas.lua'),str(inputs/(variant+'.tsv')),str(prep/'invalid-gas.lua')],
            out,out/variant,15,16,locks+[binary,out/'admission.json'],1024,1024)
        results.append(result)
        if not result['passed']:break
    report=dict(passed=len(results)==2 and all(x['passed'] for x in results),results=results,
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,CFDRuns=0)
    save(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='results'}))
    return report['passed']

if __name__=='__main__':
    if sys.argv[1:]==['build']:ok=build()
    elif sys.argv[1:]==['query']:ok=query()
    else:raise ValueError('仅具名build/query；无CFD准入')
    raise SystemExit(0 if ok else 1)
