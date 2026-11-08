"""模块职责：监督v214独立准备/构建/实际Jv控制，旧包与旧输入不覆盖。"""
from pathlib import Path
import json,sys,shutil,os,time
from runEilmerProcess import run,sha
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v214_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-frozen-nk-probe-v1';PLAN=ROOT/'docs/validation/v214/PLAN.md'

def prepare():
    BASE.mkdir(exist_ok=True)
    locks=[PLAN,Path(__file__).resolve(),ROOT/'scripts/prepareNativeNKControlV214.py',
        ROOT/'native/eilmerFrozenMapped/probe.d',ROOT/'outputs/v213_validation/build02/report.json']
    locks += [p for p in (ROOT/'outputs/v213_validation/control02').rglob('*') if p.is_file()]
    r=run([sys.executable,str(ROOT/'scripts/prepareNativeNKControlV214.py')],BASE,BASE/'preparation01',60,32,locks,1024,1024)
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'])));return r['passed']

def build():
    start=time.monotonic();assert json.loads((BASE/'preparation01/report.json').read_text())['passed']
    PACKAGE.mkdir(exist_ok=False);old=json.loads((ROOT/'outputs/v213_validation/build02/report.json').read_text())
    oldPackage=ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v2'
    shutil.copyfile(oldPackage/'buildinfo.json',PACKAGE/'buildinfo.json')
    new=BASE/'prepared01/nk_probe.d';oldProbe=ROOT/'native/eilmerFrozenMapped/probe.d'
    cmd=[str(new) if p==str(oldProbe) else p.replace(str(oldPackage),str(PACKAGE)) for p in old['argv']]
    locks=[]
    for p,h in old['inputSHA256'].items():
        if Path(p).suffix not in ('.d','.di','.a','.json'):continue
        assert sha(p)==h;locks.append(p)
    locks += [new,PLAN,Path(__file__).resolve(),ROOT/'scripts/prepareNativeNKControlV214.py',PACKAGE/'buildinfo.json']
    r=run(cmd,PACKAGE,BASE/'build01',150,512,locks,6144,6144)
    result=dict(passed=r['passed'],seconds=r['elapsedSeconds'],entrySeconds=time.monotonic()-start)
    if r['passed']:result['binarySHA256']=sha(PACKAGE/'probe')
    (BASE/'build_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return r['passed']

def query():
    start=time.monotonic();b=json.loads((BASE/'build_result.json').read_text());assert b['passed']
    binary=PACKAGE/'probe';assert sha(binary)==b['binarySHA256'];case=BASE/'control01'
    locks=[p for p in case.rglob('*') if p.is_file()]+[binary,PLAN,Path(__file__).resolve()]
    cmd=['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary)]
    r=run(cmd,case,BASE/'query01',60,32,locks,2048,2048)
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'],entrySeconds=time.monotonic()-start,exit=r['nativeExitCode'])));return r['passed']

if __name__=='__main__':
    stage=sys.argv[1]
    if stage=='prepare':ok=prepare()
    elif stage=='build':ok=build()
    elif stage=='query':ok=query()
    else:raise ValueError('stage')
    raise SystemExit(0 if ok else 1)
