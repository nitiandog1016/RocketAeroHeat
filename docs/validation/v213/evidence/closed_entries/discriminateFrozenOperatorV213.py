"""模块职责：一次分量级原生对照；保持前次二进制与结果，不改变算子。"""
from pathlib import Path
import json,sys,shutil,os,time
from runEilmerProcess import run,sha,PREFIX
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v2'
PLAN=ROOT/'docs/validation/v213/OPERATOR_DISCRIMINATOR.md'

def build():
    start=time.monotonic();PACKAGE.mkdir(exist_ok=True)
    assert sorted(p.name for p in PACKAGE.iterdir())==['buildinfo.json']
    old=json.loads((BASE/'build01/report.json').read_text());assert old['passed']
    shutil.copyfile(ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v1/buildinfo.json',PACKAGE/'buildinfo.json')
    old_package=str(ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v1')
    cmd=[x.replace(old_package,str(PACKAGE)) for x in old['argv']]
    native=str(ROOT/'native/eilmerFrozenMapped')
    locks=[]
    for p,h in old['inputSHA256'].items():
        if p.startswith(native):continue
        if p==str(ROOT/'scripts/runFrozenMappedControlV213.py'):
            # 首版包装器已因prepsim包布局修正；使用已保存且SHA相同的准确执行副本。
            saved=BASE/'build01/executed_entry.py';assert sha(saved)==h;locks.append(saved);continue
        assert sha(p)==h
        locks.append(p)
    locks += [p for p in (ROOT/'native/eilmerFrozenMapped').rglob('*') if p.is_file()]
    locks += [PLAN,Path(__file__).resolve(),PACKAGE/'buildinfo.json']
    spent=old['elapsedSeconds'];assert spent<240
    r=run(cmd,PACKAGE,BASE/'build02',240-spent,512,locks,6144,6144)
    receipt=dict(passed=r['passed'],entrySeconds=time.monotonic()-start,cumulativeBuildSeconds=spent+r['elapsedSeconds'])
    if r['passed']:
        (PACKAGE/'bin').mkdir();(PACKAGE/'lib').mkdir();os.link(PACKAGE/'probe',PACKAGE/'bin/probe')
        shutil.copyfile(PREFIX/'install/lib/prepsim.lua',PACKAGE/'lib/prepsim.lua')
        receipt['binarySHA256']=sha(PACKAGE/'probe')
    (BASE/'discriminator_build.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
    return r['passed']

def query():
    build=json.loads((BASE/'discriminator_build.json').read_text());binary=PACKAGE/'bin/probe'
    assert build['passed'] and sha(binary)==build['binarySHA256']
    old=BASE/'control02';case=BASE/'control03';case.mkdir(exist_ok=False)
    # 复用已有准备文件；只排除原失败输出，旧目录逐字保持。
    for p in old.iterdir():
        if p.name=='operator_results.jsonl':continue
        if p.is_dir():shutil.copytree(p,case/p.name)
        else:shutil.copyfile(p,case/p.name)
    locks=[p for p in case.rglob('*') if p.is_file()]+[binary,PLAN,Path(__file__).resolve()]
    spent=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob('query*/report.json'));assert spent<60
    cmd=['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary)]
    r=run(cmd,case,BASE/'query03',60-spent,32,locks,2048,2048)
    shutil.copyfile(ROOT/'native/eilmerFrozenMapped/probe.d',BASE/'query03/executed_probe.d')
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'],exit=r['nativeExitCode'],cumulativeSeconds=spent+r['elapsedSeconds'])))
    return r['passed']

if __name__=='__main__':
    stage=sys.argv[1]
    if stage=='build':raise SystemExit(0 if build() else 1)
    elif stage=='query':raise SystemExit(0 if query() else 1)
    else:raise ValueError('stage')
