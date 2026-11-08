"""模块职责：独立监督文件桥接的准备/构建/导出/导入；没有 CFD 发展入口。"""
from pathlib import Path
import json,sys,shutil,os,time
from runEilmerProcess import run,sha
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v215_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-frozen-file-probe-v1';PLAN=ROOT/'docs/validation/v215/PLAN.md'
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def invoke(cmd,cwd,name,seconds,memory,output,locks):
    start=time.monotonic();r=run(cmd,cwd,BASE/name,seconds,output,locks,memory,memory)
    print(json.dumps(dict(stage=name,passed=r['passed'],seconds=r['elapsedSeconds'],entrySeconds=time.monotonic()-start)))
    return r
def main(stage):
    BASE.mkdir(exist_ok=True)
    common=[PLAN,Path(__file__).resolve(),ROOT/'scripts/prepareNativeFileBridgeV215.py',
        ROOT/'native/eilmerFrozenMapped/lmr/bc/ghost_cell_effect/rah_frozen_file.d']
    if stage=='prepare':
        locks=common+[ROOT/'outputs/v214_validation/build01/report.json',ROOT/'native/eilmerFrozenMapped/probe.d']
        locks += [p for p in (ROOT/'outputs/v213_validation/control02').rglob('*') if p.is_file()]
        r=invoke([sys.executable,str(common[2])],BASE,'preparation01',60,1024,64,locks)
    elif stage=='recover':
        failed=json.loads((BASE/'build01/report.json').read_text());assert not failed['passed']
        out=BASE/'prepared02';out.mkdir(exist_ok=False)
        old=BASE/'prepared01/file_probe.d';new=out/'file_probe.d'
        text=old.read_text(encoding='utf8').replace('auto out=File("export.json","w");out.writefln','auto receipt=File("export.json","w");receipt.writefln')
        text=text.replace('InterpolateOption.pT','InterpolateOption.pt').replace('auto expected=meta;','auto expected=parseJSON(meta.toString());')
        new.write_text(text,encoding='utf8')
        template=json.loads((BASE/'prepared01/binding-template.json').read_text());template['sourceRecipeSHA256']=sha(new)
        save(out/'binding-template.json',template)
        cmd=[str(new) if p==str(old) else p for p in failed['argv']]
        locks=[p for p in failed['inputSHA256'] if Path(p).suffix in ('.d','.di','.a','.json')]+common+[new]
        r=invoke(cmd,PACKAGE,'build02',150-failed['elapsedSeconds'],6144,512,locks)
        save(BASE/'build_result.json',dict(passed=r['passed'],binarySHA256=sha(PACKAGE/'probe') if r['passed'] else None))
        if r['passed']:
            template['nativeBuildSHA256']=sha(PACKAGE/'probe');save(BASE/'control01/binding.json',template)
    elif stage=='build':
        assert json.loads((BASE/'preparation01/report.json').read_text())['passed'];PACKAGE.mkdir(exist_ok=False)
        old=json.loads((ROOT/'outputs/v214_validation/build01/report.json').read_text())
        oldpkg=ROOT/'tools/runtime/eilmer-frozen-nk-probe-v1'
        shutil.copyfile(oldpkg/'buildinfo.json',PACKAGE/'buildinfo.json')
        oldmain=ROOT/'outputs/v214_validation/prepared01/nk_probe.d';newmain=BASE/'prepared01/file_probe.d'
        cmd=[str(newmain) if p==str(oldmain) else p.replace(str(oldpkg),str(PACKAGE)) for p in old['argv']]
        cmd.append(str(common[3]));locks=[]
        for p,h in old['inputSHA256'].items():
            if Path(p).suffix not in ('.d','.di','.a','.json'):continue
            assert sha(p)==h;locks.append(p)
        locks += common+[newmain,PACKAGE/'buildinfo.json']
        r=invoke(cmd,PACKAGE,'build01',150,6144,512,locks)
        save(BASE/'build_result.json',dict(passed=r['passed'],binarySHA256=sha(PACKAGE/'probe') if r['passed'] else None))
        if r['passed']:
            m=json.loads((BASE/'prepared01/binding-template.json').read_text());m['nativeBuildSHA256']=sha(PACKAGE/'probe')
            save(BASE/'control01/binding.json',m)
    elif stage in ('export','export-recover','import'):
        assert json.loads((BASE/'build_result.json').read_text())['passed'];binary=PACKAGE/'probe'
        assert sha(binary)==json.loads((BASE/'build_result.json').read_text())['binarySHA256']
        case=BASE/'control01'
        if stage=='export-recover':
            assert not json.loads((BASE/'export01/report.json').read_text())['passed']
            shutil.copyfile(ROOT/'outputs/v213_validation/control02/phase-air.lua',case/'phase-air.lua')
        locks=common+[binary]+[p for p in case.rglob('*') if p.is_file()]
        remaining=60
        for name in ['export01','export02']:
            if (BASE/name/'report.json').exists():remaining-=json.loads((BASE/name/'report.json').read_text())['elapsedSeconds']
        mode='export' if stage=='export-recover' else stage
        record='export02' if stage=='export-recover' else stage+'01'
        r=invoke(['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary),mode],case,record,remaining,2048,64,locks)
    elif stage=='audit':
        case=BASE/'control01';locks=common+[ROOT/'tests/testNativeFrozenFile.py',ROOT/'scripts/nativeFrozenSnapshot.py']
        locks += [p for p in case.rglob('*') if p.is_file()]
        remaining=60-json.loads((BASE/'preparation01/report.json').read_text())['elapsedSeconds']
        r=invoke([sys.executable,str(ROOT/'tests/testNativeFrozenFile.py'),str(case)],BASE,'audit01',remaining,1024,64,locks)
    else:raise ValueError('stage')
    return r['passed']
if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
