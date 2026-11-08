"""定向保存原生剖面/边界/交换与调用链；只读取，不启动求解器。"""
from pathlib import Path
import json,shutil,time,subprocess,re
from runEilmerProcess import PREFIX,COMMIT,sha
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v212_validation'

def main():
    start=time.monotonic();out=BASE/'source01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    src=PREFIX/'source';prior=json.loads((ROOT/'outputs/v197_validation/build01/report.json').read_text())
    assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()==COMMIT
    selected=[];hits=[]
    # 模块一：原成功构建的准确边界源，避免拿其他版本同名接口作依据。
    for name,h in prior['inputSHA256'].items():
        p=Path(name)
        if not str(p).startswith(str(src)) or not p.is_file():continue
        if '/lmr/bc/' not in name and p.name not in ['flowstate.d','simcore_exchange.d','simcore_gasdynamic_step.d','ufluidblock.d']:continue
        assert sha(p)==h
        q=out/'source'/p.relative_to(src);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        selected.append(dict(path=str(p.relative_to(src)),sha256=h))
    # 模块二：Lua构造组合、官方profile/交换说明和算例入口，按内容命中保存小文件。
    tracked=subprocess.check_output(['git','-C',str(src),'ls-files'],text=True).splitlines()
    pattern=re.compile(r'InFlowBC_FromProfile|InOutFlowBC|ExchangeBC_MappedCell|StaticFlowProfile|frozen.?flow|Schwarz',re.I)
    for name in tracked:
        if not (name.startswith(('src/lmr/','doc/lmr/','examples/lmr/')) and name.endswith(('.lua','.adoc','.md','.d'))):continue
        p=src/name
        if not p.is_file() or p.stat().st_size>1024**2:continue
        text=p.read_text(errors='replace');matching=[dict(line=i,text=s) for i,s in enumerate(text.splitlines(),1) if pattern.search(s)]
        if not matching:continue
        hits.append(dict(path=name,lines=matching))
        if name.endswith('.lua') or name.startswith('doc/'):
            q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
            selected.append(dict(path=name,sha256=sha(p)))
    assert all(sha(src/x['path'])==x['sha256'] for x in selected)
    report=dict(passed=True,sourceCommit=COMMIT,files=selected,seconds=time.monotonic()-start,builds=0,CFDRuns=0)
    assert report['seconds']<120
    for name,obj in [('report.json',report),('hits.json',hits)]:
        (out/name).write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps(dict(passed=True,files=len(selected),seconds=report['seconds'])))
    print('\n'.join(x['path'] for x in hits))

if __name__=='__main__':main()
