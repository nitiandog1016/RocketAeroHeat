"""一次定向保存初始化/边界登记与原生NK入口；不编译、不执行流动。"""
from pathlib import Path
import json,subprocess,time,shutil
from runEilmerProcess import PREFIX,COMMIT,sha
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.monotonic();src=PREFIX/'source';out=ROOT/'outputs/v213_validation/source02';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(),out/'executed_entry.py')
    assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()==COMMIT
    prior=json.loads((ROOT/'outputs/v197_validation/build01/report.json').read_text())['inputSHA256']
    files=[];hits=[]
    for name,h in prior.items():
        p=Path(name)
        if not str(p).startswith(str(src/'src/lmr')) or not p.is_file() or p.suffix!='.d':continue
        text=p.read_text();words=('get_mapped_cell(', 'set_up_cell_mapping(', 'initMappedCellDataExchange', 'new GhostCellMappedCellCopy', 'initNewtonKrylovSimulation(')
        match=[dict(line=i,text=line) for i,line in enumerate(text.splitlines(),1) if any(w in line for w in words)]
        if not match and p.name not in ('ghost_cell_effect.d','bc.d','main.d','newtonkrylovsolver.d','globalconfig.d','globaldata.d'):continue
        assert sha(p)==h
        rel=p.relative_to(src);q=out/'source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        files.append(dict(path=str(rel),sha256=h));hits.append(dict(path=str(rel),lines=match))
    shutil.copyfile(Path(__file__).resolve(),out/'executed_entry.py')
    for name,data in [('report.json',dict(passed=True,files=files,seconds=time.monotonic()-start,sourceCommit=COMMIT)),('hits.json',hits)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    assert time.monotonic()-start<60
    print(json.dumps(dict(files=len(files),seconds=time.monotonic()-start)));print(json.dumps(hits,indent=2))

if __name__=='__main__':raise ValueError('v213定向源提取已完成；复用归档，禁止自动重做。')
