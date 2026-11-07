"""未通过小例数值门槛的实数v1不能准备大域或继续增加步数；不触发CFD。"""
from pathlib import Path
import argparse,copy,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import audit_run_scope
from prepareEilmerDevelopmentCase import prepare
from runPreparedEilmerCase import execute


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'outputs/v181_validation';m=json.loads((base/'small01/case_manifest.json').read_text())
    audit_run_scope(m);rejected=[]
    for label,cells,steps in [('full-domain',705490,8),('extended-small',4,400)]:
        bad=copy.deepcopy(m);bad['mesh']['cells']=cells;bad['condition']['maxNewtonSteps']=steps
        try:audit_run_scope(bad)
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    # 模块：实际运行入口在访问大网格/创建回执/启动二进制以前拒绝。
    fake=out/'fake-full';fake.mkdir();bad=copy.deepcopy(m);bad['mesh']['cells']=705490
    (fake/'case_manifest.json').write_text(json.dumps(bad))
    for stage in ['prepare','run']:
        try:execute(fake,stage,out/('must-not-'+stage),30,256)
        except ValueError:rejected.append(stage+'-entry-before-native')
        else:raise AssertionError('原生入口未拒绝')
        assert not (out/('must-not-'+stage)).exists()
    # 原复数发展入口不受此新数值方案范围限制。
    legacy=copy.deepcopy(bad);legacy['condition']['development']['profile']='rah-nk-development-v1';audit_run_scope(legacy)
    r=dict(passed=True,rejected=rejected,legacyScopeUnchanged=True,nativeRuns=0,
           sourceSHA256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                         ['scripts/readEilmerNative.py','scripts/runPreparedEilmerCase.py','scripts/prepareEilmerDevelopmentCase.py','tests/testEilmerRealAdmission.py']},
           elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);test(p.parse_args().out)
