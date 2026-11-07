"""实际原生对照：修正四面体时间尺度，并保持六面体数值不变。"""
from pathlib import Path
import json
import re
import shutil
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,sha
from readEilmerNative import read_field
from observeEilmerCellABI import GDB
from hybridRocketEuler import SavedFieldGuard

BASE=ROOT/'outputs/v196_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-geometric-length-v1'

def main():
    """模块一：准确包/来源、独立几何真值和原生单步，禁止静态替身。"""
    start=time.monotonic();out=BASE/'small_validation01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    package=json.loads((PACKAGE/'runtime.json').read_text())
    assert package['profile']=='rah-geometric-length-v1' and all(sha(p)==h for p,h in package['filesSHA256'].items())
    rows=[]
    for kind in ['tet','hex','hex_baseline']:
        source=BASE/('tet02' if kind=='tet' else 'hex01')
        case=BASE/(kind+'_fixed01' if kind!='hex_baseline' else 'hex_baseline01')
        shutil.copytree(source,case)
        # 输入复制不携带旧观测，以免旧值冒充本次证据。
        old=case/'time-observation.json';old.rename(case/'prior-time-observation.json')
        script=GDB.replace('quit\n',
            "python\ngdb.events.exited.connect(lambda e: open('inferior-exit.json','w').write(json.dumps(dict(exitCode=e.exit_code))))\nend\ncontinue\nquit\n")
        (case/'observe.gdb').write_text(script)
        binary=PACKAGE/'lmr-run' if kind!='hex_baseline' else ROOT/'tools/runtime/eilmer-explicit-workspace-v1/lmr-run'
        inputs=[p for p in case.rglob('*') if p.is_file()]
        p=run(['gdb','--batch','-x','observe.gdb',str(binary)],case,out/(kind+'_native'),
            min(20,120-(time.monotonic()-start)),32,[Path(__file__),PACKAGE/'runtime.json',binary,*inputs],4096,4096)
        assert p['passed'] and json.loads((case/'inferior-exit.json').read_text())['exitCode']==0
        obs=json.loads((case/'time-observation.json').read_text())
        assert np.isfinite(obs['lengths']).all() and min(obs['lengths'])>0
        expected=(max(abs(np.dot(obs['velocity'],obs['normals'][i])) for i in obs['selectedFaces'])+obs['soundSpeed'])/min(obs['lengths'])
        assert abs(expected/obs['signal']-1)<2e-14
        if kind=='tet':
            z=np.load(ROOT/'outputs/v194_validation/netgen_recovered02/mesh.npz');v=z['vertices'][z['tetrahedra'][0]]
            vol=abs(np.linalg.det((v[1:]-v[0]).T))/6
            length=min(vol**(1/3),*(2*np.linalg.norm(v[[j for j in range(4) if j!=i]].mean(axis=0)-v.mean(axis=0)) for i in range(4)))
            assert np.max(np.abs(np.array(obs['lengths'])/length-1))<2e-13
        guard=SavedFieldGuard(case,1,out/(kind+'_guard.json'),maximum_dt=.001)
        assert guard() is None
        f0=read_field(case,0,1);f1=read_field(case,1,1)
        change=max(float(np.max(abs(f1[k]-v)/np.maximum(abs(v),1))) for k,v in f0.items() if k!='dt_local')
        assert change<2e-12 and np.all(f1['dt_local']>0)
        times=(case/'lmrsim/snapshots/snapshot-times-metadata').read_text()
        entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)\s*\n\s+dt:\s*(\S+)",times,re.M)
        assert len(entries)==2 and int(entries[-1][1])==1
        dt=float(entries[-1][3]);assert abs(dt/(.3/expected)-1)<2e-14
        rows.append(dict(kind=kind,signal=obs['signal'],lengths=obs['lengths'],dt=dt,
                         uniformStateMaxScaledChange=change,process=p))
    hexa=BASE/'hex_fixed01';old=BASE/'hex_baseline01'
    equal={f:sha(hexa/f)==sha(old/f) for f in ['lmrsim/snapshots/0000/fluid-0000','lmrsim/snapshots/0001/fluid-0000']}
    assert all(equal.values())
    report=dict(passed=True,rows=rows,hexFieldsExactlyEqual=equal,elapsedSeconds=time.monotonic()-start,
        physicalQualification=False,scope='one uniform ambient cell, first signal and one update; no wall heat or rocket qualification',
        sourceSHA256={str(p):sha(p) for p in [Path(__file__),PACKAGE/'runtime.json',BASE/'observation03/report.json']})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
    return True

if __name__=='__main__':raise SystemExit(0 if main() else 1)
