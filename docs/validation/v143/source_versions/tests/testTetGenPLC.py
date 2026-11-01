"""同一固定空腔表面的唯一PLC真值测试，不继承旧偏心内部点。"""
from pathlib import Path
import sys,json,time,traceback
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareTetGenRuntime import ROOT,sha
from cgalTetrahedralContract import faces_and_neighbors
from tetgenPLCInterface import execute
from reviewSU2ExplicitQuality import quality

def main():
    start=time.monotonic();folder=ROOT/'docs/validation/v143/hollow_test';folder.mkdir(parents=True,exist_ok=False)
    fixture=ROOT/'docs/validation/v142/hollow_test/fixture.npz'
    paths=[Path(__file__).resolve(),fixture,ROOT/'docs/validation/v143/plan.md',ROOT/'scripts/reviewSU2ExplicitQuality.py']
    locked={str(p):sha(p) for p in paths};r=dict(passed=False,physicsQualified=False,newQualifiedArea=0)
    try:
        old=np.load(str(fixture));P=old['vertices'];T=old['tetrahedra'];F,_=faces_and_neighbors(T);ids=np.unique(F)
        assert len(ids)==64 and len(F)==120;mapping=np.full(len(P),-1,np.int64);mapping[ids]=np.arange(len(ids))
        surface=P[ids];faces=mapping[F];hole=np.array([1.5,1.5,1.5]);fluid=P[T[0]].mean(axis=0)
        Q,U,inter=execute(folder/'native',surface,faces,hole,fluid,60,paths)
        q=quality(Q[U]);q0=quality(P[T]);r.update(interface=inter,minimumSICNBefore=float(min(q0)),minimumSICNAfter=float(min(q)),tetsBelowSICN01=int(sum(q<.1)))
        assert abs(inter['geometry']['volume']-26)<1e-11 and min(q)>=.1 and min(q)>min(q0),'quality gate failed'
        r['passed']=True
    except Exception as e:
        r['error']=repr(e);(folder/'failure.log').write_text(traceback.format_exc())
    r.update(seconds=time.monotonic()-start,sourceSHA256=locked,sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()})
    assert all(r['sourceStableAfter'].values());(folder/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2));raise SystemExit(0 if r['passed'] else 1)

if __name__=='__main__':main()
