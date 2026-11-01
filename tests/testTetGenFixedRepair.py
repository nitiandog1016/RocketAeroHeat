"""唯一带孔真值小例：准确导入、固定边界、体积26及真正改善SICN。"""
from pathlib import Path
import json,sys,time,traceback
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareTetGenRuntime import ROOT,sha
from tetgenFixedInterface import execute,tet_audit
from reviewSU2ExplicitQuality import quality
from testCGALLocalRemeshing import hollow,orient

def main():
    start=time.monotonic();folder=ROOT/'docs/validation/v142/hollow_test'
    folder.mkdir(parents=True,exist_ok=False)
    paths=[Path(__file__).resolve(),ROOT/'tests/testCGALLocalRemeshing.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',
        ROOT/'docs/validation/v142/plan.md']
    locked={str(p):sha(p) for p in paths};result=dict(passed=False,physicsQualified=False,newQualifiedArea=0)
    try:
        P,T=hollow();tet=T[13];point=P[tet].T@np.array([.333,.333,.333,.001]);nid=len(P)
        P=np.r_[P,point[None]];split=np.tile(tet,(4,1));split[np.arange(4),np.arange(4)]=nid
        T=np.r_[np.delete(T,13,axis=0),orient(P,split)]
        assert len(T)==159
        np.savez_compressed(str(folder/'fixture.npz'),vertices=P,tetrahedra=T)
        _,before=tet_audit(P,T);assert abs(before['volume']-26)<1e-12
        q0=quality(P[T]);Q,U,inter=execute(folder/'native',P,T,60,paths)
        _,after=tet_audit(Q,U);q1=quality(Q[U]);assert abs(after['volume']-26)<1e-11
        result.update(interface=inter,minimumSICNBefore=float(min(q0)),minimumSICNAfter=float(min(q1)),
            tetsBelowSICN01Before=int(sum(q0<.1)),tetsBelowSICN01After=int(sum(q1<.1)),
            expectedVolume=26,actualVolume=after['volume'])
        assert min(q1)>min(q0) and min(q1)>=.1,'quality improvement or SICN gate failed'
        result['passed']=True
    except Exception as e:
        result['error']=repr(e);(folder/'failure.log').write_text(traceback.format_exc())
    result.update(seconds=time.monotonic()-start,sourceSHA256=locked,sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()})
    assert all(result['sourceStableAfter'].values())
    (folder/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)

if __name__=='__main__':main()
