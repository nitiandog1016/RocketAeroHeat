"""独立核对真实线性增量/冻结状态；没有证据不授权600步候选。"""
from pathlib import Path
import json,hashlib,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v117_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    folder=OUT/'diagnostic';target=OUT/'diagnostic_gate.json';assert not target.exists()
    receipt=json.loads((folder/'engine_process.json').read_text())
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut']
    assert receipt['frozenEnergyMode']=='diagnostic' and all(receipt['inputsStable'].values())
    for p,h in receipt['inputSHA256'].items():assert sha(Path(p))==h
    for p,h in receipt['outputSHA256'].items():assert sha(folder/p)==h
    initial=np.loadtxt(str(folder/'initial.csv'),delimiter=',',skiprows=1)
    control=np.loadtxt(str(ROOT/'outputs/v116_validation/study01/control/restart_flow.csv'),delimiter=',',skiprows=1)
    with (folder/'rah_frozen_directions.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)<=4096 and all(float(r['nativeAlpha'])==float(r['appliedAlpha'])==0 for r in rows)
    records=[]
    for node in [29214,31040,32866,34692]:
        q=initial[node,4:9];assert initial[node,0]==node and control[node,0]==node
        assert np.array_equal(q,control[node,4:9]),'目标不再符合原冻结证据'
        matched=[r for r in rows if int(r['iteration'])==0 and np.array_equal(q,np.array([float(r['q'+str(j)]) for j in range(5)]))]
        assert len(matched)==1,('不能将局部编号当原节点',node,len(matched))
        r=matched[0];delta=np.array([float(r['delta'+str(j)]) for j in range(5)])
        allowance=float(r['allowable']);assert allowance==.2
        # 独立表达式：先计算比动能再乘回密度，与C++逐项扣除有不同舍入。
        e=q[4]-np.dot(q[1:4],q[1:4])/(2*q[0])
        density_alpha=min(1,allowance*q[0]/abs(delta[0])) if delta[0]!=0 else 1
        trial=q+density_alpha*delta
        et=trial[4]-np.dot(trial[1:4],trial[1:4])/(2*trial[0])
        eligible=bool(np.isfinite(trial).all() and q[0]>0 and e>0 and trial[0]>0 and et>e and density_alpha>=1e-10)
        assert eligible==(r['eligible']=='1')
        assert abs(trial[0]/q[0]-1)<=.2+1e-14
        records.append(dict(globalOriginalNode=node,localRuntimeNode=int(r['localNode']),X_m=float(initial[node,1]),
            conservativeState=q.tolist(),actualLinearIncrement=delta.tolist(),internalBefore=float(e),
            internalTrial=float(et),densityBoundAlpha=float(density_alpha),eligible=eligible,
            frozenV114toV116=True))
    counts=np.atleast_1d(np.genfromtxt(str(folder/'rah_frozen_counts.csv'),names=True,delimiter=','))
    assert np.sum(counts['applied'])==0 and set(counts['iteration'])=={0,1}
    sums=[dict(iteration=i,checked=int(sum(counts['checked'][counts['iteration']==i])),
               frozen=int(sum(counts['frozen'][counts['iteration']==i])),eligible=int(sum(counts['eligible'][counts['iteration']==i]))) for i in [0,1]]
    result=dict(passed=all(r['eligible'] for r in records),eligibleTargetCount=sum(r['eligible'] for r in records),
        directionsSHA256=sha(folder/'rah_frozen_directions.csv'),receiptSHA256=sha(folder/'engine_process.json'),
        records=records,counts=sums,diagnosticUpdatesUnchanged=True,physicsQualified=False,
        sources={str(p):sha(p) for p in [Path(__file__),folder/'initial.csv',folder/'rah_frozen_directions.csv',folder/'rah_frozen_counts.csv']})
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
