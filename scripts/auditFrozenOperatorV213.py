"""模块职责：从逐面/逐单元输出独立复算原生对照误差及控制体收支，不重跑求解。"""
from pathlib import Path
import json,time,csv
import numpy as np
from readEilmerNative import read_grid,compare_source,sha
from auditFrozenMappedGeometryV213 import volume
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'

def scaled(a,b):return np.abs(a-b)/np.maximum(1.,np.maximum(np.abs(a),np.abs(b)))

def main():
    start=time.monotonic();case=BASE/'control03';out=BASE/'operator_audit01';out.mkdir(exist_ok=False)
    grid=read_grid(case/'lmrsim/grid/grid-0000');compare_source(grid,case/'mesh.su2')
    mapping=np.load(str(case/'mapping.npz'));native=np.load(str(case/'native_connectivity.npz'))
    assert np.array_equal(grid['points'],mapping['points']) and np.array_equal(grid['cells'][:,:6],native['cells'])
    assert np.all(volume(grid['points'][grid['cells'][:,:6]])>0)
    a=np.genfromtxt(str(case/'samples.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    assert len(a)==3*15918 and set(a['regime'])=={0,1,2}
    for f in ['resident','frozen','residentDerivative','frozenDerivative','plus','minus']:assert np.isfinite(a[f]).all()
    geom=np.genfromtxt(str(case/'geometry.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    cells=geom[geom['kind']=='cell'];assert np.array_equal(cells['id'],np.arange(416))
    assert np.all(cells['volumeOrArea']>0) and np.all(cells['wallDistance']>0)
    grad=np.genfromtxt(str(case/'gradient.csv'),delimiter=',',names=True)
    expected=np.array([3,-5,7,20,-10,3,-5,8,12,9,4,-7])
    actual=np.stack([grad[x] for x in ['Tx','Ty','Tz','ux','uy','uz','vx','vy','vz','wx','wy','wz']],axis=1)
    gradient=float(np.max(np.abs(actual-expected)/np.maximum(1,np.abs(expected))))
    normal=np.stack([grad[x] for x in ['nx','ny','nz']],axis=1)
    nref=normal.dot(expected[:3]);nact=np.sum(normal*actual[:,:3],axis=1)
    normal_error=float(np.max(np.abs(nref-nact)/np.maximum(1,np.abs(nref))))
    assert gradient<=1e-7 and normal_error<=1e-7
    reports=[];worst=[]
    for regime in range(3):
        r=a[a['regime']==regime];assert np.array_equal(r['index'],np.arange(15918))
        fd=(r['plus']-r['minus'])/2e-4;err=scaled(r['frozenDerivative'],fd)
        base=scaled(r['resident'],r['frozen']);deriv=scaled(r['residentDerivative'],r['frozenDerivative'])
        groups={}
        for kind in np.unique(r['kind']):
            choose=r['kind']==kind;indices=np.flatnonzero(choose);i=indices[np.argmax(err[choose])]
            ids=r['id'][choose];comp=r['component'][choose]
            width=dict(residual=6,temperatureGradient=3,velocityGradient=9,saGradient=3,integratedFlux=6)[kind]
            count=1197 if kind=='integratedFlux' else 416
            assert set(zip(ids,comp))=={(j,k) for j in range(count) for k in range(width)}
            groups[kind]=dict(count=int(choose.sum()),baseMax=float(base[choose].max()),derivativeMax=float(deriv[choose].max()),
                finiteDifferenceMax=float(err[choose].max()),finiteDifferenceFailed=int((err[choose]>1e-5).sum()))
            row=dict(regime=regime,kind=kind,id=int(r[i]['id']),component=int(r[i]['component']),
                originalCellID=None if kind=='integratedFlux' else int(mapping['motherCellRows'][r[i]['id']]),
                frozen=float(r[i]['frozen']),derivative=float(r[i]['frozenDerivative']),
                plus=float(r[i]['plus']),minus=float(r[i]['minus']),finiteDifference=float(fd[i]),scaledError=float(err[i]))
            worst.append(row)
        # 前五守恒式无体源：单元体积分残差＋所有边界的有向通量应抵消。
        residual=r[r['kind']=='residual']['frozen'].reshape(416,6)
        flux=r[r['kind']=='integratedFlux']['frozen'].reshape(1197,6)
        iv=(residual*cells['volumeOrArea'][:,None]).sum(axis=0)
        boundary=np.zeros(6);abs_boundary=np.zeros(6)
        for b in grid['boundaries'].values():
            q=flux[b['ids']]*b['signs'][:,None];boundary+=q.sum(axis=0);abs_boundary+=np.abs(q).sum(axis=0)
        balance=(iv+boundary)[:5]/np.maximum(1,np.maximum(np.abs(iv[:5]),abs_boundary[:5]))
        assert np.max(np.abs(balance))<=1e-11
        reports.append(dict(regime=regime,groups=groups,scaledConservationBalance=balance.tolist(),
            baseMax=float(base.max()),derivativeMax=float(deriv.max()),finiteDifferenceMax=float(err.max()),
            passed=bool(base.max()<=1e-11 and deriv.max()<=1e-11 and err.max()<=1e-5)))
    original=[json.loads(s) for s in (BASE/'control02/operator_results.jsonl').read_text().splitlines()]
    current=[json.loads(s) for s in (case/'operator_results.jsonl').read_text().splitlines()]
    assert original==current[:2]  # 原失败不因新增观测改变。
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/k.replace('\\','/'))==h for k,h in protected.items())
    result=dict(auditPassed=True,nativeAcceptance=all(r['passed'] for r in reports),
        gradientError=gradient,normalGradientError=normal_error,regimes=reports,worstComponents=worst,
        sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in [case/'samples.csv',case/'geometry.csv',case/'gradient.csv',case/'operator_results.jsonl',Path(__file__).resolve()]},
        protectedFilesUnchanged=29,originalFailureUnchanged=True,CFDAdvanceRuns=0,newQualifiedAreaM2=0,
        limitation='frozen resident donors; not coupled Jacobian or converged regional solution',seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(auditPassed=True,nativeAcceptance=result['nativeAcceptance'],gradient=gradient,
        cases=[dict(regime=r['regime'],fd=r['finiteDifferenceMax'],balance=max(map(abs,r['scaledConservationBalance']))) for r in reports],seconds=result['seconds'])))

if __name__=='__main__':main()
