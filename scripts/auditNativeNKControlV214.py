"""模块职责：独立核对守恒变量切向、实际Jv分量、同阶实差分及守恒账。"""
from pathlib import Path
import json,time
import numpy as np
from readEilmerNative import read_grid,compare_source,sha
from phaseAwareDryAir import PhaseAwareDryAir
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v214_validation'

def scaled(a,b):return np.abs(a-b)/np.maximum(1,np.maximum(np.abs(a),np.abs(b)))

def main():
    start=time.monotonic();case=BASE/'control01';out=BASE/'audit01';out.mkdir(exist_ok=False)
    z=np.load(str(case/'mapping.npz'));grid=read_grid(case/'lmrsim/grid/grid-0000');compare_source(grid,case/'mesh.su2')
    old=ROOT/'outputs/v213_validation/control02'
    prep=json.loads((BASE/'prepared01/report.json').read_text())
    for name,h in prep['caseInputSHA256'].items():assert sha(case/name)==h==sha(old/name)
    data=np.genfromtxt(str(case/'directions.csv'),delimiter=',',names=True)
    assert len(data)==1248
    gas=PhaseAwareDryAir();thermo=gas.evaluate(data['p'],data['T'])
    rho=thermo['densityKgm3'];cv=thermo['cvJkgK'];u=thermo['internalEnergyJkg']
    assert max(scaled(rho,data['rho']).max(),scaled(cv,data['Cv']).max(),scaled(u,data['u']).max())<1e-11
    drho=rho*(1/data['p']-.2/data['T']);vel=np.stack([data[k] for k in ['vx','vy','vz']],axis=1);dv=np.array([.3,-.4,.2])
    momentum=vel*drho[:,None]+rho[:,None]*dv
    de=(u+np.sum(vel*vel,axis=1)/2)*drho+rho*(cv*.2+vel.dot(dv))
    expected=np.column_stack([drho,momentum,de,drho*data['nuhat']+rho*1e-7])
    direction=np.column_stack([data[k] for k in ['drho','dmx','dmy','dmz','dE','dsa']])
    direction_error=float((np.abs(direction-expected)/np.maximum(1e-20,np.abs(expected))).max());assert direction_error<1e-11
    a=np.genfromtxt(str(case/'nk_samples.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    assert len(a)==47754
    fields=['resident','frozen','residentDerivative','frozenDerivative','residentPlus','residentMinus','frozenPlus','frozenMinus']
    assert all(np.isfinite(a[f]).all() for f in fields)
    g=np.genfromtxt(str(ROOT/'outputs/v213_validation/control03/geometry.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    volume=g[g['kind']=='cell']['volumeOrArea'];results=[];witnesses=[]
    native=[json.loads(s) for s in (case/'nk_results.jsonl').read_text().splitlines()];assert len(native)==3
    for regime in range(3):
        r=a[a['regime']==regime];assert np.array_equal(r['index'],np.arange(15918))
        eq=max(float(scaled(r[x],r[y]).max()) for x,y in [('resident','frozen'),('residentDerivative','frozenDerivative'),('residentPlus','frozenPlus'),('residentMinus','frozenMinus')])
        fdr=(r['residentPlus']-r['residentMinus'])/2e-4;fdf=(r['frozenPlus']-r['frozenMinus'])/2e-4
        er=scaled(fdr,r['residentDerivative']);ef=scaled(fdf,r['frozenDerivative'])
        groups={}
        for kind in np.unique(r['kind']):
            mask=r['kind']==kind;ids=np.flatnonzero(mask);k=int(ids[np.argmax(ef[mask])])
            width={'residual':6,'integratedFlux':6,'temperatureGradient':3,'velocityGradient':9,'saGradient':3}[kind]
            count=1197 if kind=='integratedFlux' else 416
            assert set(zip(r['id'][mask],r['component'][mask]))=={(i,j) for i in range(count) for j in range(width)}
            groups[kind]=dict(count=int(mask.sum()),residentFDMax=float(er[mask].max()),frozenFDMax=float(ef[mask].max()),failed=int((ef[mask]>1e-5).sum()))
            witnesses.append(dict(regime=regime,kind=kind,id=int(r[k]['id']),component=int(r[k]['component']),
                motherCellID=None if kind=='integratedFlux' else int(z['motherCellRows'][r[k]['id']]),
                base=float(r[k]['frozen']),derivative=float(r[k]['frozenDerivative']),
                plus=float(r[k]['frozenPlus']),minus=float(r[k]['frozenMinus']),finiteDifference=float(fdf[k]),error=float(ef[k])))
        rr=r[r['kind']=='residual']['frozen'].reshape(416,6);ff=r[r['kind']=='integratedFlux']['frozen'].reshape(1197,6)
        integrated=(rr*volume[:,None]).sum(axis=0);boundary=np.zeros(6);scale=np.zeros(6)
        for b in grid['boundaries'].values():
            flux=ff[b['ids']]*b['signs'][:,None];boundary+=flux.sum(axis=0);scale+=np.abs(flux).sum(axis=0)
        balance=(integrated+boundary)[:5]/np.maximum(1,np.maximum(np.abs(integrated[:5]),scale[:5]))
        assert np.max(np.abs(balance))<1e-11
        nr=native[regime];assert nr['regime']==regime and nr['components']==15918 and nr['equivalenceError']==eq
        assert nr['zedError']<=1e-11 and nr['residentFDError']==float(er.max()) and nr['frozenFDError']==float(ef.max())
        passed=eq<=1e-11 and float(er.max())<=1e-5 and float(ef.max())<=1e-5;assert passed==nr['passed']
        results.append(dict(regime=regime,equivalenceError=eq,zedError=nr['zedError'],groups=groups,conservationBalance=balance.tolist(),passed=passed))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/k.replace('\\','/'))==h for k,h in protected.items())
    report=dict(auditPassed=True,nativeAcceptance=all(r['passed'] for r in results),directionRelativeError=direction_error,
        regimes=results,witnesses=witnesses,protectedFilesUnchanged=29,originalInputsUnchanged=True,
        CFDAdvanceRuns=0,newQualifiedAreaM2=0,fullKrylovOrCoupledJacobianVerified=False,
        sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in [case/'nk_samples.csv',case/'directions.csv',case/'nk_results.jsonl',Path(__file__).resolve()]},
        seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(auditPassed=True,nativeAcceptance=report['nativeAcceptance'],directionRelativeError=direction_error,
        results=native,seconds=report['seconds'])))

if __name__=='__main__':main()
