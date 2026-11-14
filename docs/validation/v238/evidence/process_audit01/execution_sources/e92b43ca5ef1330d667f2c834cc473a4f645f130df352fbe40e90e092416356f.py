"""模块职责：独立解析Euler通量/一致态Jacobian、接触及能量变换审核；不实现第二套HLLC。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runHLLCCompatibilityV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def audit(case='interface01',output='audit01'):
    out=BASE/output;out.mkdir(exist_ok=False);air=PhaseAwareViscousAir();checks=[]
    with (BASE/case/'flux.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==24
    def check(name,a,b,tol,identity):
        a=np.asarray(a);b=np.asarray(b);error=float(np.max(abs(a-b)/np.maximum(1,abs(b))))
        checks.append(dict(name=name,case=identity,error=error,tolerance=tol,passed=error<=tol));return error
    def physical(r):
        p=float(r['pL']);t=float(r['TL']);v=np.array([float(r[x]) for x in ('uL','vL','wL')]);q=air.evaluate(p,t)
        rho=float(q['densityKgm3']);e=float(q['internalEnergyJkg']);kappa=air.R/float(q['cvJkgK']);chi=air.R*t-kappa*e
        u=np.r_[rho,rho*v,rho*(e+.5*sum(v*v))];pressure=np.r_[chi+.5*kappa*sum(v*v),-kappa*v,kappa]
        speed=np.array([-v[0]/rho,1/rho,0,0,0]);flux=np.r_[u[1],u[1:4]*v[0]+np.array([p,0,0]),(u[4]+p)*v[0]]
        jac=np.zeros((5,5));jac[0,1]=1
        for i in range(3):
            jac[i+1]=u[i+1]*speed; jac[i+1,i+1]+=v[0]
            if i==0:jac[i+1]+=pressure
        jac[4]=(u[4]+p)*speed+v[0]*pressure;jac[4,4]+=v[0]
        return flux,jac
    index={}
    for r in rows:
        name=r['case'];offset=float(r['offset']);reverse=int(r['reverse'])
        f=np.array([float(r['F%d'%j]) for j in range(5)]);jl=np.array([[float(r['JL%d%d'%(i,j)]) for j in range(5)] for i in range(5)]);jr=np.array([[float(r['JR%d%d'%(i,j)]) for j in range(5)] for i in range(5)])
        transform=np.eye(5);transform[4,0]=offset;inverse=np.eye(5);inverse[4,0]=-offset
        f=inverse@f;jl=inverse@jl@transform;jr=inverse@jr@transform
        assert np.isfinite(f).all() and np.isfinite(jl).all() and np.isfinite(jr).all()
        index[(name,offset,reverse)]=(f,jl,jr)
        if not reverse and r['exactLeft']=='1':
            exact,jac=physical(r);check('exactFlux',f,exact,2e-10,name+' offset '+str(offset))
            if name.startswith('uniform'):check('consistentStateTotalJacobian',jl+jr,jac,5e-9,name+' offset '+str(offset))
    for name in sorted({r['case'] for r in rows}):
        for offset in (0.,300000.):
            a=index[(name,offset,0)];b=index[(name,offset,1)]
            check('pairedFlux',a[0],-b[0],2e-10,name+' offset '+str(offset))
            check('pairedJacobian',a[1],-b[2],2e-8,name+' offset '+str(offset))
            check('pairedJacobian',a[2],-b[1],2e-8,name+' offset '+str(offset))
        for reverse in (0,1):
            a=index[(name,0.,reverse)];b=index[(name,300000.,reverse)]
            for i,kind in enumerate(('Flux','LeftJacobian','RightJacobian')):check('energyReference'+kind,b[i],a[i],2e-8,name+' reverse '+str(reverse))
    with (BASE/case/'rejected.csv').open() as f:rejects=list(csv.DictReader(f))
    assert len(rejects)==3
    passed=all(c['passed'] for c in checks);save(out/'checks.json',checks)
    report=dict(passed=passed,nativeCalls=24,cases=6,structuralRejections=3,independentEulerAnalyticOracle=True,
        failures=[c for c in checks if not c['passed']],maxError=max(c['error'] for c in checks),
        actualFactoryDispatchTested=False,wholeCombinationAccepted=False,flowUpdates=0,physicalQualification=False,
        sourceSHA256=sha(BASE/case/'flux.csv'))
    save(out/'report.json',report);print(json.dumps(report));return passed


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 interface audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':raise SystemExit(0 if audit(*sys.argv[2:]) else 2)
    run.configure();case=sys.argv[1] if len(sys.argv)>1 else 'interface01';output=sys.argv[2] if len(sys.argv)>2 else 'audit01'
    locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
           ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json']+list((BASE/case).iterdir())
    raise SystemExit(0 if run.guard.supervise(output,[sys.executable,str(Path(__file__).resolve()),'task',case,output],60,1024,locks,BASE,'audit') else 1)
