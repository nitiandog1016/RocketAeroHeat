"""模块职责：整体Euler控制的独立特征/滑移/能量/重启审核，不推进流场。"""
from pathlib import Path
import csv,copy,json,math,sys
import numpy as np
from runPhaseEulerCompatibilityV229 import ROOT,BASE,BUILD,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseEulerContract import bind,read_restart

def rows(path):
    with path.open(encoding='utf8') as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]

def main(stage):
    origin=BASE/'control01';out=BASE/('audit01' if stage=='audit' else 'restart_audit01');out.mkdir(exist_ok=False)
    air=PhaseAwareViscousAir();errors={}
    def check(name,a,b,limit=2e-9):
        a,b=np.asarray(a),np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all()
        e=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))));errors[name]=max(errors.get(name,0),e)
        assert e<=limit,(name,e)
    def state(p,t,v):
        q=air.evaluate(p,t);rho=float(q['densityKgm3']);e=float(q['internalEnergyJkg'])+300000
        return np.r_[rho,rho*v,rho*(e+.5*(v@v))],q
    metadata=bind(origin/'restart_flow.csv',origin/'mesh.su2',origin/'case.cfg',BUILD/'SU2_CFD.exe',BUILD/'compatibility.exe','manufactured-control')
    frame=read_restart(origin/'restart_flow.csv',metadata,metadata,origin/'mesh.su2')
    def mapped(path):
        data=rows(path);mapping={tuple(r[k] for k in ('x','y','z')):np.array([r[k] for k in ('rho','rhou','rhov','rhow','rhoE')]) for r in data}
        assert len(mapping)==len(data)==24;return mapping
    saved=mapped(origin/'saved.csv')
    for xyz,U in zip(frame['xyz'],frame['nativeConservative']):assert U.tobytes()==saved[tuple(xyz)].tobytes()
    if stage=='restart_audit':
        assert read(BASE/'process_restart/report.json')['passed'];second=mapped(BASE/'restart01/restarted.csv');assert set(saved)==set(second)
        for xyz in saved:assert saved[xyz].tobytes()==second[xyz].tobytes(),xyz
        cfg=(origin/'case.cfg').read_text(encoding='utf8')
        assert (BASE/'restart01/case.cfg').read_text(encoding='utf8')==cfg.replace('RESTART_SOL= NO','RESTART_SOL= YES')
        assert read(origin/'restart_flow.rah.json')==metadata
        restarted=read_restart(BASE/'restart01/restart_flow.csv',metadata,metadata,BASE/'restart01/mesh.su2')
        assert restarted['wallHeatFlux'] is None and restarted['solidTemperature'] is None
        raw=(origin/'restart_flow.csv').read_text(encoding='utf8').splitlines(True)
        first=raw[1].split(',');first[1]=str(float(first[1])+.001)
        variants={'missing_node':raw[:-1],'duplicate_node':raw[:-1]+[raw[1]],'wrong_coordinate':raw[:1]+[','.join(first)]+raw[2:]}
        refusals=[]
        for name,lines in variants.items():
            p=out/(name+'.csv');p.write_text(''.join(lines),encoding='utf8');changed=dict(metadata,payloadSHA256=sha(p))
            try:read_restart(p,changed,changed,origin/'mesh.su2')
            except ValueError as exc:refusals.append(dict(case=name,reason=str(exc)))
            else:raise AssertionError(name)
        save(out/'report.json',dict(passed=True,actualRestartNodesBitwise=24,conservativeValuesBitwise=120,
          wholeNativeInterfacePassed=True,manufacturedControlOnly=True,CFDIterations=0,newHeatAreaM2=0,newQualifiedAreaM2=0,
          actualEulerEngineBuilt=True,fullJacobianVerified=False,physicalSteadyBalanceVerified=False,
          rejectedFiles=refusals,control=read(BASE/'audit01/report.json')))
        print(json.dumps(dict(passed=True,restartValuesBitwise=120,fileRefusals=len(refusals))));return
    assert read(BASE/'process_control/report.json')['passed'] and read(origin/'dispatch.json')['passed']
    q=air.evaluate(1200,350.123);a=float(q['soundSpeedMps']);v=a*np.array([2.,.4,.3]);Ui,q=state(1200,350.123,v)
    qe=air.evaluate(1197,350);ve=float(qe['soundSpeedMps'])*np.array([2.001,.401,.299]);Ue,_=state(1197,350,ve)
    rho=Ui[0];kappa=air.R/float(q['cvJkgK']);chi=air.R*350.123-kappa*(float(q['internalEnergyJkg'])+300000)
    classes=set();energy=[]
    for r in rows(origin/'boundary.csv'):
        normal=np.array([r[k] for k in ('nx','ny','nz')]);n=normal/np.linalg.norm(normal);vn=v@n
        d=Ue-Ui;dv=(d[1:4]-v*d[0])/rho;dvn=dv@n
        dp=chi*d[0]+kappa*(d[4]-v@d[1:4]+.5*(v@v)*d[0]);dr=0.;vel=np.zeros(3);pres=0.
        if vn<0:dr+=d[0]-dp/a**2;vel+=dv-dvn*n
        for sign in (-1,1):
            if vn+sign*a<0:
                wave=.5*(dp+sign*rho*a*dvn);pres+=wave;dr+=wave/a**2;vel+=sign*wave/(rho*a)*n
        dm=rho*vel+v*dr;de=(pres-chi*dr)/kappa+v@dm-.5*(v@v)*dr;Ub=Ui+np.r_[dr,dm,de]
        observed=np.r_[r['rho'],r['rho']*np.array([r[k] for k in ('u','v','w')]),r['rho']*r['H']-r['p']]
        check('characteristicState',observed,Ub)
        qb=air.evaluate(r['p'],r['T']);check('boundaryEOS',r['rho'],float(qb['densityKgm3']))
        ub=np.array([r[k] for k in ('u','v','w')]);check('boundaryH',r['H'],float(qb['enthalpyJkg'])+300000+.5*(ub@ub))
        flux=np.r_[r['rho']*(ub@normal),r['rho']*(ub@normal)*ub+r['p']*normal,r['rho']*r['H']*(ub@normal)]
        native=np.array([r['R%d'%j] for j in range(5)]);check('boundaryFlux',native,flux)
        physical=native.copy();physical[4]-=300000*physical[0]
        energy.append(dict(node=int(r['node']),native=native.tolist(),physical=physical.tolist()))
        classes.add('all-in' if vn<-a else 'all-out' if vn>a else 'mixed-in' if vn<0 else 'mixed-out')
    assert classes=={'all-in','all-out','mixed-in','mixed-out'}
    initial=rows(origin/'initial.csv');coords={int(r['point']):np.array([r[k] for k in ('x','y','z')]) for r in initial}
    residual=np.array([.17,.31,-.27,.19,17.]);wallrows=rows(origin/'wall.csv');assert len(wallrows)==6
    for r in wallrows:
        n=np.array([r[k] for k in ('nx','ny','nz')]);check('unitNormal',n@n,1.)
        direction=coords[int(r['node'])];direction/=np.linalg.norm(direction)
        check('geometryNormal',abs(n@direction),1.)
        P=np.eye(5);P[1:4,1:4]-=np.outer(n,n);expected=P@Ui
        observed=np.array([r['U%d'%j] for j in range(5)]);check('slipMomentumAndEnergy',observed,expected)
        check('slipResidual',[r['R%d'%j] for j in range(5)],P@residual)
        check('slipDiagonalJacobian',[r['J%d'%j] for j in range(25)],P.ravel())
        vp=observed[1:4]/observed[0];qb=air.evaluate(r['p'],r['T']);check('projectedEOS',observed[0],float(qb['densityKgm3']))
        check('projectedEnergy',observed[4]/observed[0],float(qb['internalEnergyJkg'])+300000+.5*(vp@vp))
    q0=air.evaluate(1197.0279007383137,226.50908361133003)
    v0=5*float(q0['soundSpeedMps'])*np.array([-math.cos(math.radians(20)),0,math.sin(math.radians(20))]);U0,_=state(1197.0279007383137,226.50908361133003,v0)
    for r in initial:check('initialConservative',[r[k] for k in ('rho','rhou','rhov','rhow','rhoE')],U0)
    with (origin/'refusals.csv').open(encoding='utf8') as f:native_refusals=list(csv.DictReader(f))
    assert len(native_refusals)==2 and all(r['rejected']=='1' for r in native_refusals)
    count=0
    for key in metadata:
        bad=copy.deepcopy(metadata);bad[key]='wrong'
        try:read_restart(origin/'restart_flow.csv',bad,metadata,origin/'mesh.su2')
        except ValueError:count+=1
        else:raise AssertionError(key)
    for bad in (None,{},dict(metadata,extra=1)):
        try:read_restart(origin/'restart_flow.csv',bad,metadata,origin/'mesh.su2')
        except ValueError:count+=1
        else:raise AssertionError('schema')
    save(origin/'restart_flow.rah.json',metadata);save(out/'boundary_energy.json',energy)
    save(out/'report.json',dict(passed=True,errors=errors,characteristicClasses=sorted(classes),boundaryRows=6,wallRows=6,
      nativeStateRefusals=native_refusals,identityRefusals=count,initialNodes=24,savedValuesBitwise=120,
      actualRestartVerified=False,fullFlowJacobianVerified=False,CFDIterations=0,newQualifiedAreaM2=0))
    print(json.dumps(dict(passed=True,errors=errors,identityRefusals=count)))

if __name__=='__main__':main(sys.argv[1])
