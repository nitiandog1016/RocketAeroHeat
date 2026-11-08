"""模块职责：独立特征投影/TMR源项/热壁与能量保存审核，不重写流场内核。"""
from pathlib import Path
import csv,json,math,copy,sys
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2EnergyReference import flow_packet,convert_flow,SOLVER_REFERENCE,PHYSICAL_REFERENCE
from su2PhaseAirContract import bind,read_restart
from runSU2CompatibilityV221 import ROOT,BASE,BUILD,sha,read,save

def rows(name):
    with (BASE/'control01'/name).open(encoding='utf8') as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def main():
    out=BASE/'audit01';out.mkdir(exist_ok=False);air=PhaseAwareViscousAir();err={}
    def check(name,a,b,limit=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all()
        e=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))));err[name]=max(err.get(name,0.),e)
        assert e<=limit,(name,e,a,b)
    def state(p,t,u):
        q=air.evaluate(p,t);rho=float(q['densityKgm3']);e=float(q['internalEnergyJkg'])+300000
        return np.r_[rho,rho*u,rho*(e+.5*np.dot(u,u))],q
    qi=air.evaluate(1200,350.123);a=float(qi['soundSpeedMps']);u=a*np.array([2.,.4,.3]);Ui,qi=state(1200,350.123,u)
    qe=air.evaluate(1197,350);ue=float(qe['soundSpeedMps'])*np.array([2.001,.401,.299]);Ue,qe=state(1197,350,ue)
    rho=Ui[0];cv=float(qi['cvJkgK']);kappa=air.R/cv;e=float(qi['internalEnergyJkg'])+300000
    chi=air.R*350.123-kappa*e;mu=float(air.interpolated_transport(1200,350.123)['dynamicViscosityPas']);nu=mu/rho
    classes=set();records=[]
    for row in rows('boundary.csv'):
        normal=np.array([row[k] for k in ('nx','ny','nz')]);area=np.linalg.norm(normal);n=normal/area;un=u@n
        d=Ue-Ui;du=(d[1:4]-u*d[0])/rho;dun=du@n
        dp=chi*d[0]+kappa*(d[4]-u@d[1:4]+.5*(u@u)*d[0]);dr=0.;dv=np.zeros(3);pres=0.
        if un<0:dr+=d[0]-dp/a**2;dv+=du-dun*n
        for sign in (-1,1):
            if un+sign*a<0:
                wave=.5*(dp+sign*rho*a*dun);pres+=wave;dr+=wave/a**2;dv+=sign*wave/(rho*a)*n
        dm=rho*dv+u*dr;dE=(pres-chi*dr)/kappa+u@dm-.5*(u@u)*dr;Ub=Ui+np.r_[dr,dm,dE]
        ub=np.array([row[k] for k in ('u','v','w')]);observed=np.r_[row['rho'],row['rho']*ub,row['rho']*row['H']-row['p']]
        check('characteristicConservative',observed,Ub)
        qb=air.evaluate(row['p'],row['T']);check('boundaryEOS',row['rho'],float(qb['densityKgm3']))
        check('boundaryShiftedH',row['H'],float(qb['enthalpyJkg'])+300000+.5*(ub@ub))
        qbtr=air.interpolated_transport(row['p'],row['T']);check('boundaryCp',row['cp'],float(qb['cpJkgK']))
        check('boundaryTransportRelative',[row['mu']/float(qbtr['dynamicViscosityPas']),row['k']/float(qbtr['conductivityWmK'])],[1.,1.])
        qn=.5*(u+ub)@normal;outside=3*float(qbtr['dynamicViscosityPas'])/row['rho'];inside=row['sign']*3*nu
        check('SAUpwind',row['saResidual'],max(qn,0)*inside+min(qn,0)*outside)
        classes.add('all-in' if un<-a else 'all-out' if un>a else 'mixed-in' if un<0 else 'mixed-out')
        native=np.array([row['flowR%d'%j] for j in range(5)])
        physical=convert_flow(flow_packet(native,SOLVER_REFERENCE,'residual'),PHYSICAL_REFERENCE)['values']
        check('boundaryMassFlux',native[0],row['rho']*(ub@normal))
        records.append(dict(sign=row['sign'],localNode=int(row['node']),native=native.tolist(),physical=physical.tolist()))
    assert classes=={'all-in','all-out','mixed-in','mixed-out'}
    for row in rows('wall.csv'):
        keff=row['k']+row['cp']*row['mut']/.9;factor=keff*row['area']/row['distance']
        check('wallHeatResidual',row['RE'],factor*(row['Tnear']-300))
        check('wallJ40',row['J40'],factor*(row['u2']-row['E'])/(row['rho']*row['cv']),3e-6)
        check('wallJ44',row['J44'],factor/(row['rho']*row['cv']),3e-6)
        check('wallZero',[row[k] for k in ('oldMx','oldMy','oldMz','oldSA','saResidual')],np.zeros(5))
        check('wallSADiagonal',row['saDiag'],1.)
    # TMR SA-neg，正支采用Note 1 option c的S-hat；本例零梯度，二次扩散项为0。
    def source_term(nt,r):
        om=r['omega'];dist=r['distance'];vol=r['volume'];nu0=r['mu']/r['rho'];cw1=.1355/.41**2+1.622/(2/3)
        if nt<=0:return (.1355*(1-1.2)*om*nt+cw1*(nt/dist)**2)*vol
        x=nt/nu0;fv1=x**3/(x**3+7.1**3);fv2=1-x/(1+x*fv1);sbar=nt*fv2/(.41*dist)**2
        st=om+sbar if sbar>=-.7*om else om+om*(.7**2*om+.9*sbar)/((.9-1.4)*om-sbar)
        st=max(st,1e-10);rr=min(nt/(st*(.41*dist)**2),10);g=rr+.3*(rr**6-rr);fw=g*(65/(g**6+64))**(1/6)
        ft2=1.2*math.exp(-.5*x*x)
        return (.1355*(1-ft2)*st*nt-(cw1*fw-.1355/.41**2*ft2)*(nt/dist)**2)*vol
    for row in rows('sa_source.csv'):
        nt=row['nutilde'];x=nt/(row['mu']/row['rho']);mut=0 if nt<0 else row['rho']*nt*x**3/(x**3+7.1**3)
        check('SAEddyViscosityRelative',row['mut']/mut if mut else row['mut'],1 if mut else 0)
        check('SASource',row['residual'],source_term(nt,row));step=abs(nt)*1e-6
        check('SASourceJacobian',row['jacobian'],(source_term(nt+step,row)-source_term(nt-step,row))/(2*step),3e-6)
    initial=rows('initial.csv');q0=air.evaluate(1197.0279007383137,226.50908361133003)
    v0=5*float(q0['soundSpeedMps'])*np.array([-math.cos(math.radians(20)),0,math.sin(math.radians(20))])
    U0,_=state(1197.0279007383137,226.50908361133003,v0)
    for row in initial:
        check('actualInitialization',[row[k] for k in ('rho','rhou','rhov','rhow','rhoE')],U0)
    data=BASE/'control01';metadata=bind(data/'restart_flow.csv',data/'mesh.su2',data/'case.cfg',
        ROOT/'outputs/v220_validation/build02/SU2_CFD.exe',BUILD/'compatibility.exe','manufactured-control')
    frame=read_restart(data/'restart_flow.csv',metadata,metadata,data/'mesh.su2')
    saved=rows('saved.csv');mapping={tuple(row[k] for k in ('x','y','z')):row for row in saved};assert len(mapping)==24
    for xyz,U,nt in zip(frame['xyz'],frame['nativeConservative'],frame['nuTilde']):
        row=mapping[tuple(xyz)];expected=np.array([row[k] for k in ('rho','rhou','rhov','rhow','rhoE')])
        assert U.tobytes()==expected.tobytes() and np.float64(nt).tobytes()==np.float64(row['nutilde']).tobytes()
    rejection=0
    for key in ('schema','model','meshSHA256','engineSHA256','producerSHA256','payloadSHA256','sourceThermoSHA256','role','turbulence','saVariable'):
        bad=copy.deepcopy(metadata);bad[key]='wrong'
        try:read_restart(data/'restart_flow.csv',bad,metadata,data/'mesh.su2')
        except ValueError:rejection+=1
        else:raise AssertionError(key)
    for bad in (None,{},dict(metadata,extra=1),dict(metadata,energy=dict(metadata['energy'],solverMinusPhysicalJkg=0))):
        try:read_restart(data/'restart_flow.csv',bad,metadata,data/'mesh.su2')
        except ValueError:rejection+=1
        else:raise AssertionError('schema/energy accepted')
    save(data/'restart_flow.rah.json',metadata);save(out/'boundary_energy.json',records)
    save(out/'report.json',dict(passed=True,errors=err,characteristicClasses=sorted(classes),boundaryRows=12,wallRows=12,
        sourceBranches=2,initialNodes=24,savedConservativeAndSARowsBitwise=24,identityRejections=rejection,
        restartAdmission=True,actualRestartVerified=False,CFDAdmission=False,newQualifiedAreaM2=0))
    print(json.dumps(dict(passed=True,errors=err,identityRejections=rejection)))
if __name__=='__main__':main()
