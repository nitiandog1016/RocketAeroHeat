"""模块职责：恢复观察与独立温度变量沿程复算；积分账区分离散截断和求解误差。"""
from pathlib import Path
import sys, time
import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT, read, save, sha

BASE=ROOT/'outputs/v252_validation'
SOURCE=BASE/'matlab01'


def recover():
    """模块一：只记当前可观察事实，不补造丢失的父监督退出码/时间/峰值。"""
    import psutil
    out=BASE/'recovery01';out.mkdir(exist_ok=False)
    admission=read(BASE/'process_matlab01/admission.json')
    live=read(BASE/'process_matlab01/live.json')
    assert not psutil.pid_exists(live['pid']), 'Original PID still exists; inspect before recovery'
    active=[]
    for p in psutil.process_iter(['pid','name','cmdline']):
        if 'matlab' in (p.info['name'] or '').lower():active.append(p.info)
    assert not active,active
    matches={p:sha(p)==h for p,h in admission['inputSHA256'].items()}
    assert all(matches.values()),matches
    first={p.name:sha(p) for p in SOURCE.glob('*.json')}
    assert len(first)==8 and read(SOURCE/'report.json')['passed']
    second={p.name:sha(p) for p in SOURCE.glob('*.json')};assert first==second
    save(out/'report.json',dict(passed=True,observationOnly=True,originalSupervisorReceiptPresent=False,
        originalExitCode=None,originalElapsedSeconds=None,originalPeakRSSBytes=None,
        originalWholeRunInputStabilityKnown=False,originalProcessCurrentlyAbsent=True,
        inputsMatchAdmissionNow=matches,completeMatlabTestReportPresent=True,outputSHA256=first,
        outputsStableAtRecovery=True,matlabBudgetConservativelyDebitedSeconds=300,
        noMatlabRerun=True,physicalQualification=False))


def decoded(r,air):
    y=np.array(r['numericalProfile']['states'],float)
    y[3]=60+1940/(1+np.exp(-y[3]))
    return PchipInterpolator(r['numericalProfile']['eta'],y,axis=1)


def parameters(r,air):
    e=r['edge'];s=r['station'];Te=e['temperatureK'];U=e['speedMps'];p=e['pressurePa']
    cp,h,_=air.thermal(60.);offset=float(h-60*cp);H=e['totalEnthalpyJkg']-offset
    mu=float(np.interp(Te,air.T,air.mu));rho=p/(air.R*Te)
    K=(np.sqrt(mu/rho/(2*s['velocityGradientPerSecond'])) if U==0 else
       np.sqrt(2*s['xi'])/(rho*U*s['metricM']))
    return Te,U,p,offset,H,mu,K


def properties(y,air,Te,U,offset,H,mu_e):
    T=y[3];cp,h,_=air.thermal(T)
    mu=np.interp(T,air.T,air.mu);k=np.interp(T,air.T,air.kappa)
    ell=mu/mu_e*Te/T;Pr=cp*mu/k;g=(h+.5*U*U*y[1]**2-offset)/H
    return cp,mu,k,ell,Pr,g


def independent_chain(chain,air):
    """模块二：SciPy直接解物理T，无MATLAB logit/下游解；每站继承自己的独立上游。"""
    rows=[];records=chain['profiles'];past=decoded(records[0],air)
    for j,r in enumerate(records[1:],1):
        prev=records[j-1];Te,U,p,offset,H,mu_e,K=parameters(r,air)
        A=2*r['station']['xi']/(r['station']['xi']-prev['station']['xi']);beta=r['station']['beta']
        def ode(x,y):
            cp,mu,k,ell,Pr,g=properties(y,air,Te,U,offset,H,mu_e)
            old=past(x);_,oldh,_=air.thermal(old[3])
            oldg=(oldh+.5*prev['edge']['speedMps']**2*old[1]**2-offset)/H
            fp=y[2]/ell;gp=Pr/ell*(y[4]-(1-1/Pr)*U*U/H*y[1]*y[2])
            tp=Pr/(ell*cp)*(H*y[4]-U*U*y[1]*y[2]);df=y[0]-old[0]
            return np.array([y[1],fp,-y[0]*fp+beta*(y[1]**2-y[3]/Te)+A*(y[1]*(y[1]-old[1])-fp*df),
                             tp,-y[0]*gp+A*(y[1]*(g-oldg)-gp*df)])
        def bc(a,b):return np.array([a[0],a[1],(a[3]-r['wallTemperature'])/Te,b[1]-1,(b[3]-Te)/Te])
        x=np.linspace(0,12,151)
        sol=solve_bvp(ode,bc,x,past(x),tol=1e-7,max_nodes=12000)
        assert sol.success,(j,sol.message)
        nodes=np.linspace(0,12,1201);Y=sol.sol(nodes);air.evaluate(p,Y[3])
        cp,mu,k,ell,Pr,g=properties(Y,air,Te,U,offset,H,mu_e)
        fp=Y[2]/ell;gp=Pr/ell*(Y[4]-(1-1/Pr)*U*U/H*Y[1]*Y[2]);fs=max(abs(fp));gs=max(.005*max(abs(gp)),1e-8)
        def defect(x):
            y=sol.sol(x);cp,mu,k,l,Pr,g=properties(y,air,Te,U,offset,H,mu_e)
            gp=Pr/l*(y[4]-(1-1/Pr)*U*U/H*y[1]*y[2])
            return np.max(np.array([abs(1-y[1])/.005,abs(1-g)/.005,abs(y[3]/Te-1)/.005,
                abs(Te/y[3]-1)/.005,abs(y[2]/l)/fs/.005,abs(gp)/gs]),axis=0)
        bad=np.flatnonzero(defect(nodes)>1);assert len(bad) and nodes[bad[-1]+1]<=11
        ee=brentq(lambda x:float(defect(x)-1),nodes[bad[-1]],nodes[bad[-1]+1],xtol=1e-12)
        xx=np.linspace(0,ee,2401);yy=sol.sol(xx);F=yy[1]
        integ=lambda v:float(K*simpson(v,x=xx))
        delta=integ(yy[3]/Te);mass=integ(F)
        wall=sol.sol(0);der=ode(0,wall);_,mu,k,_,_,_=properties(wall,air,Te,U,offset,H,mu_e)
        calc=dict(heatFlux=float(k*der[3]/(K*wall[3]/Te)),wallShearStress=float(mu*U*der[1]/(K*wall[3]/Te)),
            delta=delta,deltaStar=delta-mass,massThickness=mass,thetaPhysical=integ(F*(1-F)),energyThickness=integ(F*(1-F*F)))
        err={key:abs(val-r[key])/max(abs(val),1e-12) for key,val in calc.items()}
        xp=np.array(r['profile']['eta']);py=sol.sol(xp)
        err['temperatureScaled']=float(max(abs(py[3]-np.array(r['profile']['temperature'])))/Te)
        err['velocityRatio']=float(max(abs(py[1]-np.array(r['profile']['states'])[1])))
        assert all(v<(.002 if k in ('delta','deltaStar','massThickness') else .0002) for k,v in err.items()),(j,err)
        rows.append(dict(station=j,thetaRad=r['station']['thetaRad'],values=calc,relativeErrors=err,nodes=len(sol.x)))
        # 上游物理T的连续解成为下一站入口，不借用下游MATLAB剖面。
        past=sol.sol
    return rows


def balances(chain,air):
    """模块三：对完整eta域独立积分；显式列出后退差分乘积法则缺项。"""
    x=np.linspace(0,12,6001);integral=lambda a:float(simpson(a,x=x));rows=[]
    for j,r in enumerate(chain['profiles'][1:],1):
        prev=chain['profiles'][j-1];Te,U,p,offset,H,mu_e,K=parameters(r,air)
        y=decoded(r,air)(x);old=decoded(prev,air)(x);F=y[1];Fo=old[1]
        _,_,_,ell,Pr,g=properties(y,air,Te,U,offset,H,mu_e)
        _,ho,_=air.thermal(old[3]);go=(ho+.5*prev['edge']['speedMps']**2*Fo**2-offset)/H
        A=2*r['station']['xi']/(r['station']['xi']-prev['station']['xi']);beta=r['station']['beta']
        # 完整边界项保留F(L)/g(L)，不以数值近1的假设消除外边界通量。
        dfL=y[0,-1]-old[0,-1]
        mom=y[2,-1]+y[0,-1]*F[-1]-integral(F*F)+beta*integral(y[3]/Te-F*F)-A*(integral(F*F-Fo*Fo)-F[-1]*dfL)
        energy=y[4,-1]+y[0,-1]*g[-1]-integral(F*g)-A*(integral(F*g-Fo*go)-g[-1]*dfL)
        tm=A*integral((F-Fo)**2);te=A*integral((F-Fo)*(g-go))
        scaleM=max(abs(y[2,0]),1e-12);scaleE=max(abs(y[4,0]),1e-12)
        rows.append(dict(station=j,momentumRawRelative=abs(mom-y[2,0])/scaleM,
            energyRawRelative=abs(energy-y[4,0])/scaleE,momentumTruncationRelative=tm/scaleM,
            energyTruncationRelative=te/scaleE,momentumCorrectedRelative=abs(mom-tm-y[2,0])/scaleM,
            energyCorrectedRelative=abs(energy-te-y[4,0])/scaleE,
            wallMomentumFlux=float(y[2,0]),wallEnthalpyFlux=float(y[4,0])))
    return rows


def audit():
    """模块四：独立守恒/原面映射/源身份；数值通过与物理资格分别出具。"""
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    recovery=read(BASE/'recovery01/report.json');assert recovery['passed']
    assert all(sha(SOURCE/name)==h for name,h in recovery['outputSHA256'].items())
    chains={k:read(SOURCE/(k+'.json')) for k in ('base','fine','history')};air=PhaseAwareViscousAir()
    independent=independent_chain(chains['base'],air);save(out/'independent_chain.json',independent)
    allbalances={name:balances(c,air) for name,c in chains.items()};save(out/'integral_balances.json',allbalances)
    maxima={name:{k:max(row[k] for row in rows) for k in ('momentumRawRelative','energyRawRelative',
            'momentumCorrectedRelative','energyCorrectedRelative')} for name,rows in allbalances.items()}
    corrected=all(v[k]<=.01 for v in maxima.values() for k in ('momentumCorrectedRelative','energyCorrectedRelative'))
    raw=all(v[k]<=.01 for v in maxima.values() for k in ('momentumRawRelative','energyRawRelative'))
    geometry=read(SOURCE/'geometry_viewer.json');V=np.array(geometry['vertices']);F=np.array(geometry['faces'])
    tri=V[F];lengths=np.linalg.norm(tri-np.roll(tri,1,axis=1),axis=2);s=lengths.sum(axis=1)/2
    area=np.sqrt(s*np.prod(s[:,None]-lengths,axis=1));stored=np.array(geometry['faceArea'])
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<2e-13
    assert max(abs(area-stored)/np.maximum(area,1e-30))<1e-10
    mappings={}
    for name,c in chains.items():
        h=read(SOURCE/('surface_'+name+'.json'));valid=np.array(h['validFace'],bool)
        lo=np.array(h['lowerStation'])[valid]-1;hi=np.array(h['upperStation'])[valid]-1;w=np.array(h['upperWeight'])[valid]
        assert valid.sum()==1306 and np.all(hi==lo+1) and np.all((w>=0)&(w<=1))
        checked=0
        for key in ('heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness','energyThickness'):
            values=np.array(h[key],float);assert np.all(np.isnan(values[~valid]))
            table=np.array([r[key] for r in c['profiles']]);ref=(1-w)*table[lo]+w*table[hi]
            assert np.array_equal(ref,values[valid]),(name,key)
            checked+=len(ref)
        power=float(sum(area[valid]*np.array(h['heatFlux'],float)[valid]));assert abs(power-h['summary']['signedPowerW'])<1e-8
        mappings[name]=dict(candidateFaces=int(valid.sum()),candidateAreaM2=float(sum(area[valid])),signedPowerW=power,checkedValues=checked)
    report=dict(passed=True,independentChainPassed=True,independentStations=len(independent),
        independentMaximumRelativeError=max(v for row in independent for v in row['relativeErrors'].values()),
        integralBalances=maxima,correctedDiscreteIntegralGate=corrected,rawConservativeIntegralGate=raw,
        integralGateThreshold=.01,allNumericalGatesPassed=corrected and raw and read(SOURCE/'report.json')['streamwiseRefinementGate'],
        sameEquationsNumericalVerificationOnly=True,experimentalValidation=False,
        originalMatlabSupervisorTerminalKnown=False,noMatlabRerun=True,
        originalFaces=len(F),totalAreaM2=float(sum(area)),mappings=mappings,physicalQualification=False,newQualifiedAreaM2=0,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(report)


if __name__=='__main__':
    {'recover':recover,'audit':audit}[sys.argv[1]]()
