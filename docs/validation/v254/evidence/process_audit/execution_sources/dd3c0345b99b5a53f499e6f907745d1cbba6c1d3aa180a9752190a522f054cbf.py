"""模块职责：独立DOP853连续积分及物性/质量账，区别数值闭合与物理资格。"""
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v254_validation';CLOSED=False


def audit():
    started=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    src=read(BASE/'prepare01/source.json');air=PhaseAwareViscousAir()
    production={b:read(BASE/'matlab02'/f'{b}.json') for b in ('base','fine','reference','scaled')}
    assert read(BASE/'matlab02/report.json')['passed']
    p0=src['stagnationPressurePa'];pinf=src['pressurePa'];T0=src['stagnationTemperatureK']
    R=src['radiusM'];Uinf=src['speedMps'];H=src['totalEnthalpyJkg'];cm=5.27*.4696
    stag=air.evaluate(p0,T0);wall=air.evaluate(p0,300.)
    muw=float(air.interpolated_transport(p0,300.)['dynamicViscosityPas'])
    mus=float(air.interpolated_transport(p0,T0)['dynamicViscosityPas'])
    kw=float(air.interpolated_transport(p0,300.)['conductivityWmK'])
    rw=float(wall['densityKgm3']);rs=float(stag['densityKgm3']);hw=float(wall['enthalpyJkg'])
    C0=rw*muw*Uinf;Pr=float(wall['cpJkgK'])*muw/kw
    a0=np.sqrt(2*(p0-pinf)/(rs*R*R))
    coeff=-2*T0*src['entropyDerivativeByMassAtZero']*cm*np.sqrt(rw*muw/2)
    a=brentq(lambda v:v*v-a0*a0-coeff*np.sqrt(v),a0,16*a0,xtol=1e-10)
    cp60,h60,_=air.thermal(60.);offset=float(h60-60*cp60);Hw=H-offset;zeta=(hw-offset)/Hw
    G7=(1.116+.411*.5**.686)/(1+.527*.5**.686)
    q0=.768*1.033/G7*Pr**(-.6)*np.sqrt(a0)*(rw*muw)**.1*(rs*mus)**.4*(H-hw)
    lookup=PchipInterpolator(src['sourceMassCoordinate'],src['sourceEntropyJkgK'],extrapolate=False)
    upper=src['sourceMassCoordinate'][-1];theta0=.001

    def state(angle,j,mode='coupled'):
        # 与MATLAB的对数I隐式梯形不同：无量纲J=I/R^3，显式连续ODE。
        psi=cm*np.sqrt(2*C0*R**3*j);assert 0<=psi<=upper
        entropy=float(lookup(psi)) if mode=='coupled' else src['sourceEntropyJkgK'][0]
        p=pinf+(p0-pinf)*np.cos(angle)**2
        wanted=entropy+air.R*np.log(p/100000.)
        te=brentq(lambda t:float(air.thermal(t)[2])-wanted,60.,T0,xtol=2e-11)
        st=air.evaluate(p,te);u=np.sqrt(2*(H-float(st['enthalpyJkg'])));rho=float(st['densityKgm3'])
        h=R*np.sin(angle);integ=j*R**3;proxy=2*(p0-pinf)*np.sin(angle)*np.cos(angle)/(R*rho*u)
        enthalpyRatio=(float(st['enthalpyJkg'])-offset)/Hw
        beta=2*integ*proxy*Uinf/(enthalpyRatio*(p/p0)*u*u*h*h)
        assert np.isfinite(beta) and beta>=0
        shape=G7*(1+.527*beta**.686)/(1.116+.411*beta**.686)*(1.1-.1625*enthalpyRatio+.0625*enthalpyRatio**2)*(.85+.15*enthalpyRatio-zeta)/(1-zeta)
        q=q0*(p/p0)*(u/Uinf)*h*shape/np.sqrt(4*a0/Uinf*integ)
        return dict(edgeTemperatureK=te,edgeSpeedMps=float(u),massThicknessM=float(psi/(rho*u*h)),
                    heatFluxWm2=float(q),massCoordinate=float(psi),pressurePa=float(p),edgeEntropyJkgK=entropy,betaPressureProxy=float(beta))

    target=np.asarray(src['targetAnglesRad']);results={};errors={};branches={}
    for mode,key in [('coupled','base'),('normal-shock-reference','reference')]:
        gradient=a if mode=='coupled' else a0
        initial=gradient*R*theta0**4/(4*Uinf)
        def rhs(theta,y):
            e=state(theta,y[0],mode)
            return [e['pressurePa']/p0*e['edgeSpeedMps']/Uinf*np.sin(theta)**2]
        sol=solve_ivp(rhs,[theta0,float(target.max())],[initial],method='DOP853',rtol=2e-10,atol=1e-18,dense_output=True,max_step=np.deg2rad(2))
        assert sol.success,sol.message
        rr=[state(t,j,mode) for t,j in zip(target,sol.sol(target)[0])];results[key]=rr
        for branch in (['base','fine'] if key=='base' else ['reference']):
            actual=production[branch];idx=np.array(actual['targetIndices'])-1;error={}
            for name in ('heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps'):
                ref=np.array([r[name] for r in rr]);y=np.array(actual[name])[idx]
                error[name]=float(np.max(abs(y/ref-1)))
            errors[branch]=error
        branches[key]=dict(acceptedODESteps=len(sol.t)-1,evaluations=sol.nfev,terminalAngleDeg=float(np.rad2deg(sol.t[-1])))
    save(out/'independent_targets.json',dict(targetAnglesRad=target.tolist(),branches=results))
    # 在实际保存态上独立核对H/S/代数质量，不把厚度伪装成完整剖面的数值积分。
    balance={}
    for branch,r in production.items():
        assert r['available'] and r['acceptedStations']==len(r['angleRad']) and not r['physicalQualification'] and not r['profileReadyForCFD']
        th=np.array(r['angleRad']);t=np.array(r['edgeTemperatureK']);p=np.array(r['pressurePa']);u=np.array(r['edgeSpeedMps'])
        mass=np.array(r['massCoordinate']);scale=r['options']['metricScale'];st=air.evaluate(p,t)
        h=np.array(r['metricM']);delta=np.array(r['massThicknessM'])
        energy=float(np.max(abs(st['enthalpyJkg']+.5*u*u-H)/abs(H)))
        entropy=float(np.max(abs(st['entropyJkgK']-r['edgeEntropyJkgK'])/max(1,abs(src['stagnationEntropyJkgK']))))
        matching=float(np.max(abs(mass[1:]/(st['densityKgm3'][1:]*u[1:]*h[1:]*delta[1:])-1)))
        expected=lookup(mass/scale) if branch!='reference' else np.full(len(th),src['sourceEntropyJkgK'][0])
        sourceError=float(np.max(abs(expected-r['edgeEntropyJkgK'])/max(1,abs(src['stagnationEntropyJkgK']))))
        leadingError=abs(r['jointLeadingVelocityGradient']/(a0 if branch=='reference' else a)-1)
        balance[branch]=dict(energyRelative=energy,entropyRelative=entropy,massRelative=matching,sourceEntropyRelative=sourceError,leadingGradientRelative=float(leadingError))
    # 只对首失败已接受前缀核对；失败试探站没有成为已接受结果。
    failed=read(BASE/'matlab01/base.json');n=failed['acceptedStations'];base=production['base'];prefix={}
    for key in ('heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps','integralI'):
        x=np.array(failed[key][:n]);y=np.array(base[key][:n]);prefix[key]=float(np.max(abs(x-y)/np.maximum(abs(x),1e-30)))
    # 原环真实标签尺度与拟合球帽尺度分别报告，实际几何未被替换。
    old=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'))
    fitted=R*np.sin(target[:46])*np.array(src['targetLabelJacobian']);actual=np.array(src['targetActualMetricM'])
    metricDifference=float(np.max(abs(fitted/actual-1)))
    ids=np.array(base['targetIndices'][:46])-1;refs=np.array(production['reference']['targetIndices'][:46])-1
    correction={key:dict(min=float(np.min(np.array(base[key])[ids]/np.array(production['reference'][key])[refs]-1)),
                         max=float(np.max(np.array(base[key])[ids]/np.array(production['reference'][key])[refs]-1)))
                for key in ('heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps')}
    numerical=all(v<=2e-3 for e in errors.values() for v in e.values()) and all(v<=1e-8 for e in balance.values() for v in e.values())
    report=dict(passed=bool(numerical),independentMethod='DOP853 dimensionless continuous integral; Brent entropy recovery; shared previously audited mass/entropy source',
        comparisonRelativeErrors=errors,balanceErrors=balance,solver=branches,preservedPrefixRelativeDifference=prefix,
        actualRingPositions=46,extraTerminalPosition=1,allRingIntegralStatesAvailable=True,missingOldProfilePositions=(np.flatnonzero(~np.isfinite(old['base_yM'][:,0]))).tolist(),
        actualProfilesFilled=0,maximumFittedVsActualMetricRelativeDifference=metricDifference,
        ringCorrectionVsFixedEntropy=correction,betaRange=[min(base['betaPressureProxy']),max(base['betaPressureProxy'])],
        maximumRingMomentumDefectPam=float(np.max(abs(np.array(base['alongMomentumDefectPam'])[ids]))),
        sourceSHA256=src['sourceSHA256'],matlabReportSHA256=sha(BASE/'matlab02/report.json'),
        originalFaces=19752,originalAreaM2=1.9646292565715553,physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,
        elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report)
    assert numerical,errors
    assert max(prefix.values())<1e-8
    print(json.dumps({k:report[k] for k in ('passed','comparisonRelativeErrors','ringCorrectionVsFixedEntropy','elapsedSeconds')}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v254 fixed independent audit closed; follow NEXT_DECISION')
    audit()
