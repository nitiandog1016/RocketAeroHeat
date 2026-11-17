"""模块职责：直接J(theta)的DOP853与独立稳定量热恢复，审核全原面/起步/环及资格边界。"""
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v255_validation';CLOSED=False


def audit():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    source=read(BASE/'prepare01/source.json');geo=read(BASE/'prepare01/geometry.json');air=PhaseAwareViscousAir()
    src=BASE/('matlab02' if (BASE/'matlab02/report.json').exists() else 'matlab01');assert read(src/'report.json')['passed']
    data={k:read(src/f'{k}.json') for k in ['base','fine','seed','scaled']};base=data['base']
    p0=source['stagnationPressurePa'];pinf=source['pressurePa'];T0=source['stagnationTemperatureK'];Uinf=source['speedMps'];R=source['radiusM'];H=source['totalEnthalpyJkg']
    wall=air.evaluate(p0,300.);stag=air.evaluate(p0,T0);trw=air.interpolated_transport(p0,300.);trs=air.interpolated_transport(p0,T0)
    rw=float(wall['densityKgm3']);rs=float(stag['densityKgm3']);muw=float(trw['dynamicViscosityPas']);mus=float(trs['dynamicViscosityPas']);Pr=float(wall['cpJkgK']*muw/trw['conductivityWmK']);hw=float(wall['enthalpyJkg'])
    cm=5.27*.4696;C0=rw*muw*Uinf;a0=np.sqrt(2*(p0-pinf)/(rs*R*R));A=1-pinf/p0
    entropy=np.asarray(source['sourceEntropyJkgK']);lookup=PchipInterpolator(source['sourceMassCoordinate'],entropy-entropy[0],extrapolate=False)
    slope=float(lookup.derivative()(0));coef=-2*T0*slope*cm*np.sqrt(rw*muw/2)
    a=brentq(lambda x:x*x-a0*a0-coef*np.sqrt(x),a0,16*a0,xtol=1e-10)
    cp60,h60,_=air.thermal(60.);offset=float(h60-60*cp60);Hw=H-offset;zeta=(hw-offset)/Hw;G7=(1.116+.411*.5**.686)/(1+.527*.5**.686)
    q0=.768*1.033/G7*Pr**(-.6)*np.sqrt(a0)*(rw*muw)**.1*(rs*mus)**.4*(H-hw)
    # 与MATLAB向下累计表不同：逐个原Cp区间直接计算与[T0-drop,T0]的交集。
    def thermal_drop(drop):
        lo=np.maximum(air.T[:-1],T0-drop);hi=np.minimum(air.T[1:],T0);mask=hi>lo;lo=lo[mask];hi=hi[mask]
        cpLo=air.Cp[:-1][mask]+air.slopes[mask]*(lo-air.T[:-1][mask]);width=hi-lo
        # 温降极小时hi-lo会消去；同一个原区间时直接用传入drop。
        if len(width)==1:width=np.array([drop])
        cpHi=air.Cp[:-1][mask]+air.slopes[mask]*(hi-air.T[:-1][mask]);slopeT=air.slopes[mask]
        dh=np.sum(width*(cpHi-.5*slopeT*width))
        ds=np.sum(-(cpHi-slopeT*hi)*np.log1p(-width/hi)+slopeT*width)
        return float(dh),float(ds)
    def edge(theta,j):
        psi=cm*np.sqrt(2*C0*R**3*j);assert 0<=psi<=source['sourceMassCoordinate'][-1]
        de=float(lookup(psi));lp=np.log1p(-A*np.sin(theta)**2);p=p0*np.exp(lp);wanted=-de-air.R*lp
        dt=brentq(lambda v:thermal_drop(v)[1]-wanted,0,T0-60,xtol=2e-14)
        dh,_=thermal_drop(dt);te=T0-dt;state=air.evaluate(p,te);u=np.sqrt(2*dh);rho=float(state['densityKgm3']);metric=R*np.sin(theta)
        integral=R**3*j;proxy=(p0-pinf)*np.sin(2*theta)/(R*rho*u);teh=(H-dh-offset)/Hw
        beta=2*integral*proxy*Uinf/(teh*(p/p0)*u*u*metric*metric);assert np.isfinite(beta) and beta>=0
        shape=G7*(1+.527*beta**.686)/(1.116+.411*beta**.686)*(1.1-.1625*teh+.0625*teh**2)*(.85+.15*teh-zeta)/(1-zeta)
        q=q0*(p/p0)*(u/Uinf)*metric*shape/np.sqrt(4*a0/Uinf*integral)
        return dict(edgeTemperatureK=te,edgeSpeedMps=float(u),heatFluxWm2=float(q),massThicknessM=float(psi/(rho*u*metric)),
                    pressurePa=float(p),massCoordinate=float(psi),integralI=float(integral),edgeEntropyJkgK=float(entropy[0]+de))
    # 独立较早有限种子，直接积分J，不使用W/logθ或MATLAB任何积分轨迹。
    thetaSeed=5e-5;initial=a*R*thetaSeed**4/(4*Uinf)
    def rhs(theta,y):
        e=edge(theta,y[0]);return [e['pressurePa']/p0*e['edgeSpeedMps']/Uinf*np.sin(theta)**2]
    sol=solve_ivp(rhs,[thetaSeed,np.pi/3],[initial],method='DOP853',rtol=1e-11,atol=1e-28,max_step=np.deg2rad(1),dense_output=True)
    assert sol.success and sol.t[-1]==np.pi/3,sol.message
    theta=np.asarray(base['angleRad']);rows=[edge(t,j) for t,j in zip(theta[1:],sol.sol(theta[1:])[0])]
    # 驻点解析值独立重建；零速度直接检查为零而不使用相对除法。
    b0=a0*a0/(2*a*a);origin=dict(edgeTemperatureK=T0,edgeSpeedMps=0.,
        heatFluxWm2=q0*np.sqrt(a/a0)*G7*(1+.527*b0**.686)/(1.116+.411*b0**.686),
        massThicknessM=cm*np.sqrt(rw*muw/(2*a))/rs,pressurePa=p0,massCoordinate=0.,integralI=0.,edgeEntropyJkgK=float(entropy[0]))
    rows=[origin]+rows;comparisons={}
    for name in ['base','fine','seed','scaled']:
        r=data[name];indexes=np.asarray(r['targetIndices'])-1;targets=np.asarray(r['targetAnglesRad'])
        refIndex=np.searchsorted(theta,targets);assert np.array_equal(theta[refIndex],targets)
        err={}
        for key in ['heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps']:
            values=np.array([x[key] for x in rows])[refIndex];actual=np.asarray(r[key])[indexes];good=values!=0
            assert np.all(actual[~good]==0);err[key]=float(np.max(abs(actual[good]/values[good]-1)))
        comparisons[name]=err
    save(out/'independent_samples.json',dict(angleRad=theta.tolist(),states=rows,source='direct J(theta) DOP853'))
    # 保存每个原生步和连续扩展系数，后续定向审核无需重跑已过ODE。
    segments=sol.sol.interpolants
    np.savez_compressed(out/'independent_continuous.npz',theta=sol.t,J=sol.y[0],segmentStart=np.array([s.t_old for s in segments]),
        segmentStep=np.array([s.h for s in segments]),segmentValue=np.array([s.y_old for s in segments]),denseF=np.array([s.F for s in segments]))
    balances={}
    for name,r in data.items():
        p=np.asarray(r['pressurePa']);t=np.asarray(r['edgeTemperatureK']);u=np.asarray(r['edgeSpeedMps']);mass=np.asarray(r['massCoordinate']);sc=r['options']['metricScale']
        st=air.evaluate(p,t);metric=np.asarray(r['metricM']);delta=np.asarray(r['massThicknessM'])
        errH=float(np.max(abs(st['enthalpyJkg']+.5*u*u-H)/abs(H)));errS=float(np.max(abs(st['entropyJkgK']-r['edgeEntropyJkgK'])/abs(entropy[0])))
        errM=float(np.max(abs(mass[1:]/(st['densityKgm3'][1:]*u[1:]*metric[1:]*delta[1:])-1)))
        errSource=float(np.max(abs(st['entropyJkgK']-(entropy[0]+lookup(mass/sc)))/abs(entropy[0])))
        balances[name]=dict(energy=errH,entropy=errS,mass=errM,sourceEntropy=errSource)
    # 独立几何：半空间/二阶锥不等式验证顶点，再检查固定内部重心样本。
    raw=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'));meta=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json')
    ids=np.array(geo['candidateFaceIDs'])-1;tri=raw['originalVertices'][raw['originalFaces'][ids]]-np.array(meta['sphereCenter']);axis=-np.asarray(meta['condition']['flowDirection']);axis/=np.linalg.norm(axis)
    projection=tri@axis;perpendicular=np.linalg.norm(np.cross(tri,axis),axis=2)
    admitted=np.all((projection>0)&(perpendicular<=np.sqrt(3)*projection),axis=1)
    assert np.array_equal(ids[admitted]+1,geo['coneAdmittedFaceIDs'])
    for weights in [[1/3]*3,[.8,.1,.1],[.1,.8,.1],[.1,.1,.8],[.5,.5,0],[0,.5,.5],[.5,0,.5]]:
        points=np.einsum('ijk,j->ik',tri[admitted],weights);proj=points@axis
        assert np.all(np.linalg.norm(np.cross(points,axis),axis=1)<=np.sqrt(3)*proj)
    assert len(raw['originalFaces'])==19752 and abs(raw['originalFaceArea'].sum()-1.9646292565715553)<2e-13
    old=read(ROOT/'outputs/v254_validation/surface_audit01/report.json');recovered=sorted(set(old['withheldFaceIDs'])&set(geo['coneAdmittedFaceIDs']))
    numerical=all(v<=.002 for e in comparisons.values() for v in e.values()) and all(v<=1e-8 for e in balances.values() for v in e.values())
    report=dict(passed=bool(numerical),comparisonRelativeErrors=comparisons,balanceErrors=balances,fullSamplesCompared=len(theta),
        admittedOriginalFaces=int(admitted.sum()),recoveredPreviouslyWithheldFaceIDs=recovered,originalFaces=19752,originalAreaM2=1.9646292565715553,
        nativeODEAcceptedSteps=len(sol.t)-1,nativeEvaluations=sol.nfev,independentMethod='direct J(theta) DOP853 / independent Cp-interval drop integrals / Brent / original mesh cone inequalities',
        seedAngleRad=thetaSeed,zeroAngleAnalyticLimitOnly=True,belowPositiveSeedNumericallySolved=False,
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,actualProfileSolved=False,
        matlabFolder=src.name,sourceSHA256=sha(BASE/'prepare01/source.json'),matlabSHA256=sha(src/'base.json'),elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);assert numerical,comparisons
    print(json.dumps({k:report[k] for k in ['passed','comparisonRelativeErrors','admittedOriginalFaces','elapsedSeconds']}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v255 fixed full-domain audit closed; follow NEXT_DECISION')
    audit()
