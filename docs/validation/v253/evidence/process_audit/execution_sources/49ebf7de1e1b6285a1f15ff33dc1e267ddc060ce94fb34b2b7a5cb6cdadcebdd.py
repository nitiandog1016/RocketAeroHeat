"""模块职责：独立压力比激波/直接积分、原拓扑与实际剖面质量审核，资格不升级。"""
import sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v253_validation';SRC=BASE/'inlet02'


def audit():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    air=PhaseAwareViscousAir();meta=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json');c=meta['condition']
    p1=c['pressurePa'];T1=c['temperatureK'];U1=c['speedMps'];rho1=p1/(air.R*T1)
    h1=float(air.thermal(T1)[1]);H=h1+.5*U1*U1
    table=read(BASE/'inlet01/shock_mass_source.json');rows=read(SRC/'inlets.json');a=np.load(str(SRC/'inlet_profiles.npz'));r=read(SRC/'report.json')
    def shock(u):
        # 与生产z变量不同，直接以压强比求Hugoniot压缩根；不调用生产波系函数。
        v=U1*u;upper=1+rho1*v*v/p1
        def values(pr):
            un=v-p1*(pr-1)/(rho1*v);T=p1*pr*un/(rho1*v*air.R)
            return un,T
        def residual(pr):
            un,T=values(pr)
            if not 60<=T<=2000 or un<=0:return np.nan
            return float(air.thermal(T)[1])+.5*un*un-h1-.5*v*v
        x=np.linspace(1.000001,upper*(1-1e-8),100)
        pair=next(((lo,hi) for lo,hi in zip(x[:-1],x[1:]) if residual(lo)>0 and residual(hi)<0),None)
        if pair is None:raise ValueError('INDEPENDENT_SHOCK_BRACKET')
        pr=brentq(residual,*pair,xtol=1e-12);un,T=values(pr);p=p1*pr;s=air.evaluate(p,T)
        return p,T,un,float(s['entropyJkgK'])
    p2,_,_,_=shock(1.);R=table['radiusM'];Q=rho1*U1*U1
    probe=np.array([1e-3,2.5e-4,6.25e-5]);d=np.array([(p2-shock(np.sqrt(1-w))[0])/w for w in probe]);D=(4*d[2]-d[1])/3
    A=(p2-p1)/(Q+2*D);assert abs(A/table['A']-1)<1e-8
    q0=1e-6;init=[1-A*q0,q0/2-(A+1/6)*q0*q0/4];last=table['thetaRad'][-1]
    def rhs(q,y):
        th=np.sqrt(q);u,psi=y;pb=p1+(p2-p1)*np.cos(th)**2
        return [(pb-shock(u)[0])/Q*np.sin(th)/(2*th*psi),u*np.sin(th)/(2*th)]
    ref=solve_ivp(rhs,[q0,last*last],init,method='RK45',rtol=1e-10,atol=[1e-12,1e-14],dense_output=True)
    assert ref.success
    th=np.array(table['thetaRad'][1:]);values=ref.sol(th*th)
    tableU=np.array(table['sinShockAngle'][1:]);tableMass=np.array(table['massCoordinate'][1:])
    tableError=max(float(max(abs(values[0]-tableU))),float(max(abs(values[1]*rho1*U1*R*R/tableMass-1))))
    assert tableError<=2e-5,tableError
    # 全原面独立构造边-面关联，不调用生产common_ring。
    V=a['originalVertices'];F=a['originalFaces'];region=read(ROOT/'outputs/v228_validation/prepare03/region.json')
    allEdges=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);owners=np.tile(np.arange(len(F)),3)
    keys,inv,count=np.unique(np.sort(allEdges,axis=1),axis=0,return_inverse=True,return_counts=True);assert np.all(count==2)
    order=np.argsort(inv,kind='stable');facePairs=owners[order].reshape(-1,2)
    nose=np.isin(facePairs,np.array(region['retainedNoseCandidateIDs'])-1);target=np.isin(facePairs,np.array(region['originalFaceIDs'])-1)
    shared=keys[np.any(nose,axis=1)&np.any(target,axis=1)]
    assert {tuple(e) for e in shared}=={tuple(sorted(e)) for e in a['ringEdges']} and len(shared)==46
    lengths=np.linalg.norm(V[F]-np.roll(V[F],1,axis=1),axis=2);semi=lengths.sum(1)/2;areas=np.sqrt(semi*np.prod(semi[:,None]-lengths,axis=1))
    assert abs(areas.sum()-1.9646292565715553)<2e-13
    h=np.einsum('ij,ij->i',np.cross(a['ringLabelDerivative'],a['unitNormal']),a['unitTangent'])
    assert np.max(abs(h-a['metricM']))<1e-15
    # 小段精确原边方向上的风轴方位差分，独立核对标签Jacobian。
    center=np.array(meta['sphereCenter']);flow=np.array(c['flowDirection']);e1=np.array([0.,1.,0.]);e2=np.cross(-flow,e1)
    xp=a['wallXYZ']-center+1e-5*a['ringLabelDerivative'];xm=a['wallXYZ']-center-1e-5*a['ringLabelDerivative']
    angleP=np.arctan2(xp@e2,xp@e1);angleM=np.arctan2(xm@e2,xm@e1)
    derivative=np.angle(np.exp(1j*(angleP-angleM)))/(2e-5)
    labelError=float(max(abs(derivative/a['labelJacobian']-1)));assert labelError<1e-8
    checks=[];missing=[];massErrors=[];waveErrors=[];edgeErrors=[];coordinate=np.array(table['massCoordinate']);shockLookup=PchipInterpolator(coordinate,table['sinShockAngle'])
    for row in rows:
        branch=row['branch'];i=row['index']
        if not row['available']:
            missing.append([branch,i]);assert all(np.isnan(a[branch+'_'+k][i]).all() for k in ('yM','temperatureK','speedMps','densityKgm3','physicalConservative'))
            assert row['matchedEdge'] is None and not row['profileReadyForCFD'];continue
        y=a[branch+'_yM'][i];T=a[branch+'_temperatureK'][i];u=a[branch+'_speedMps'][i];rho=a[branch+'_densityKgm3'][i];cons=a[branch+'_physicalConservative'][i]
        pressure=row['pressurePa'];assert max(abs(rho*air.R*T/pressure-1))<1e-14
        # 第二种法向积分（梯形）与生产Simpson比较，严格保持原1e-4门限。
        mass=float(np.trapz(rho*u,x=y));err=abs(mass/row['massFluxPerWidthKgMs']-1);massErrors.append(err);assert err<=1e-4
        pvec=np.outer(rho*u,a['unitTangent'][i]);assert np.array_equal(pvec,cons[:,1:4])
        internal=cons[:,4]/rho-.5*np.sum((cons[:,1:4]/rho[:,None])**2,axis=1)
        assert max(abs(internal-(air.thermal(T)[1]-air.R*T)))<1e-8
        matched=row['matchedEdge'];sinA=matched['sourceSinShockAngle'];pj,Tj,un,sj=shock(sinA)
        v=U1*sinA;massjump=pj/(air.R*Tj)*un;momjump=pj+pj/(air.R*Tj)*un*un;energy=float(air.thermal(Tj)[1])+.5*un*un
        wave=max(abs(massjump/(rho1*v)-1),abs(momjump/(p1+rho1*v*v)-1),abs((energy-h1-.5*v*v)/H));waveErrors.append(wave);assert wave<1e-8
        target=sj+air.R*np.log(pressure/1e5);Te=brentq(lambda t:float(air.thermal(t)[2])-target,60.,2000.,xtol=1e-10);Ue=np.sqrt(2*(H-float(air.thermal(Te)[1])))
        edge=max(abs(Te/matched['temperatureK']-1),abs(Ue/matched['speedMps']-1));edgeErrors.append(edge);assert edge<1e-8
        # 先独立倒查实际插值算子的质量，再比较不使用压力插值表的直接ODE。
        recoveredMass=brentq(lambda psi:float(shockLookup(psi))-sinA,coordinate[0],coordinate[-1],xtol=1e-16)
        lookupError=abs(recoveredMass/matched['sourceMassCoordinate']-1);assert lookupError<1e-8
        thetaRef=brentq(lambda th:float(ref.sol(th*th)[0])-sinA,.001,last,xtol=1e-12)
        directMass=float(ref.sol(thetaRef*thetaRef)[1])*rho1*U1*R*R
        directError=abs(directMass/matched['sourceMassCoordinate']-1);assert directError<2e-5,directError
        checks.append(dict(branch=branch,index=i,massIntegrationRelativeDifference=err,lookupMassRelativeError=lookupError,
            directODEMassRelativeDifference=directError,independentEdgeRelativeError=edge,waveConservationRelativeError=wave))
    assert len(checks)==123 and len(missing)==15
    # 以真实原始源重新核验插值身份，保证缺口源于范围而非丢行。
    for branch in ('base','fine','history'):
        chain=read(ROOT/'outputs/v252_validation/matlab01'/(branch+'.json'));assert sha(ROOT/'outputs/v252_validation/matlab01'/(branch+'.json'))==r['sourceSHA256']['outputs\\v252_validation\\matlab01\\'+branch+'.json']
        lo=chain['stations'][0]['thetaRad'];hi=chain['stations'][-1]['thetaRad']
        assert np.flatnonzero((a['thetaRad']<lo)|(a['thetaRad']>hi)).tolist()==[i for b,i in missing if b==branch]
    save(out/'checks.json',checks)
    report=dict(passed=True,commonRingEdges=46,availableProfileInlets=123,missingProfileInlets=15,
        physicalConservativeSamples=123*601,sourceTableDirectODEMaximumDifference=tableError,labelJacobianRelativeDifference=labelError,
        maximumMassIntegrationRelativeDifference=max(massErrors),maximumWaveConservationRelativeError=max(waveErrors),
        maximumIndependentEdgeRelativeDifference=max(edgeErrors),maximumDirectODEMassRelativeDifference=max(x['directODEMassRelativeDifference'] for x in checks),
        onePassJointInletAccepted=False,profileReadyForCFD=False,physicalQualification=False,newQualifiedAreaM2=0,
        missingProfilePositions=missing,inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [SRC/'inlets.json',SRC/'inlet_profiles.npz',BASE/'inlet01/shock_mass_source.json']},
        sourceKind='same-equations independent pressure-coordinate shock / RK45 direct source / original topology',
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(report)


if __name__=='__main__':audit()
