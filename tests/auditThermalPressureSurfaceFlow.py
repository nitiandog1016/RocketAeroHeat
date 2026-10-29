"""独立复算连续NASA7积分、平面守恒轨迹、真实切锥及锚点约定转换。"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from scipy.io import loadmat
from scipy.optimize import brentq

from auditPressureDrivenSurfaceFlow import audit_case, rows

ROOT=Path(__file__).resolve().parents[1]
MASS=np.array([.755,.232,.013])
SPECIES_R=8.31446261815324/np.array([.0280134,.0319988,.039948])
LOW=np.array([[3.53100528,-1.23660987e-4,-5.02999433e-7,2.43530612e-9,-1.40881235e-12,-1046.97628,2.96747468],
 [3.78245636,-2.99673416e-3,9.84730201e-6,-9.68129509e-9,3.24372837e-12,-1063.94356,3.65767573],
 [2.5,0,0,0,0,0,0]])
HIGH=np.array([[2.95257626,1.39690040e-3,-4.92631603e-7,7.86010367e-11,-4.60755321e-15,-923.948688,5.87188762],
 [3.28253784,1.48308754e-3,-7.57966669e-7,2.09470555e-10,-2.16717794e-14,-1088.45772,5.45323129],
 [2.5,0,0,0,0,0,0]])


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def polynomial(T,coefficients):
    """独立展开各组分多项式，再质量加权；不调用项目MATLAB物性。"""
    cp=coefficients@np.array([1,T,T*T,T**3,T**4,0,0])
    h=coefficients@np.array([T,T*T/2,T**3/3,T**4/4,T**5/5,1,0])
    s=coefficients@np.array([np.log(T),T,T*T/2,T**3/3,T**4/4,0,1])
    return np.array([np.dot(MASS*SPECIES_R,x) for x in [cp,h,s]])


JUMP=polynomial(1000.,HIGH)-polynomial(1000.,LOW)


def air(T,continuous=True):
    assert 80<=T<=2500
    if T<150:
        cp,h,s=polynomial(150.,LOW)
        value=np.array([cp,h+cp*(T-150),s+cp*np.log(T/150)])
    else:value=polynomial(T,HIGH if T>=1000 else LOW)
    if continuous and T>=1000:value[1:]-=JUMP[1:]
    return value


def from_pressure(p,p0,T0,R):
    _,h0,s0=air(T0);wanted=s0+R*np.log(p/p0)
    T=brentq(lambda t:air(t)[2]-wanted,80,T0,xtol=1e-11,rtol=1e-14)
    cp,h,s=air(T);U=np.sqrt(2*(h0-h))
    return T,p/(R*T),U,cp,h,s


def plane_reference(spec,arc):
    """用横向速度积分守恒构造独立ODE，不复用曲面压力转向方程。"""
    p0,T0,R=[float(spec[k]) for k in ['totalPressurePa','totalTemperatureK','R']]
    entry_p=float(spec['entryPressurePa']);k=float(spec['pressureExponentPerM'])
    C=float(spec['entrySpeedMS'])*np.sin(spec['entryAngleRad'])
    def rhs(_,state):
        U=from_pressure(entry_p*np.exp(k*state[0]),p0,T0,R)[2]
        return [-np.sqrt(1-(C/U)**2),C/U]
    result=solve_ivp(rhs,[0,float(arc[-1])],[0.,0.],method='DOP853',t_eval=arc,
                     rtol=2e-12,atol=1e-13,max_step=.01)
    assert result.success and result.y.shape==(2,len(arc))
    tangent=np.array([rhs(s,result.y[:,i]) for i,s in enumerate(arc)])
    return np.c_[result.y.T,np.zeros(len(arc))],np.c_[tangent,np.zeros(len(arc))],C


def cone_reference(source,cone):
    """独立法向跃迁根和Taylor–Maccoll积分，核对实际内部锥面源。"""
    R=float(source['R']);T1=float(source['temperature']);p1=float(source['pressure'])
    rho1=float(source['density']);U1=float(source['speed'])
    beta=np.deg2rad(float(cone['shockAngleDeg']));delta=np.deg2rad(float(cone['coneHalfAngleDeg']))
    un=U1*np.sin(beta);h1=air(T1,False)[1];H0=h1+U1*U1/2
    def jump_energy(ratio):
        p2=p1+rho1*un*un*(1-1/ratio);T2=p2/(rho1*ratio*R)
        return air(T2,False)[1]+(un/ratio)**2/2-h1-un*un/2
    grid=1+np.geomspace(1e-5,9.,120);bracket=None;last=None
    for ratio in grid:
        try:value=jump_energy(ratio)
        except AssertionError:continue
        if last is not None and value*last[1]<0:bracket=(last[0],ratio);break
        last=(ratio,value)
    assert bracket is not None
    ratio=brentq(jump_energy,*bracket,xtol=1e-13,rtol=1e-13)
    p2=p1+rho1*un*un*(1-1/ratio);T2=p2/(rho1*ratio*R);S2=air(T2,False)[2]
    def temperature(velocity):
        h=H0-np.dot(velocity,velocity)/2
        return brentq(lambda t:air(t,False)[1]-h,80,2500,xtol=1e-11,rtol=1e-14)
    def rhs(theta,velocity):
        radial,polar=velocity;T=temperature(velocity);cp=air(T,False)[0];a2=cp/(cp-R)*R*T
        return [polar,(polar*polar*radial-a2*(2*radial+polar/np.tan(theta)))/(a2-polar*polar)]
    result=solve_ivp(rhs,[beta,delta],[U1*np.cos(beta),-un/ratio],method='DOP853',
                     rtol=2e-11,atol=1e-9,max_step=.001)
    assert result.success
    velocity=result.y[:,-1];T=temperature(velocity);U=float(np.linalg.norm(velocity))
    p=p2*np.exp((air(T,False)[2]-S2)/R)
    assert abs(velocity[1])<1e-3
    expected=np.array([p,T,U]);actual=np.array([cone[k] for k in ['pressure','temperature','speed']])
    error=np.max(abs(actual/expected-1));assert error<5e-7
    return dict(pressurePa=p,temperatureK=T,speedMS=U,wallPolarSpeedMS=float(velocity[1]),
                sourceRelativeStateError=float(error),experimentalValidation=False)


def audit_thermal_ray(ray,result,spec,kind):
    arc=np.asarray(result['arcLengthM']);valid=np.asarray(ray['valid'],dtype=bool);assert np.all(valid)
    R=float(result['gas']['R']);refT=float(result['gas']['referenceTemperatureK']);refP=float(result['gas']['referencePressurePa'])
    coupling=result['thermalPressureCoupling'];p0=float(coupling['totalPressure']);T0=float(coupling['totalTemperature'])
    assert coupling['thermalPropertyConvention']=='continuous-project-cp-v1'
    assert not coupling['totalPressureResetAlongRay'] and not coupling['fullExternalContinuitySolved']
    H0=air(T0)[1];S0=air(T0)[2]-air(refT)[2]-R*np.log(p0/refP)
    for i,sample in enumerate(rows(ray['samples'])):
        p,T,rho,U=[float(sample[k]) for k in ['pressure','temperature','density','speed']]
        cp,h,S=air(T)
        np.testing.assert_allclose(p,rho*R*T,rtol=2e-14)
        np.testing.assert_allclose(h+U*U/2,H0,rtol=2e-12,atol=2e-8)
        np.testing.assert_allclose(S-air(refT)[2]-R*np.log(p/refP),S0,rtol=2e-12,atol=2e-8)
        np.testing.assert_allclose(ray['totalEnthalpy'][i],h+U*U/2,rtol=2e-12,atol=2e-8)
        np.testing.assert_allclose(ray['entropy'][i],S0,rtol=2e-12,atol=2e-8)
        dp=np.asarray(sample['pressureGradient']);dT=np.asarray(sample['temperatureGradient'])
        dr=np.asarray(sample['densityGradient']);dU=np.asarray(sample['speedGradient'])
        np.testing.assert_allclose(cp*dT+U*dU,0,rtol=0,atol=2e-8)
        np.testing.assert_allclose(cp*dT/T-R*dp/p,0,rtol=0,atol=2e-10)
        np.testing.assert_allclose(dp,R*(T*dr+rho*dT),rtol=2e-12,atol=2e-9)
        np.testing.assert_allclose(dp+rho*U*dU,0,rtol=0,atol=2e-8)
    assert np.max(ray['unitTangentError'])<=result['options']['maximumUnitTangentError']
    if kind=='NASA7_exponential_pressure_plane':
        position,tangent,C=plane_reference(spec,arc);position[:,1]-=float(ray['label'])
        pe=np.max(np.linalg.norm(ray['position']-position,axis=1));te=np.max(np.linalg.norm(ray['tangent']-tangent,axis=1))
        assert pe<2e-8 and te<2e-8
        momentum=float(np.max(abs(ray['speed']*ray['tangent'][:,1]/C-1)));assert momentum<2e-7
        return float(pe),float(te),momentum
    beta=float(ray['label']);angle=float(spec['coneAngleRad']);meridian=float(spec['startMeridianM'])+arc
    tangent=np.tile([-np.cos(angle),np.sin(angle)*np.cos(beta),np.sin(angle)*np.sin(beta)],(len(arc),1))
    position=meridian[:,None]*tangent
    pe=np.max(np.linalg.norm(ray['position']-position,axis=1));te=np.max(np.linalg.norm(ray['tangent']-tangent,axis=1))
    assert pe<2e-10 and te<2e-10
    return float(pe),float(te),0.


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');parser.add_argument('process')
    parser.add_argument('--report');parser.add_argument('--failure');args=parser.parse_args()
    folder=Path(args.folder).resolve();folder.relative_to(ROOT)
    target=Path(args.report).resolve() if args.report else folder/'independent_thermal_pressure_surface_flow.json'
    target.relative_to(ROOT);assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True)
    source=folder/'thermal_pressure_surface_flow_evidence.mat';data=loadmat(source,simplify_cells=True)
    process_path=Path(args.process).resolve();process=json.loads(process_path.read_text(encoding='utf-8-sig'))
    assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
    matching=data['report']['continuousThermodynamicMatching']
    np.testing.assert_allclose([matching['highEnthalpyOffsetJPerKg'],matching['highEntropyOffsetJPerKgK']],-JUMP[1:],rtol=0,atol=2e-8)
    checks=[];cone_check=None;maximum_momentum=0.
    for case in rows(data['cases']):
        result=case['result'];kind=case['kind'];spec=case['specification']
        assert np.all(result['validStation']) and not result['productionQualified'] and not result['newProductionSolvedAreaM2']
        if kind.startswith('caloric_'):
            checks.append(audit_case(result,spec));continue
        assert result['thermalPropertyMode']=='thermallyPerfect'
        errors=[audit_thermal_ray(ray,result,spec,kind) for ray in rows(result['rays'])]
        rays=rows(result['rays']);transverse=np.cross(result['center']['normal'],result['center']['tangent'])
        J=(rays[3]['position']-rays[1]['position'])/float(result['labelStep'])
        expected=np.sum(transverse*J,axis=1)
        np.testing.assert_allclose(result['metric'],expected,rtol=2e-12,atol=2e-12)
        if kind=='NASA7_exponential_pressure_plane':
            np.testing.assert_allclose(result['metric'],-result['center']['tangent'][:,0],rtol=2e-12,atol=2e-12)
        else:
            anchor=case['sourceCone'];converted=case['convertedTotalState'];R=float(result['gas']['R'])
            H=air(anchor['temperature'])[1]+anchor['speed']**2/2
            T0=brentq(lambda t:air(t)[1]-H,anchor['temperature'],2500,xtol=1e-11,rtol=1e-14)
            p0=anchor['pressure']*np.exp((air(T0)[2]-air(anchor['temperature'])[2])/R)
            np.testing.assert_allclose([converted['totalPressure'],converted['totalTemperature']],[p0,T0],rtol=2e-12)
            assert not converted['sourcePhysicalStateChanged'] and not converted['totalPressureResetAlongRay']
            for key in ['pressure','temperature','density','speed']:
                assert converted['sourceAnchor'][key]==anchor[key]
                np.testing.assert_allclose(result['center'][key],anchor[key],rtol=2e-12,atol=2e-10)
            if cone_check is None:cone_check=cone_reference(case['sourceFreestream'],anchor)
        maximum_momentum=max(maximum_momentum,max(e[2] for e in errors))
        checks.append(dict(kind=kind,rays=len(errors),maximumPositionErrorM=max(e[0] for e in errors),
                           maximumTangentError=max(e[1] for e in errors),maximumRelativeYMomentumError=max(e[2] for e in errors)))
    prefix=data['prefix'];valid=np.asarray(prefix['validStation'],dtype=bool);first=np.flatnonzero(~valid)[0]
    assert first>0 and np.all(~valid[first:])
    for key in ['position','temperature','density','speed']:assert np.all(np.isnan(prefix['center'][key][~valid]))
    assert np.all(np.isnan(prefix['metric'][~valid]))
    assert all(ray['failureIdentifier']=='SurfaceFlow:Entropy' for ray in rows(prefix['rays']))
    inputs={source.relative_to(ROOT).as_posix():sha(source),process_path.relative_to(ROOT).as_posix():sha(process_path)}
    sources=['setupSurfacePressureFlow.m','airThermalProperties.m','airRecoveryTemperature.m','selectAirThermalModel.m',
             'solveThermallyPerfectConicalShock.m','tests/testThermalPressureSurfaceFlow.m','tests/pressureSurfaceBenchmarkProblem.m']
    sources+=['regionalFlow/surface/'+name+'.m' for name in ['createSurfaceFlowGas','surfaceFlowAirThermalProperties',
              'surfaceFlowThermodynamicState','surfaceFlowTotalStateFromAnchor','bindIsentropicSurfacePressure',
              'marchPressureDrivenSurfaceFlow','pressureSurfaceFlowRHS','pressureSurfaceFlowOptions']]
    for rel in sources:
        assert sha(ROOT/rel)==process['sourceSHA256'][rel],rel;inputs[rel]=sha(ROOT/rel)
    original_air=subprocess.check_output(['git','show','f056d11805f2a366c1f0ddc806f8f8cb6d07c686:airThermalProperties.m'],cwd=ROOT)
    assert hashlib.sha256(original_air).hexdigest()==sha(ROOT/'airThermalProperties.m')
    old_error=None
    if args.failure:
        failure=Path(args.failure).resolve();failure.relative_to(ROOT);old=loadmat(failure,simplify_cells=True)
        old_error=float(old['momentumRelativeError']);assert old_error>2e-7 and maximum_momentum<2e-7
        inputs[failure.relative_to(ROOT).as_posix()]=sha(failure)
    report=dict(status='passed-independent-continuous-thermodynamics-plane-invariant-cone-and-anchor-audit',
                cases=len(checks),checks=checks,independentCone=cone_check,
                originalNASA7HighMinusLowCpHStdS=JUMP.tolist(),maximumRelativeYMomentumError=maximum_momentum,
                originalFailedRelativeYMomentumError=old_error,originalProductionAirLibraryUnchanged=True,
                totalPressureResetPrefixRejected=True,inputSHA256=inputs,scriptSHA256=sha(__file__),
                caloricReferenceHelperSHA256=sha(ROOT/'tests/auditPressureDrivenSurfaceFlow.py'),
                independentExperimentalValidation=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:report[k] for k in ['status','cases','independentCone','maximumRelativeYMomentumError','originalFailedRelativeYMomentumError']}))


if __name__=='__main__':main()
