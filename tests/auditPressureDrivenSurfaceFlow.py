"""独立解析核对曲面流线、双间距横向度量和失败前缀；不调用MATLAB求解器。"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def rows(value):
    return [value] if isinstance(value, dict) else list(value)


def reference(spec, arc, beta=0.):
    """直接在世界坐标写直线、圆、螺旋及球面大圆，不复用曲面联络公式。"""
    s = np.asarray(arc).reshape(-1); z = np.zeros_like(s); kind = spec['kind']
    if kind in ('plane', 'unsupported'):
        t = np.array([-np.cos(spec['directionAngleRad']), np.sin(spec['directionAngleRad']), 0.])
        b = np.cross([0., 0., 1.], t)
        p = np.outer(s, t)+beta*b; tangent = np.tile(t, (len(s), 1)); normal = np.tile([0., 0., 1.], (len(s), 1))
    elif kind == 'focus':
        t = np.array([-1., beta/spec['focusDistanceM'], 0.]); t /= np.linalg.norm(t)
        p = np.outer(s, t)+[0., -beta, 0.]; tangent = np.tile(t, (len(s), 1)); normal = np.tile([0., 0., 1.], (len(s), 1))
    elif kind == 'vortex':
        r = spec['vortexRadiusM']-beta; a = s/r
        p = r*np.c_[np.cos(a), np.sin(a), z]; tangent = np.c_[-np.sin(a), np.cos(a), z]
        normal = np.tile([0., 0., 1.], (len(s), 1))
    elif kind == 'cylinder':
        r = spec['radius']; psi = spec['directionAngleRad']; a = beta+s*np.sin(psi)/r
        p = np.c_[-s*np.cos(psi), r*np.cos(a), r*np.sin(a)]
        tangent = np.c_[-np.cos(psi)*np.ones_like(s), -np.sin(psi)*np.sin(a), np.sin(psi)*np.cos(a)]
        normal = np.c_[z, np.cos(a), np.sin(a)]
    elif kind == 'cone':
        a = spec['coneAngleRad']; distance = spec['startMeridianM']+s
        p = distance[:, None]*np.array([-np.cos(a), np.sin(a)*np.cos(beta), np.sin(a)*np.sin(beta)])
        tangent = np.tile([-np.cos(a), np.sin(a)*np.cos(beta), np.sin(a)*np.sin(beta)], (len(s), 1))
        normal = np.tile([np.sin(a), np.cos(a)*np.cos(beta), np.cos(a)*np.sin(beta)], (len(s), 1))
    elif kind == 'sphere_sine':
        a = spec['startPolarRad']+s/spec['radius']
        normal = np.c_[np.cos(a), np.sin(a)*np.cos(beta), np.sin(a)*np.sin(beta)]
        p = spec['radius']*normal
        tangent = np.c_[-np.sin(a), np.cos(a)*np.cos(beta), np.cos(a)*np.sin(beta)]
    elif kind == 'sphere_geodesic':
        phi = spec['startPolarRad']; psi = spec['directionAngleRad']; a = s/spec['radius']
        n0 = np.array([np.cos(phi), np.sin(phi), 0.])
        t0 = np.cos(psi)*np.array([-np.sin(phi), np.cos(phi), 0.])+np.sin(psi)*np.array([0., 0., 1.])
        p = spec['radius']*(np.outer(np.cos(a), n0)+np.outer(np.sin(a), t0))
        tangent = -np.outer(np.sin(a), n0)+np.outer(np.cos(a), t0)
        rotate_label = np.array([[1., 0., 0.], [0., np.cos(beta), -np.sin(beta)], [0., np.sin(beta), np.cos(beta)]])
        p = p@rotate_label.T; tangent = tangent@rotate_label.T; normal = p/spec['radius']
    else:
        raise AssertionError(kind)
    rotation = spec['rotation']; origin = spec['origin']
    return p@rotation.T+origin, tangent@rotation.T, normal@rotation.T


def infinite_metric(spec, arc):
    s = np.asarray(arc); kind = spec['kind']
    if kind in ('plane', 'vortex'): return np.ones_like(s)
    if kind == 'cylinder': return np.full_like(s, spec['radius']*np.cos(spec['directionAngleRad']))
    if kind == 'cone': return (spec['startMeridianM']+s)*np.sin(spec['coneAngleRad'])
    if kind == 'sphere_sine': return spec['radius']*np.sin(spec['startPolarRad']+s/spec['radius'])
    p, t, n = reference(spec, s); local=(p-spec['origin'])@spec['rotation']
    J = np.c_[np.zeros_like(s), -local[:, 2], local[:, 1]]@spec['rotation'].T
    return np.sum(np.cross(n, t)*J, axis=1)


def audit_case(result, spec):
    arc = np.atleast_1d(result['arcLengthM']); rays = rows(result['rays'])
    assert len(rays) == 5 and result['centralRayIndexOneBased'] == 3
    assert not result['physicalStateClipped'] and not result['productionQualified']
    assert not result['metricDifferenceIsRigorousBound'] and not result['massContinuityOfFullExternalFlowSolved']
    position_errors = []; tangent_errors = []; maximum_input_residual = 0.
    gamma, R = float(result['gas']['gamma']), float(result['gas']['R']); cp = gamma*R/(gamma-1)
    for ray in rays:
        valid = np.asarray(ray['valid'], dtype=bool); p, t, n = reference(spec, arc, float(ray['label']))
        pe = np.max(np.linalg.norm(ray['position'][valid]-p[valid], axis=1))
        te = np.max(np.linalg.norm(ray['tangent'][valid]-t[valid], axis=1))
        position_errors.append(float(pe)); tangent_errors.append(float(te))
        assert pe < 2e-8 and te < 2e-8
        np.testing.assert_allclose(ray['normal'][valid], n[valid], rtol=0, atol=2e-8)
        assert np.all(ray['unitTangentError'][valid] <= result['options']['maximumUnitTangentError'])
        samples = rows(ray['samples'])
        for k in np.flatnonzero(valid):
            v = samples[k]; A = v['firstDerivatives']; state = ray['state'][k]
            actual_t = A@state[2:]; actual_t /= np.linalg.norm(actual_t)
            np.testing.assert_allclose(actual_t, ray['tangent'][k], rtol=0, atol=4e-15)
            assert v['isAttached'] and v['isSmooth']
            p_value, rho, temperature, speed = [float(v[name]) for name in ['pressure', 'density', 'temperature', 'speed']]
            assert 80 <= temperature <= 2000 and speed > 0
            np.testing.assert_allclose(p_value, rho*R*temperature, rtol=3e-15)
            np.testing.assert_allclose(cp*temperature+speed**2/2, cp*spec['totalTemperatureK'], rtol=3e-15)
            expected_p=spec['totalPressurePa']*(temperature/spec['totalTemperatureK'])**(gamma/(gamma-1))
            np.testing.assert_allclose(p_value,expected_p,rtol=3e-15)
            gradients={key:A@np.linalg.solve(A.T@A,np.asarray(v[key]).reshape(2)) for key in
                       ['pressureGradient','densityGradient','temperatureGradient','speedGradient']}
            gp, gr, gt, gu = [gradients[key] for key in ['pressureGradient','densityGradient','temperatureGradient','speedGradient']]
            scaled = np.linalg.norm(gp-R*(temperature*gr+rho*gt))/max(1.,np.linalg.norm(gp),R*temperature*np.linalg.norm(gr),R*rho*np.linalg.norm(gt))
            maximum_input_residual=max(maximum_input_residual,float(scaled));assert scaled < 2e-12
            assert abs(actual_t@(rho*speed*gu+gp)) <= 1e-8*max(1.,np.linalg.norm(gp),rho*speed*np.linalg.norm(gu))
            assert np.linalg.norm(cp*gt+speed*gu) <= 1e-8*max(1.,cp*np.linalg.norm(gt),speed*np.linalg.norm(gu))
            curvature=-(gp-actual_t*(actual_t@gp))/(rho*speed**2)
            if spec['kind']=='vortex':
                radial=(ray['position'][k]-spec['origin']);expected_curvature=-radial/np.dot(radial,radial)
            else: expected_curvature=np.zeros(3)
            assert np.linalg.norm(curvature-expected_curvature)<2e-8
            for name in ['pressure','temperature','density','speed']:
                assert float(ray[name][k])==float(v[name])
        extra=int(ray['RHSCalls'])-4*int(ray['substeps'])-int(np.count_nonzero(valid))
        assert (extra==0 if np.all(valid) else 1<=extra<=4)
    valid=np.asarray(result['validStation'],dtype=bool);h=float(result['labelStep']);center=float(result['labelCenter'])
    pp,_,_=reference(spec,arc,center+h/2);pm,_,_=reference(spec,arc,center-h/2)
    coarse_p,_,_=reference(spec,arc,center+h);coarse_m,_,_=reference(spec,arc,center-h)
    _,t,n=reference(spec,arc,center);b=np.cross(n,t);J=(pp-pm)/h;Jc=(coarse_p-coarse_m)/(2*h)
    finite=np.sum(b*J,axis=1);coarse=np.sum(b*Jc,axis=1);difference=abs(coarse-finite)/3
    np.testing.assert_allclose(np.asarray(result['metric'])[valid],finite[valid],rtol=2e-8,atol=2e-9)
    np.testing.assert_allclose(np.asarray(result['fullLabelSeparationNorm'])[valid],np.linalg.norm(J[valid],axis=1),rtol=2e-8,atol=2e-9)
    np.testing.assert_allclose(np.asarray(result['metricDifferenceIndicator'])[valid],difference[valid],rtol=2e-4,atol=3e-10)
    assert np.all(np.isnan(np.asarray(result['metric'])[~valid]))
    for key in ['coordinate','position','tangent','normal','pressure','temperature','density','speed','totalEnthalpy','entropy']:
        assert np.all(np.isnan(result['center'][key][~valid]))
    assert result['entryValidationRHSCalls']==5
    return dict(kind=spec['kind'],rayCount=5,maximumPositionErrorM=max(position_errors),
                maximumTangentError=max(tangent_errors),maximumEOSGradientResidual=maximum_input_residual,
                finiteMetricCheckedStations=int(np.count_nonzero(valid)),pressureIsInput=True,
                fullExternalContinuityOrHeatingValidated=False)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');parser.add_argument('process')
    parser.add_argument('--report',help='将复核结果写入新的项目内路径，避免改写冻结输入或历史回执。')
    parser.add_argument('--source-ref',help='显式核对历史Git提交中的MATLAB原字节；不会将旧输入冒充当前代码重算。')
    args=parser.parse_args()
    folder=Path(args.folder).resolve();folder.relative_to(ROOT);source=folder/'pressure_surface_flow_evidence.mat'
    target=Path(args.report).resolve() if args.report else folder/'independent_pressure_surface_flow.json'
    target.relative_to(ROOT);assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True)
    process_path=Path(args.process).resolve();process=json.loads(process_path.read_text(encoding='utf-8-sig'))
    assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
    d=loadmat(source,simplify_cells=True);cases=rows(d['cases']);checks=[]
    for case in cases:
        checks.append(audit_case(case['result'],case['specification']))
        p,t,_=reference(case['specification'],case['result']['arcLengthM'])
        np.testing.assert_allclose(case['referencePosition'],p,rtol=0,atol=2e-15)
        np.testing.assert_allclose(case['referenceTangent'],t,rtol=0,atol=2e-15)
        np.testing.assert_allclose(case['referenceMetric'],infinite_metric(case['specification'],case['result']['arcLengthM']),rtol=2e-15,atol=2e-15)
    for item in rows(d['rk']): checks.append(audit_case(item['result'],cases[1]['specification']))
    for item in rows(d['label']): checks.append(audit_case(item['result'],cases[3]['specification']))
    mapped_spec=dict(cases[1]['specification']);a=.4
    mapped_spec.update(rotation=np.array([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]]),origin=np.array([-3.,1.,2.]))
    checks.append(audit_case(d['mapped'],mapped_spec));checks.append(audit_case(d['scaled'],cases[5]['specification']))
    prefix_spec=dict(cases[0]['specification']);prefix_spec.update(kind='unsupported',directionAngleRad=0.)
    checks.append(audit_case(d['prefix'],prefix_spec));focus_spec=dict(cases[0]['specification']);focus_spec['kind']='focus'
    checks.append(audit_case(d['focus'],focus_spec))
    prefix=d['prefix'];focus=d['focus'];s=focus['arcLengthM'];bad=np.flatnonzero(~focus['validStation'].astype(bool))[0]
    assert abs(s[bad]-.8)<2e-15 and np.all(~focus['validStation'][bad:].astype(bool))
    assert np.all(prefix['validStation'].astype(bool)==(prefix['arcLengthM']<=.44+1e-12))
    errors=np.array([item['positionErrorM'] for item in rows(d['rk'])]);ratios=errors[:-1]/errors[1:]
    label_errors=np.array([item['relativeMetricError'] for item in rows(d['label'])]);label_ratios=label_errors[:-1]/label_errors[1:]
    assert np.all((ratios>8)&(ratios<24)) and np.all((label_ratios>3.9)&(label_ratios<4.1))
    bindings={source.relative_to(ROOT).as_posix():sha(source),process_path.relative_to(ROOT).as_posix():sha(process_path)}
    source_reference=None;historical_modules={}
    if args.source_ref:
        assert re.fullmatch('[0-9a-fA-F]{7,40}',args.source_ref),'必须显式指定Git提交SHA。'
        source_reference=subprocess.check_output(['git','rev-parse','--verify',args.source_ref+'^{commit}'],cwd=ROOT,text=True).strip()
    for rel in ['setupSurfacePressureFlow.m','regionalFlow/surface/marchPressureDrivenSurfaceFlow.m',
                'regionalFlow/surface/pressureSurfaceFlowRHS.m','regionalFlow/surface/pressureSurfaceFlowOptions.m',
                'tests/testPressureDrivenSurfaceFlow.m','tests/pressureSurfaceBenchmarkProblem.m']:
        if source_reference:
            content=subprocess.check_output(['git','show',source_reference+':'+rel],cwd=ROOT)
            digest=hashlib.sha256(content).hexdigest();assert digest==process['sourceSHA256'][rel],rel
            historical_modules[rel]=digest
        else:
            assert sha(ROOT/rel)==process['sourceSHA256'][rel],rel;bindings[rel]=sha(ROOT/rel)
    report=dict(status='passed-independent-analytic-pressure-turn-geometry-dual-metric-and-prefix-audit',
                checks=checks,RKErrorRatios=ratios.tolist(),labelErrorRatios=label_ratios.tolist(),
                focusFirstUnresolvedArcM=float(s[bad]),metricDifferenceIsRigorousBound=False,
                inputSHA256=bindings,scriptSHA256=sha(__file__),independentExperimentalValidation=False,
                fullExternalFlowOrHeatingSolved=False,newActualNSEvaluations=0,newProductionSolvedAreaM2=0)
    if source_reference:
        report['historicalSourceIdentity']={'gitCommit':source_reference,'modules':historical_modules,
                                          'currentMATLABCodeRerunClaimed':False}
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(status=report['status'],auditedResults=len(checks),rays=5*len(checks),
                         RKErrorRatios=report['RKErrorRatios'],labelErrorRatios=report['labelErrorRatios'],focusM=float(s[bad]))))


if __name__=='__main__':main()
