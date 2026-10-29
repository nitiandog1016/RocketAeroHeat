# -*- coding: utf-8 -*-
"""独立核对逐标签H/S、轨迹、物理横向梯度、Crocco关系与真实斜激波入口。"""
from pathlib import Path
import argparse
import json
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from scipy.io import loadmat
from scipy.optimize import brentq

from auditThermalPressureSurfaceFlow import air, sha
from auditPressureDrivenSurfaceFlow import rows

ROOT = Path(__file__).resolve().parents[1]


def thermal(T, gas):
    if gas['thermalPropertyMode'] == 'caloricallyPerfect':
        cp = gas['gamma'] * gas['R'] / (gas['gamma'] - 1)
        return cp, cp*T, cp*np.log(T)
    return air(float(T))


def recover(p, total, gas):
    """独立Brent熵根，不使用MATLAB二分结果作为温度参考。"""
    R = gas['R']; T0 = total['totalTemperature']; p0 = total['totalPressure']
    cp0, H, S = thermal(T0, gas)
    if gas['thermalPropertyMode'] == 'caloricallyPerfect':
        T = T0*(p/p0)**((gas['gamma']-1)/gas['gamma'])
    else:
        wanted = S+R*np.log(p/p0)
        T = brentq(lambda t: thermal(t, gas)[2]-wanted, 80, T0, xtol=1e-11, rtol=1e-14)
    return T, p/(R*T), np.sqrt(2*(H-thermal(T, gas)[1]))


def oblique(upstream, theta_deg):
    """独立密度比能量根和弱支激波角根；不读取保存的激波角作为输入。"""
    R = upstream['R']; T1 = upstream['temperature']; p1 = upstream['pressure']; U = upstream['speed']
    rho1 = p1/(R*T1); H1 = air(T1, False)[1]

    def jump(beta):
        un = U*np.sin(beta)

        def energy(r):
            p2 = p1+rho1*un*un*(1-1/r); T2 = p2/(rho1*r*R)
            return air(T2, False)[1]+(un/r)**2/2-H1-un*un/2

        previous = None; bracket = None
        for r in 1+np.geomspace(1e-5, 9., 100):
            value = energy(r)
            if previous is not None and value*previous[1] < 0:
                bracket = (previous[0], r); break
            previous = (r, value)
        assert bracket is not None
        ratio = brentq(energy, *bracket, xtol=1e-13, rtol=1e-13)
        p2 = p1+rho1*un*un*(1-1/ratio); T2 = p2/(rho1*ratio*R)
        turn = beta-np.arctan(np.tan(beta)/ratio)
        speed = np.hypot(U*np.cos(beta), un/ratio)
        return turn, p2, T2, rho1*ratio, speed

    cp1 = air(T1, False)[0]; gamma1 = cp1/(cp1-R)
    mach_angle = np.arcsin(max(np.sqrt(gamma1*R*T1), np.sqrt(upstream['gamma']*R*T1))/U)
    beta = brentq(lambda b: jump(b)[0]-np.deg2rad(theta_deg), mach_angle+.08, 1.1, xtol=1e-13, rtol=1e-13)
    _, p, T, rho, speed = jump(beta)
    return np.array([p, T, rho, speed])


def reference_ray(spec, label, arc, total, gas):
    rotation = np.asarray(spec['rotation'])
    if spec['kind'] == 'cylinderShock':
        position = np.column_stack([-(spec['startAxialM']+spec['entrySkew']*label+arc),
                                    np.full(len(arc), spec['radiusM']*np.cos(label)),
                                    np.full(len(arc), spec['radiusM']*np.sin(label))])
        tangent = np.tile([-1., 0., 0.], (len(arc), 1))
    else:
        angle = spec['directionAngleRad']; t = np.array([-np.cos(angle), np.sin(angle)])
        b = np.array([-np.sin(angle), -np.cos(angle)])
        initial = np.r_[label*(b+spec['entrySkew']*t), t]
        if spec['pressureExponentPerM'] == 0:
            position = np.column_stack([initial[:2]+arc[:, None]*t, np.zeros(len(arc))])
            tangent = np.tile(np.r_[t, 0.], (len(arc), 1))
        else:
            def rhs(_, y):
                p = spec['pressurePa']*np.exp(spec['pressureExponentPerM']*y[0])
                _, rho, U = recover(p, total, gas)
                gradient = np.array([spec['pressureExponentPerM']*p, 0.])
                turn = -(gradient-y[2:]*np.dot(y[2:], gradient))/(rho*U*U)
                return np.r_[y[2:], turn]
            solution = solve_ivp(rhs, [0, arc[-1]], initial, t_eval=arc, method='DOP853',
                                 rtol=2e-12, atol=2e-13, max_step=.01)
            assert solution.success
            position = np.column_stack([solution.y[:2].T, np.zeros(len(arc))])
            tangent = np.column_stack([solution.y[2:].T, np.zeros(len(arc))])
    return position@rotation.T, tangent@rotation.T


def audit_case(case):
    spec = case['specification']; result = case['result']; gas = result['gas']; R = gas['R']
    arc = np.asarray(result['arcLengthM']); rays = rows(result['rays']); entries = rows(result['entryTotalStates'])
    records = rows(case['entryRecords']); labels = np.asarray(result['labels']); step = result['labelStep']
    assert len(rays) == len(entries) == len(records) == 5 and np.all(result['validStation'])
    assert not result['auxiliaryRayDerivativesArePhysical'] and not result['entryTotalStateResetAlongRay']
    assert not any(result[k] for k in ['pressureSolved', 'fullExternalContinuitySolved', 'heatFluxSolved',
                                     'wallTemperatureSolved', 'productionQualified', 'newProductionSolvedAreaM2'])
    H = []; S = []; max_position = 0.; max_tangent = 0.; shock_errors = []
    for ray, total, record, label in zip(rays, entries, records, labels):
        assert ray['label'] == record['label'] == label and np.all(ray['valid'])
        assert '辅助' in ray['sampleDerivativeRole']
        if spec['kind'] == 'cylinderShock':
            expected = oblique(record['upstream'], record['deflectionDeg'])
            actual = np.array([record['shock'][k] for k in ['pressure', 'temperature', 'density', 'speed']])
            error = float(np.max(abs(actual/expected-1))); assert error < 2e-8; shock_errors.append(error)
            anchor = total['sourceAnchor']
            assert not total['sourcePhysicalStateChanged'] and not total['totalPressureResetAlongRay']
            assert all(anchor[k] == record['shock'][k] for k in ['pressure', 'temperature', 'density', 'speed'])
        else:
            np.testing.assert_allclose(total['totalTemperature'], spec['totalTemperatureK']+
                                       spec['temperatureSlopeK']*label+spec['temperatureCurvatureK']*label**2, rtol=0, atol=0)
            np.testing.assert_allclose(total['totalPressure'], spec['totalPressurePa']*
                                       np.exp(spec['logPressureSlope']*label+spec['logPressureCurvature']*label**2), rtol=3e-15)
        T0 = total['totalTemperature']; p0 = total['totalPressure']; H.append(thermal(T0, gas)[1])
        S.append(thermal(T0, gas)[2]-thermal(gas['referenceTemperatureK'], gas)[2]-R*np.log(p0/gas['referencePressurePa']))
        position, tangent = reference_ray(spec, label, arc, total, gas)
        max_position = max(max_position, float(np.max(np.linalg.norm(ray['position']-position, axis=1))))
        max_tangent = max(max_tangent, float(np.max(np.linalg.norm(ray['tangent']-tangent, axis=1))))
        assert max_position < 2e-8 and max_tangent < 2e-8
        np.testing.assert_allclose(ray['totalEnthalpy'], H[-1], rtol=2e-11, atol=2e-7)
        np.testing.assert_allclose(ray['entropy'], S[-1], rtol=2e-11, atol=2e-8)
        for i, sample in enumerate(rows(ray['samples'])):
            T, rho, U = recover(sample['pressure'], total, gas)
            np.testing.assert_allclose([sample['temperature'], sample['density'], sample['speed']], [T, rho, U], rtol=2e-11, atol=2e-8)
            assert '辅助' in sample['sampleDerivativeRole']
    H = np.asarray(H); S = np.asarray(S)
    np.testing.assert_allclose(result['entryTotalEnthalpy'], H, rtol=2e-11, atol=2e-7)
    np.testing.assert_allclose(result['entryEntropy'], S, rtol=2e-11, atol=2e-8)
    J = (rays[3]['position']-rays[1]['position'])/step
    center = result['center']; t = np.asarray(center['tangent']); b = np.cross(center['normal'], t)
    h = np.sum(J*b, axis=1); shift = np.sum(J*t, axis=1)
    np.testing.assert_array_equal(result['labelPositionDerivative'], J)
    np.testing.assert_allclose(result['metric'], h, rtol=3e-14, atol=3e-14)
    np.testing.assert_allclose(result['labelAlongStreamShift'], shift, rtol=3e-14, atol=3e-14)
    dH = (H[3]-H[1])/step; dS = (S[3]-S[1])/step
    resolution = result['entryDerivativeResolution']
    savedH = np.asarray(result['entryTotalEnthalpy']); savedS = np.asarray(result['entryEntropy'])
    exactDH = (savedH[3]-savedH[1])/step; exactDS = (savedS[3]-savedS[1])/step
    errH = abs((savedH[4]-savedH[0])/(2*step)-exactDH)/3
    errS = abs((savedS[4]-savedS[0])/(2*step)-exactDS)/3
    initial_cp = [thermal(rows(ray['samples'])[0]['temperature'], gas)[0] for ray in rays]
    scaleH = 4*result['options']['relativeStateTolerance']*max(1., np.max(abs(savedH)))/step
    scaleS = 4*result['options']['relativeStateTolerance']*max(R, np.max(abs(savedS)), max(initial_cp))/step
    np.testing.assert_allclose([resolution['enthalpyDifferenceIndicator'], resolution['entropyDifferenceIndicator'],
                               resolution['enthalpyConsistencyScale'], resolution['entropyConsistencyScale']],
                              [errH, errS, scaleH, scaleS], rtol=3e-12, atol=1e-12)
    assert errH <= resolution['relativeDifferenceLimit']*abs(exactDH)+scaleH
    assert errS <= resolution['relativeDifferenceLimit']*abs(exactDS)+scaleS
    assert bool(resolution['enthalpyVariationResolved']) == (abs(exactDH) > scaleH)
    assert bool(resolution['entropyVariationResolved']) == (abs(exactDS) > scaleS)
    assert not resolution['consistencyScalesAreRigorousErrorBounds']
    if spec['kind'] == 'cylinderShock':
        assert not resolution['enthalpyVariationResolved'] and resolution['entropyVariationResolved']
    g = result['physicalGradients']; crocco = []; eos = []; chain = []
    for i, sample in enumerate(rows(center['samples'])):
        assert '物理' in sample['sampleDerivativeRole'] or '曲面两向' in sample['sampleDerivativeRole']
        p = sample['pressure']; T = sample['temperature']; rho = sample['density']; U = sample['speed']; cp = thermal(T, gas)[0]
        rotation = np.asarray(spec['rotation'])
        gradP = rotation@np.array([spec['pressureExponentPerM']*p, 0., 0.]) if spec['kind'] != 'cylinderShock' else np.zeros(3)
        gradH = dH/h[i]*b[i]; gradS = dS/h[i]*b[i]
        gradT = T/cp*(gradS+R*gradP/p); gradU = (gradH-cp*gradT)/U; gradRho = rho*(gradP/p-gradT/T)
        for key, expected in [('totalEnthalpyGradient', gradH), ('entropyGradient', gradS),
                              ('temperatureGradient', gradT), ('speedGradient', gradU),
                              ('densityGradient', gradRho), ('pressureGradient', gradP)]:
            np.testing.assert_allclose(g[key][i], expected, rtol=2e-8, atol=.003 if key == 'totalEnthalpyGradient' else 3e-6)
        A = sample['firstDerivatives']
        for key, expected in [('temperatureGradient', gradT), ('speedGradient', gradU), ('densityGradient', gradRho)]:
            np.testing.assert_allclose(sample[key], A.T@expected, rtol=2e-8, atol=3e-6)
        eos.append(float(np.linalg.norm(gradP-R*(T*np.asarray(g['densityGradient'][i])+rho*np.asarray(g['temperatureGradient'][i])))))
        velocity_beta = (rays[3]['speed'][i]*rays[3]['tangent'][i]-rays[1]['speed'][i]*rays[1]['tangent'][i])/step
        dUds = -np.dot(t[i], gradP)/(rho*U); curvature_b = -np.dot(b[i], gradP)/(rho*U*U)
        normal_vorticity = U*curvature_b-(np.dot(t[i], velocity_beta)-shift[i]*dUds)/h[i]
        expected_crocco = (dH-T*dS)/h[i]
        crocco.append(abs(-U*normal_vorticity-expected_crocco)/max(1., abs(dH/h[i]), abs(T*dS/h[i])))
        chain.append([(rays[3][k][i]-rays[1][k][i])/step-np.dot(g[v][i], J[i])
                      for k, v in [('temperature', 'temperatureGradient'), ('speed', 'speedGradient'), ('density', 'densityGradient')]])
    np.testing.assert_allclose(result['thermalLabelChainRuleDifference'], chain, rtol=2e-7, atol=2e-8)
    assert max(crocco) < 5e-5
    assert max(eos) < 2e-7
    return dict(kind=case['kind'], maximumPositionErrorM=max_position, maximumTangentError=max_tangent,
                maximumEOSGradientResidualPaPerM=max(eos), maximumRelativeCroccoDifference=max(crocco),
                maximumThermalLabelChainRuleDifference=np.max(np.abs(chain), axis=0).tolist(),
                independentShockMaximumRelativeError=max(shock_errors) if shock_errors else None,
                inletEnthalpyVariationResolved=bool(resolution['enthalpyVariationResolved']),
                inletEntropyVariationResolved=bool(resolution['entropyVariationResolved']))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('folder'); parser.add_argument('process'); parser.add_argument('--report')
    args = parser.parse_args(); folder = Path(args.folder).resolve(); folder.relative_to(ROOT)
    source = folder/'labelled_pressure_surface_flow_evidence.mat'; data = loadmat(source, simplify_cells=True)
    target = Path(args.report).resolve() if args.report else folder/'independent_labelled_surface_flow.json'
    target.relative_to(ROOT); assert not target.exists()
    process_path = Path(args.process).resolve(); process = json.loads(process_path.read_text(encoding='utf-8-sig'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
    cases = rows(data['cases']); assert len(cases) == data['report']['cases'] == 8
    assert 'SurfaceFlow:EntryDerivativeResolution' in rows(data['report']['rejections'])
    checks = [audit_case(c) for c in cases]
    turn_errors = np.array([c['maximumThermalLabelChainRuleDifference'] for c in checks[2:5]])
    ratios = turn_errors[:-1]/turn_errors[1:]; assert np.all((ratios > 3.9) & (ratios < 4.1))
    prefix = data['prefix']; valid = np.asarray(prefix['validStation'], dtype=bool); first = np.flatnonzero(~valid)[0]
    assert first > 0 and np.all(~valid[first:]) and np.all(np.isnan(prefix['metric'][~valid]))
    for key in ['temperature', 'density', 'speed', 'totalEnthalpy', 'entropy']:
        assert np.all(np.isnan(prefix['center'][key][~valid]))
    for gradient in prefix['physicalGradients'].values(): assert np.all(np.isnan(gradient[~valid]))
    bindings = {source.relative_to(ROOT).as_posix(): sha(source), process_path.relative_to(ROOT).as_posix(): sha(process_path)}
    leaves = ['regionalFlow/surface/marchLabelledPressureSurfaceFlow.m', 'tests/labelledSurfaceBenchmarkProblem.m',
              'tests/testLabelledPressureSurfaceFlow.m', 'regionalFlow/surface/pressureSurfaceFlowRHS.m',
              'regionalFlow/surface/bindIsentropicSurfacePressure.m', 'regionalFlow/surface/surfaceFlowAirThermalProperties.m',
              'regionalFlow/surface/surfaceFlowThermodynamicState.m', 'regionalFlow/surface/surfaceFlowTotalStateFromAnchor.m',
              'solveThermallyPerfectObliqueShock.m', 'solveThermallyPerfectNormalShock.m']
    for rel in leaves:
        assert sha(ROOT/rel) == process['sourceSHA256'][rel], rel; bindings[rel] = sha(ROOT/rel)
    unchanged = ['regionalFlow/surface/marchPressureDrivenSurfaceFlow.m', 'regionalFlow/surface/pressureSurfaceFlowRHS.m',
                 'regionalFlow/surface/bindIsentropicSurfacePressure.m', 'airThermalProperties.m']
    import hashlib
    for rel in unchanged:
        original = subprocess.check_output(['git', 'show', '24150d06e8a3e659c2e83a2a37bc987d5f7fbd81:'+rel], cwd=ROOT)
        assert hashlib.sha256(original).hexdigest() == sha(ROOT/rel)
    report = dict(status='passed-independent-labelled-inlet-rays-physical-gradients-Crocco-and-shock-anchor-audit',
                  cases=len(checks), rays=len(checks)*5, checks=checks, labelChainRuleRefinementRatios=ratios.tolist(),
                  oldV95NumericalLeavesUnchanged=unchanged, stoppedPrefixRemainsNaN=True,
                  inputSHA256=bindings, scriptSHA256=sha(__file__), thermalReferenceHelperSHA256=sha(ROOT/'tests/auditThermalPressureSurfaceFlow.py'),
                  independentExperimentalValidation=False, newProductionSolvedAreaM2=0)
    target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf8'))
    print(json.dumps({'status': report['status'], 'cases': len(checks), 'refinementRatios': ratios.tolist(),
                      'maxRelativeCroccoDifference': max(c['maximumRelativeCroccoDifference'] for c in checks)}))


if __name__ == '__main__':
    main()
