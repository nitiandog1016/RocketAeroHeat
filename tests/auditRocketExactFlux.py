# -*- coding: utf-8 -*-
"""独立核对真实三维域精确通量、公开压力曲线及保存五量散度；不代替NS物理验证。"""
from pathlib import Path
import ast, hashlib, json, math, os, sys
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq, bisect

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v82_validation'
case = sys.argv[1] if len(sys.argv) > 1 else 'integrate01'
assert case.replace('_', '').isalnum()
target = VAL / (case + '_independent_audit.json')
assert not target.exists(), '不覆盖已有独立核对'
file = VAL / case / ('exact_integration_evidence.mat' if case == 'integrate01' else 'steady_krylov_evidence.mat')
domain_file = ROOT / 'outputs/v80_validation/domain02/shell_domain_evidence.mat'
source = VAL / 'clawpack_exact_euler.py'
probe = VAL / 'probe_exact_riemann.py'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings = {str(p.relative_to(ROOT)): sha(p) for p in [file, domain_file, source, probe]}
assert sha(source) == '2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'

# 模块一：只执行公开源码的纯数学压力曲线和先前独立标量Brent/二分参考。
allowed = {'pospart', 'sound_speed', 'beta', 'integral_curve_1', 'integral_curve_3', 'hugoniot_locus_1', 'hugoniot_locus_3'}
nodes = [n for n in ast.parse(source.read_text(encoding='utf8')).body if isinstance(n, ast.FunctionDef) and n.name in allowed]
assert {n.name for n in nodes} == allowed
for definition in nodes:
    for n in ast.walk(definition):
        assert not isinstance(n, (ast.Import, ast.ImportFrom, ast.Global, ast.Nonlocal, ast.With, ast.Try))
        if isinstance(n, ast.Call):
            assert (isinstance(n.func, ast.Name) and n.func.id in allowed) or (isinstance(n.func, ast.Attribute) and isinstance(n.func.value, ast.Name) and n.func.value.id == 'np' and n.func.attr in ['maximum', 'sqrt'])
reference = {'np': np}
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), reference)
exact_node = [n for n in ast.parse(probe.read_text(encoding='utf8')).body if isinstance(n, ast.FunctionDef) and n.name == 'exact']
assert len(exact_node) == 1
namespace = dict(np=np, math=math, brentq=brentq, bisect=bisect, reference=reference)
exec(compile(ast.Module(body=exact_node, type_ignores=[]), str(probe), 'exec'), namespace)
exact = namespace['exact']

d = loadmat(domain_file, simplify_cells=True, variable_names=['domain'])['domain']
s = loadmat(file, simplify_cells=True)
r = s.get('result', s)
q, info, gas = r['Q'] if 'Q' in r else s['Q'], r['info'], s['gas']
f, w = d['openFaces'], d['walls']
ids = np.flatnonzero((f['area'] > 0) & (f['cellMinus'] > 0) & (f['cellPlus'] > 0))
m, p = f['cellMinus'].astype(int) - 1, f['cellPlus'].astype(int) - 1
a = info['exactEulerFluxAudit']
assert np.array_equal(a['faceID'].astype(int) - 1, ids) and len(ids) == 849336
normal = f['normal'][ids]
gamma, R = gas['gamma'], gas['R']
cv = R / (gamma - 1)
rho = q[:, 0]
velocity = q[:, 1:4] / rho[:, None]
pressure = (gamma - 1) * (q[:, 4] - .5 * rho * np.sum(velocity ** 2, axis=1))
pl, pr = pressure[m[ids]], pressure[p[ids]]
rl, rr = rho[m[ids]], rho[p[ids]]
ul = np.sum(velocity[m[ids]] * normal, axis=1)
ur = np.sum(velocity[p[ids]] * normal, axis=1)
ps, us = a['pressureStarPa'], a['normalVelocityStarMs']

# 模块二：全部849336面用第三方原压力曲线核对正根残差与接触速度。
assert ps.min() > 1e-15  # 公开参考pospart下限在本批真实正根没有触发。
phi_l = np.where(ps >= pl, reference['hugoniot_locus_1'](ps, rl, ul, pl, gamma), reference['integral_curve_1'](ps, rl, ul, pl, gamma))
phi_r = np.where(ps >= pr, reference['hugoniot_locus_3'](ps, rr, ur, pr, gamma), reference['integral_curve_3'](ps, rr, ur, pr, gamma))
velocity_scale = np.sqrt(gamma * pl / rl) + np.sqrt(gamma * pr / rr) + abs(ur - ul)
curve_error = float(np.max(abs(phi_l - phi_r) / velocity_scale))
contact_error = float(np.max(abs(.5 * (phi_l + phi_r) - us) / velocity_scale))
assert curve_error < 5e-12 and contact_error < 5e-12

# 模块三：固定抽样加每类波系与压力/温度极值，独立求根、扇采样及三维重建。
sample_primitive = a['samplePrimitive']
regions = np.asarray(a['sampleRegion'], dtype=int)
selected = list(np.random.default_rng(8201).choice(len(ids), 256, replace=False))
selected += [int(np.flatnonzero(regions == k)[0]) for k in np.unique(regions)]
for values in [ps, a['samplePressurePa'], sample_primitive[:, 4], a['vacuumVelocityMarginMs']]:
    selected += [int(np.argmin(values)), int(np.argmax(values))]
selected += [int(np.flatnonzero(ids == 70994)[0])]
selected = sorted(set(selected))
errors = dict(pressureStar=0., normalVelocityStar=0., samplePrimitive=0., flux=0.)
for j in selected:
    ref = exact((rl[j], ul[j], pl[j]), (rr[j], ur[j], pr[j]), gamma)
    rs, un, pressure_sample = ref['sample']
    side = m[ids[j]] if ref['side'] == 'left' else p[ids[j]]
    n = normal[j]
    tangent = velocity[side] - np.dot(velocity[side], n) * n
    v = tangent + un * n
    temperature = pressure_sample / (rs * R)
    primitive = np.r_[rs, v, temperature]
    energy = pressure_sample / (gamma - 1) + .5 * rs * np.dot(v, v)
    flux = np.r_[rs * un, rs * v * un + pressure_sample * n, (energy + pressure_sample) * un]
    actual_flux = info['convectiveFaceFlux'][ids[j]]
    errors['pressureStar'] = max(errors['pressureStar'], abs(ps[j] - ref['pressureStar']) / max(1, abs(ps[j])))
    errors['normalVelocityStar'] = max(errors['normalVelocityStar'], abs(us[j] - ref['velocityStar']) / max(1, abs(us[j])))
    errors['samplePrimitive'] = max(errors['samplePrimitive'], float(np.max(abs(primitive - sample_primitive[j])) / max(1, np.max(abs(primitive)))))
    errors['flux'] = max(errors['flux'], float(np.max(abs(flux - actual_flux)) / max(1, np.max(abs(flux)))))
    assert regions[j] == a['regionNames'].tolist().index(ref['sampleRegion']) + 1
assert max(errors.values()) < 2e-10

# 模块四：逐面五量散度和长精度边界账；这不是另一套完整NS求解器。
volume = d['cells']['volume']
inside = (m >= 0) & (p >= 0)
summed = np.zeros((len(volume), 5))
flux = f['area'][:, None] * info['faceFlux']
np.add.at(summed, m, -flux)
np.add.at(summed, p[inside], flux[inside])
owner = w['ownerCell'].astype(int) - 1
np.add.at(summed, owner, -w['area'][:, None] * info['wallFlux'])
rate = r['rate']
divergence_error = float(np.max(abs(summed / volume[:, None] - rate)) / max(1, np.max(abs(rate))))
assert divergence_error < 2e-12
exterior = ~inside
boundary = (f['area'][exterior, None] * info['faceFlux'][exterior]).astype(np.longdouble).sum(axis=0) + (w['area'][:, None] * info['wallFlux']).astype(np.longdouble).sum(axis=0)
magnitude = np.maximum(1, (f['area'][exterior, None] * abs(info['faceFlux'][exterior])).sum(axis=0) + (w['area'][:, None] * abs(info['wallFlux'])).sum(axis=0))
conservation = float(np.max(abs(summed.astype(np.longdouble).sum(axis=0) + boundary) / magnitude))
assert conservation < 1e-9
assert not a['automaticFallbackUsed'] and not a['temperatureClipped'] and not a['pressureClipped']

# 模块五：新方程初猜身份、实际接受更新、下降规则和方向导数算术分别核对。
iteration_audit = None
if 'result' in s:
    integration_file = VAL / 'integrate01/exact_integration_evidence.mat'
    bindings[str(integration_file.relative_to(ROOT))] = sha(integration_file)
    base = loadmat(integration_file, simplify_cells=True, variable_names=['Q', 'rate'])
    assert np.array_equal(r['initialQ'], base['Q']) and np.array_equal(r['initialRate'], base['rate'])
    del base
    scale = np.asarray(r['residualScale'])
    trace = r['trace']
    trace = [trace] if isinstance(trace, dict) else list(trace)
    assert len(trace) == r['acceptedNonlinearIterations']
    previous_q, previous_rate = r['initialQ'], r['initialRate']
    update_error, accepted_t = 0., []
    for item in trace:
        assert np.array_equal(item['initialQ'], previous_q) and np.array_equal(item['initialRate'], previous_rate)
        if r['options'].get('trialStateCoordinates', 'conservative') == 'primitive-tangent':
            oldq = item['initialQ']; delta = item['linearUpdateU'] * scale
            oldrho = oldq[:, 0]; oldu = oldq[:, 1:4] / oldrho[:, None]
            oldt = (oldq[:, 4] - .5 * oldrho * np.sum(oldu ** 2, axis=1)) / (oldrho * cv)
            drho = delta[:, 0]; du = (delta[:, 1:4] - oldu * drho[:, None]) / oldrho[:, None]
            dt = (delta[:, 4] - np.sum(oldu * delta[:, 1:4], axis=1) + (.5 * np.sum(oldu ** 2, axis=1) - cv * oldt) * drho) / (oldrho * cv)
            nr = oldrho + item['alpha'] * drho; nu = oldu + item['alpha'] * du; nt = oldt + item['alpha'] * dt
            expected = np.column_stack((nr, nr[:, None] * nu, nr * (cv * nt + .5 * np.sum(nu ** 2, axis=1))))
            assert item['trialStateCoordinates'] == 'primitive-tangent'
        else:
            expected = item['initialQ'] + item['alpha'] * item['linearUpdateU'] * scale
        update_error = max(update_error, float(np.max(abs(expected - item['Q'])) / max(1, np.max(abs(expected)))))
        assert update_error < 3e-15
        iq = item['Q']
        irho = iq[:, 0]
        iv = iq[:, 1:4] / irho[:, None]
        it = (iq[:, 4] - .5 * irho * np.sum(iv ** 2, axis=1)) / (irho * cv)
        assert np.all(irho > 0) and np.all((it >= gas['temperatureRange'][0]) & (it <= gas['temperatureRange'][1]))
        accepted_t += [float(it.min()), float(it.max())]
        before = float(np.sqrt(np.mean((item['initialRate'] / scale) ** 2)))
        after = float(np.sqrt(np.mean((item['rate'] / scale) ** 2)))
        assert after <= (1 - r['options']['sufficientDecrease'] * item['alpha']) * before
        assert item['linearReport']['fullSystemTargetMet'] and item['linearReport']['fullSystemRelativeDefect'] <= r['options']['linearRelativeTolerance']
        previous_q, previous_rate = item['Q'], item['rate']
    assert np.array_equal(previous_q, r['Q']) and np.array_equal(previous_rate, r['rate'])
    jvp = r['firstJvpEvidence']
    expected = (jvp['ratePlus'] - jvp['rateMinusOrBase']) / scale / jvp['actualScalarStep']
    jvp_error = float(np.max(abs(expected - jvp['value'])) / max(1, np.max(abs(expected))))
    assert jvp_error < 3e-15 and np.array_equal(jvp['baseQ'], r['initialQ'])
    logs = r['jvpLog']; logs = [logs] if isinstance(logs, dict) else list(logs)
    assert len(logs) == r['matrixVectorProducts']
    assert sum(x['actualOriginalResidualEvaluations'] for x in logs) == r['jacobianResidualEvaluations']
    assert not r['physicalSolution'] and r['physicalTimeAdvanced'] == 0 and r['BETimeStepsClaimed'] == 0
    assert not s['mapped']['eligibleForSurfaceCoverage'] and not s['lineage']['sourceHistoryAppended']
    maximum = float(np.max(abs(r['rate'] / scale)))
    rms = float(np.sqrt(np.mean((r['rate'] / scale) ** 2)))
    assert abs(maximum - r['maximumResidual']) / maximum < 3e-15 and abs(rms - r['rmsMerit']) / rms < 3e-15
    if r['steadyConverged']:
        assert maximum <= 2e-6 and r['steadyBoundaryFluxBalance'] <= 5e-10
    iteration_audit = dict(acceptedIterations=len(trace), maximumResidual=maximum, rmsMerit=rms,
                           updateArithmeticRelativeDifference=update_error, firstFullResidualJvpArithmeticDifference=jvp_error,
                           acceptedTemperatureRangeK=[min(accepted_t), max(accepted_t)] if accepted_t else None,
                           actualFullResidualJvpCalls=int(r['jacobianResidualEvaluations']),
                           sourceHistoryAppended=False, initialNewEquationQAndRateExactlyReplayed=True,
                           nonlinearTrialStateCoordinates=r['options'].get('trialStateCoordinates', 'conservative'),
                           steadyConverged=bool(r['steadyConverged']))
assert all(sha(ROOT / p) == h for p, h in bindings.items())
report = dict(status='passed-independent-public-pressure-curves-all-real-faces-and-scalar-sampling', case=case,
              processID=os.getpid(), interiorFaces=len(ids), independentlyRootedSamples=len(selected),
              publicCurveRelativeResidual=curve_error, publicContactVelocityRelativeDifference=contact_error,
              sampledMaximumRelativeDifferences=errors, sampleRegionsPresent=np.unique(regions).tolist(),
              savedFiveQuantityDivergenceRelativeDifference=divergence_error, independentSpatialConservationRelative=conservation,
              referencePressureFloorActive=False, fullIndependentNSSolver=False, physicalSolution=False,
              newProductionSolvedAreaM2=0, iterationAudit=iteration_audit, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
print(json.dumps(report))
