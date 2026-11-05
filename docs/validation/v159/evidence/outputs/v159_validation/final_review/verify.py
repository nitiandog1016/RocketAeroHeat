"""v159只读交付复核，原始计算回执和失败报告保持准确来源。"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import types

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts'))
from su2BenchmarkAcceptance import acceptance
from su2IdealAirCheckpoints import output_plan
from su2ThermalRuntime import read, sha

OUT = Path(__file__).resolve().parent
CASE = ROOT/'outputs/v159_validation/public01'
review = read(CASE/'public_review02/report.json')
receipt = read(CASE/'engine_process.json')

# 模块一。逐字核对13对真实场和输入，源码后续修订使用实际执行副本重绑定。
complete = receipt['checkpoints']['complete']
assert [row['localIteration'] for row in complete] == list(range(0, 1200, 100))+[1199]
for row in complete:
    for kind in ('volume', 'surface'):
        assert sha(CASE/row[kind]) == row[kind+'SHA256'] == receipt['outputSHA256'][row[kind]]
snapshots = [p for p in CASE.rglob('*') if p.is_file() and
             ('executed_sources' in p.parts or p.name == 'executed_audit.py')]
snapshots += list(OUT.glob('before_*.py'))
rebound = []
for name, expected in dict(receipt['inputSHA256'], **review['sourceSHA256']).items():
    path = Path(name)
    if sha(path) != expected:
        matching = [p for p in snapshots if sha(p) == expected]
        assert matching, name
        rebound.append(dict(original=name, exactSnapshot=str(matching[0]), sha256=expected))

# 模块二。直接从原CSV提取末态和历史，与审核报告独立比对。
with (CASE/complete[-1]['volume']).open(encoding='utf8') as stream:
    names = next(csv.reader(stream))
data = np.loadtxt(str(CASE/complete[-1]['volume']), delimiter=',', skiprows=1)
f = {name: data[:, index] for index, name in enumerate(names)}
c = review['currentFieldAssessment']
ids = np.array(c['wall']['ids'], dtype=int)
errors = {}
for label, field in [('p', 'Pressure'), ('q', 'Heat_Flux'), ('fluidT', 'Temperature'), ('yplus', 'Y_Plus')]:
    expected = np.array(c['wall'][label])
    errors['wall_'+label] = float(max(abs(f[field][ids]-expected)/np.maximum(abs(expected), 1)))
    assert errors['wall_'+label] < 2e-12
with (CASE/'history.csv').open(encoding='utf8') as stream:
    hist = list(csv.reader(stream, skipinitialspace=True))
head = [name.strip() for name in hist[0]]
last = list(map(float, hist[-1]))
assert len(hist)-1 == review['completedUpdates'] == 1200
assert last[head.index('Inner_Iter')] == 1199
for field, label in [('Residual_Density', 'rms[Rho]'), ('Residual_Energy', 'rms[RhoE]')]:
    value = float(np.log10(np.linalg.norm(f[field])/np.sqrt(len(data))))
    errors[label] = abs(value-last[head.index(label)])
    assert errors[label] < 1e-8

# 模块三。用两点高斯积分独立复算整周壁热功率、测站误差和21点速度剖面。
radius = np.array(c['wall']['radius'])
x = np.array(c['wall']['x'])
q = np.array(c['wall']['q'])
length = np.hypot(np.diff(radius), np.diff(x))
power = 0.
for t in ((1-1/np.sqrt(3))/2, (1+1/np.sqrt(3))/2):
    power += np.pi*np.dot(length, ((1-t)*radius[:-1]+t*radius[1:])*((1-t)*q[:-1]+t*q[1:]))
errors['revolvedPower'] = abs(power/c['revolvedHeatPowerW']-1)
assert errors['revolvedPower'] < 2e-12
contract = read(ROOT/'tests/data/tmr_mach7/contract.json')
for group, rows, field, key in (
        ('pressure', [r for r in contract['measurements']['pressure_rows'] if r['independent_comparison_point']], 'p', 'pressure_Pa'),
        ('heat', contract['measurements']['heating_rows'], 'q', 'reported_heating_magnitude_W_m2')):
    predicted = np.interp([r['S_m'] for r in rows], c['wall']['S'], c['wall'][field])
    measured = np.array([r[key] for r in rows])
    rms = float(np.linalg.norm((predicted-measured)/measured)/np.sqrt(len(rows)))
    errors[group+'RMS'] = abs(rms-c['stations'][group]['rmsRelativeError'])
    assert errors[group+'RMS'] < 2e-12
with (ROOT/'outputs/v158_validation/preparation01/profile_locations.csv').open(encoding='utf8') as stream:
    locations = list(csv.DictReader(stream))
sampled = []
for row in locations:
    sampled.append(sum(float(row['weight%d' % i])*(-f['Momentum_x'][int(row['node%d' % i])]/
        f['Density'][int(row['node%d' % i])]) for i in range(4))/1274.)
errors['upstreamVelocity'] = float(max(abs(np.array(sampled)-c['profile']['predictedU'])))
assert errors['upstreamVelocity'] < 2e-12
for quantity, unit in [('mass', 'KgS'), ('energy', 'W')]:
    vals = [r[quantity+unit] for r in c['boundaryBalance']['markers'].values()]
    value = abs(sum(r['net'] for r in vals))/sum(r['incoming'] for r in vals)
    assert abs(value-c['boundaryBalance'][quantity+'RelativeImbalance']) < 2e-12

# 模块四。只改待验原因的准确说明，原26项数值、门槛和状态逐项相等。
old_accept = acceptance(c, review['finalRecordedResiduals'], review['temporalStability'], True, True)
assert old_accept == review['acceptance']
new_accept = acceptance(c, review['finalRecordedResiduals'], review['temporalStability'], True, True,
    missing_reasons={'upstream total temperature profile': 'experimental total-temperature denominator unresolved'})
for old, new in zip(old_accept['gates'], new_accept['gates']):
    assert {k: v for k, v in old.items() if k != 'reason'} == {k: v for k, v in new.items() if k != 'reason'}
assert sum(row['passed'] is True for row in new_accept['gates']) == 4
assert len(new_accept['failed']) == 21 and len(new_accept['pending']) == 1
assert next(row for row in new_accept['gates'] if row['passed'] is None)['reason'] == 'experimental total-temperature denominator unresolved'

# 模块五。共享分析函数的旧默认调用与当前提交基线逐值比较。
old_text = subprocess.check_output(['git', 'show', 'HEAD:scripts/su2AxisymmetricAnalysis.py'], cwd=str(ROOT)).decode('utf8')
old_module = types.ModuleType('v158_analysis')
exec(compile(old_text, 'v158_analysis', 'exec'), old_module.__dict__)
import su2AxisymmetricAnalysis as current_analysis
small = ROOT/'outputs/v159_validation/small01'
points, cells, markers = old_module.mesh_data(small/'mesh.su2')
with (small/'restart_flow_000003.csv').open(encoding='utf8') as stream:
    columns = next(csv.reader(stream))
raw = np.loadtxt(str(small/'restart_flow_000003.csv'), delimiter=',', skiprows=1)
fields = {name: raw[:, index] for index, name in enumerate(columns)}
fields['velocity'] = np.column_stack((fields['Momentum_x'], fields['Momentum_y']))/fields['Density'][:, None]
conditions = dict(gamma=1.4, gasConstant=287.05, speedMS=1295., densityKgM3=.02257,
                  temperatureK=82.47, turbulenceIntensity=.00021)
for name in ('boundary_balance', 'residual_scales'):
    args = (points, cells, markers, fields, conditions) if name == 'boundary_balance' else (points, cells, fields, conditions)
    assert getattr(old_module, name)(*args) == getattr(current_analysis, name)(*args)
next_capacity = output_plan(335726, 9878, 19752, 150, 75, 1536*1024**2)
result = dict(passed=True, completedUpdates=1200, completePairs=13, errors=errors,
    sourceRebinding=rebound, originalReportSHA256=sha(CASE/'public_review02/report.json'),
    acceptanceWithSpecificPendingReason=new_accept, legacyDefaultAnalysisFunctions=2,
    nextRealDomainCapacity=next_capacity, physicsQualified=False, newQualifiedArea=0,
    sourceSHA256={str(ROOT/'scripts'/name): sha(ROOT/'scripts'/name) for name in
        ['auditSU2AxisymmetricAir.py', 'su2BenchmarkAcceptance.py', 'su2AxisymmetricAnalysis.py']})
(OUT/'report.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(passed=True, errors=errors, completePairs=13, sourceRebinding=len(rebound),
    nextCapacityBytes=next_capacity['estimatedMaximumBytes']), indent=2))
