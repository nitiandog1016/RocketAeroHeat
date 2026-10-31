"""独立检查真实箭体参考流长约束、公开热流公式、场量身份和导出值。"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def vector(value):
    return np.asarray(value).reshape(-1)


def close(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape
    assert np.array_equal(np.isfinite(actual), np.isfinite(expected))
    assert np.allclose(actual, expected, rtol=8e-12, atol=2e-12, equal_nan=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', default='outputs/v104_validation/actual01')
    parser.add_argument('--baseline', default='outputs/v104_validation/baseline02/current_case.mat')
    parser.add_argument('--out', default='docs/validation/v104/independent_audit.json')
    args = parser.parse_args()
    folder = Path(args.folder)
    loaded = loadmat(str(folder / 'current_case.mat'), simplify_cells=True)
    previous = loadmat(args.baseline, simplify_cells=True)['current']
    mesh, fs = loaded['current'], loaded['fs']
    h, edge = mesh['heating'], mesh['boundaryLayerEdge']
    k, p = h['kinneyAcreage'], mesh['surfaceStreamlines']
    lengths = k['length']
    n = len(mesh['F'])
    assert np.array_equal(mesh['V'], previous['V']) and np.array_equal(mesh['F'], previous['F'])
    triangles = mesh['V'][mesh['F'].astype(int) - 1]
    area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                   triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    close(area, mesh['localGeometry']['area'])
    close(mesh['pressure']['staticPressure'], previous['pressure']['staticPressure'])
    assert np.array_equal(mesh['pressure']['fieldStatus'], previous['pressure']['fieldStatus'])
    for key in ('temperature', 'pressure', 'density', 'dynamicViscosity', 'speed'):
        close(edge[key], previous['boundaryLayerEdge'][key])
    accepted, target = vector(k['acceptedFace']).astype(bool), vector(k['targetFace']).astype(bool)
    old_known = vector(previous['heating']['eligibleFace']).astype(bool)
    known = vector(h['eligibleFace']).astype(bool)
    assert not np.any(known & ~old_known)
    blocked = vector(previous['continuumRegime']['stateUnresolvedFace']).astype(bool)
    assert not np.any(known[blocked]) and np.all(np.isnan(vector(h['heatFlux'])[blocked]))

    # 模块一：验证每条流线的积分守恒及根边界。约束全秩时唯一确定全积分。
    rho, mu, speed = (vector(edge[x]) for x in ('density', 'dynamicViscosity', 'speed'))
    weight = rho * mu * speed
    flow = vector(edge['algebraicallyValid']).astype(bool)
    flow &= vector(edge['thermalStateValidFace']).astype(bool)
    flow &= ~vector(edge['thermalExpansionOutOfRangeFace']).astype(bool)
    for key in ('continuumUnresolvedFace', 'separationUnresolvedFace',
                'turbulentSeparationUnresolvedFace', 'turbulentNumericalUnresolvedFace'):
        flow &= ~vector(h[key]).astype(bool)
    flow &= np.isfinite(rho) & (rho > 0) & np.isfinite(mu) & (mu > 0)
    flow &= np.isfinite(speed) & (speed >= 0)
    parent = vector(p['predecessorFace']).astype(int) - 1
    distance, step = vector(p['distance']), vector(p['stepLength'])
    connected = vector(p['connectedFace']).astype(bool)
    integral = vector(lengths['integral'])
    flow &= connected
    # 独立从各根向下检查可达性，未知父节点不得因本地流场有效而恢复。
    children = [[] for _ in range(n)]
    roots = []
    for i in np.flatnonzero(connected):
        if parent[i] < 0:
            roots.append(i)
        else:
            children[parent[i]].append(i)
    queue, complete = list(roots), np.zeros(n, dtype=bool)
    position = 0
    while position < len(queue):
        i = queue[position]
        position += 1
        complete[i] = flow[i] and (parent[i] < 0 or complete[parent[i]])
        queue.extend(children[i])
    assert len(queue) == np.count_nonzero(connected)
    assert np.array_equal(np.isfinite(integral), complete)
    for i in np.flatnonzero(complete):
        if parent[i] < 0:
            close(integral[i], weight[i] * distance[i])
        else:
            j = parent[i]
            close(distance[i] - distance[j], step[i])
            close(integral[i] - integral[j], .5 * (weight[i] + weight[j]) * step[i])
    weighted = np.full(n, np.nan)
    positive = complete & (weight > 0)
    weighted[positive] = integral[positive] / weight[positive]
    close(lengths['weightedLength'], weighted)
    x = np.clip(-vector(mesh['aerothermalRegion']['incidenceAngleDeg']) / lengths['leewardBlendDeg'], 0, 1)
    blend = x * x * (3 - 2 * x)
    close(lengths['blendWeight'], blend)
    selected = np.full(n, np.nan)
    use = complete & (blend == 0)
    selected[use] = distance[use]
    use = np.isfinite(weighted) & (blend > 0)
    selected[use] = distance[use] + blend[use] * (weighted[use] - distance[use])
    close(lengths['selectedLength'], selected)

    # 模块二：按原文Cf/St和焓差另算所有选中面元；负值和零值不裁切。
    te, tw = vector(edge['temperature']), vector(h['wallTemperatureByFace'])
    expected_heat = np.full(n, np.nan)
    for family in ('Laminar', 'Turbulent'):
        mask = accepted & vector(h[family.lower() + 'Face']).astype(bool)
        recovery = vector(h['recoveryTemperature' + family])
        cp = vector(h['effectiveCp' + family])
        tstar = (te + tw) / 2 + .22 * (recovery - te)
        rho_star = vector(edge['pressure']) / (fs['R'] * tstar)
        mu_star = 1.458e-6 * tstar**1.5 / (tstar + 110.4)
        reynolds = rho_star * speed * selected / mu_star
        close(k['referenceTemperature' + family], tstar)
        close(k['referenceReynolds' + family], reynolds)
        if family == 'Laminar':
            cf = .664 * np.sqrt(k['manglerLaminar'] / reynolds[mask]) * te[mask] / tstar[mask]
        else:
            cf = .37 * np.log10(reynolds[mask] / k['manglerTurbulent'])**(-2.584) * te[mask] / tstar[mask]
        expected_heat[mask] = rho_star[mask] * speed[mask] * cf / (2 * h['Prandtl']**(2/3)) * cp[mask] * (recovery[mask] - tw[mask])
    close(vector(h['heatFlux'])[accepted], expected_heat[accepted])
    with (folder / 'faces.csv').open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == n
    for csv_key, values in (
        ('face', np.arange(1, n+1)), ('area_m2', area),
        ('heat_flux_W_m2', h['heatFlux']), ('weight_integral', integral),
        ('selected_length_m', selected), ('weighted_length_m', weighted),
        ('heat_method_id', h['faceMethodID']), ('reason_code', mesh['heatingCoverage']['reasonCode'])):
        close(np.array([float(r[csv_key]) for r in rows]), vector(values))

    # 模块三：网页只是展示当前结果；给定壁温与未知求解壁温逐一保护。
    raw = (folder / 'workbench/case-data.js').read_text(encoding='utf-8')
    payload = json.loads(raw.partition('=')[2].strip().rstrip(';'))
    fields = {f['key']: f for f in payload['fields']}
    assert all(v == 2 for v in fields['wallInput']['status'])
    assert all(v == 0 for v in fields['wallSolved']['status'])
    for key, values, valid in (
        ('heatFlux', vector(h['heatFlux']) / 1000, known),
        ('heatingReferenceLength', selected, accepted),
        ('baselineHeatFlux', vector(k['baselineHeatFlux']) / 1000, old_known)):
        f = fields[key]
        assert np.array_equal(np.asarray(f['status']) == 1, valid)
        actual = np.array([np.nan if v is None else v for v in f['values']])
        expected = np.where(valid, values, np.nan)
        close(actual, expected)
    for i in np.flatnonzero(accepted):
        q = fields['heatFlux']
        assert '20070017947' in q['methodNames'][q['methodIndex'][i]]
        tr = fields['recoveryTemperature']
        assert 'Kinney' not in tr['methodNames'][tr['methodIndex'][i]]
    checks = len(np.flatnonzero(complete)) + int(np.count_nonzero(accepted)) + 8*n + 3*n
    report = dict(status='passed', faceCount=n, acceptedSidewallFaces=int(np.count_nonzero(accepted)),
                  integralConstraintFaces=int(np.count_nonzero(complete)), checkedValues=checks,
                  defaultHeatAreaFraction=float(np.sum(area[old_known])/np.sum(area)),
                  selectedHeatAreaFraction=float(np.sum(area[known])/np.sum(area)),
                  newlyComputedHeatArea=float(np.sum(area[known & ~old_known])),
                  blockedStateFaces=int(np.count_nonzero(blocked)),
                  source='Kinney 2007 II.I printed pages 4-5; source formulas, graph constraints and snapshot identities',
                  experimentalValidation=False, completePhysicalSolution=False,
                  auditorSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False))


if __name__ == '__main__':
    main()
