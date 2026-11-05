"""公开轴对称基准的同阶段状态、壁热与旋转积分独立审核。"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json, audit_geometry, source_blocks
from su2SSTInterface import state, wall_geometry, mesh, PHASE_OUTPUT


def audit(folder):
    """保留所有壁节点，包括入口/出口端点，检查失败时保存失败原文。"""
    folder = Path(folder).resolve()
    receipt = json.loads((folder/'engine_process.json').read_text(encoding='utf-8'))
    manifest = json.loads((folder/'benchmark_manifest.json').read_text(encoding='utf-8'))
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
    assert not receipt['timedOut'] and all(receipt['inputsStable'].values())
    for path, expected in receipt['sourceSHA256'].items():
        assert sha(path) == expected
    for name, expected in receipt['outputSHA256'].items():
        assert sha(folder/name) == expected
    assert manifest['model'] == 'SST-Vm (V1994m,VORTICITY)'
    c = manifest['conditions']
    geometry = audit_geometry(folder/'mesh.su2')
    points, cells, wall = mesh(folder/'mesh.su2')
    blocks, source = source_blocks()
    expected = np.concatenate((blocks[0].reshape(-1, 2), blocks[1][:, 1:].reshape(-1, 2)))
    assert np.array_equal(points, expected)
    assert source['archiveSHA256'] == manifest['source']['archiveSHA256']
    a, columns, nd, rho, velocity, thermal = state(folder/'restart_flow.csv', True, True)
    assert nd == 2 and len(a) == len(points)
    assert np.max(np.abs(a[:, 1:3]-points)) < 2e-14
    t = a[:, columns.index('Temperature')]
    p = a[:, columns.index('Pressure')]
    assert min(t) > 0 and min(p) > 0
    eos = float(np.max(np.abs(p/(rho*c['gasConstant']*t)-1)))
    assert eos < 2e-12
    g = wall_geometry(folder/'mesh.su2')
    ids = g['wallIDs']
    wall_temperature_error = float(np.max(np.abs(t[ids]-c['wallTemperatureK'])))
    # 原启动门槛继续单列。弱能量边界使用给定 Tw，保存的流体节点温度另记。
    wall_temperature_gate = wall_temperature_error < 2e-8
    viscosity = 1.716e-5*(t[ids]/273.15)**1.5*(273.15+110.4)/(t[ids]+110.4)
    native_mu = a[ids, columns.index('Laminar_Viscosity')]
    mu_error = float(np.max(np.abs(viscosity/native_mu-1)))
    assert mu_error < 2e-12
    cp = c['gamma']*c['gasConstant']/(c['gamma']-1)
    q_gradient = viscosity*cp/.72*(t[g['normalNeighbor']]-c['wallTemperatureK'])/g['thermalDistance']
    q = a[ids, columns.index('Heat_Flux')]
    gradient_error = float(np.max(np.abs(q_gradient-q)/np.maximum(np.abs(q_gradient), 1)))
    assert gradient_error < 2e-6, '独立壁热梯度误差 %.17g' % gradient_error
    assert np.max(np.abs(q)) > 1
    assert np.max(a[ids, columns.index('Eddy_Viscosity')]/viscosity) < 1e-12
    with (folder/'history.csv').open() as f:
        rows = list(csv.reader(f, skipinitialspace=True))
    labels = [x.strip() for x in rows[0]]
    history = np.array(rows[1:], dtype=float)
    assert np.array_equal(history[:, labels.index('Inner_Iter')], [0, 1])
    native_heat = float(history[-1, labels.index('HF')])
    planar_dual_heat = float(np.dot(g['nodeArea'], q))
    planar_error = abs(planar_dual_heat-native_heat)/max(abs(native_heat), 1)
    assert planar_error < 2e-6
    # 对每条线段的线性 q 和线性 r 做精确积分，再绕轴旋转一周。
    segments = np.array([row[1:] for row in wall])
    node_q = np.zeros(len(points))
    node_q[ids] = q
    endpoints = points[segments]
    radii = endpoints[:, :, 1]
    values = node_q[segments]
    length = np.linalg.norm(endpoints[:, 1]-endpoints[:, 0], axis=1)
    weighted = (2*values[:, 0]*radii[:, 0] + values[:, 0]*radii[:, 1]
                + values[:, 1]*radii[:, 0] + 2*values[:, 1]*radii[:, 1])/6
    full_power = float(2*np.pi*np.dot(length, weighted))
    # 独立的两点高斯积分对二次 q*r 精确，检查单位与径向权重。
    gauss_power = 0.
    for s in ((1-1/np.sqrt(3))/2, (1+1/np.sqrt(3))/2):
        rq = ((1-s)*radii[:, 0]+s*radii[:, 1])*((1-s)*values[:, 0]+s*values[:, 1])
        gauss_power += float(np.pi*np.dot(length, rq))
    power_error = abs(full_power-gauss_power)/max(abs(full_power), 1)
    assert power_error < 2e-12
    # 坐标和算法列随研究结果输出，便于后续在原测站进行定量比较。
    order = np.argsort(-points[ids, 0])
    with (folder/'wall_quantities.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['PointID', 'source_x_m', 'project_X_m', 'radius_m', 'pressure_Pa',
                    'q_into_wall_W_m2', 'given_Tw_K', 'y_plus', 'method', 'status'])
        for j in order:
            i = ids[j]
            w.writerow([i, -points[i, 0], points[i, 0], points[i, 1], p[i], q[j],
                        c['wallTemperatureK'], a[i, columns.index('Y_Plus')],
                        manifest['model'], 'startup_two_updates_unconverged'])
    report = dict(interfacePassed=wall_temperature_gate, thermalReadbackPassed=True,
                  originalStartupGates=dict(savedWallTemperatureMatchesGiven=wall_temperature_gate,
                                            positiveState=True, sameStageEOS=True,
                                            nonzeroHeat=True, heatGradient=True),
                  physicsQualified=False, newQualifiedArea=0,
                  model=manifest['model'], geometry=geometry, checkpointLocalIteration=1,
                  thermalStatePhase=PHASE_OUTPUT, fullPrimitiveEOSRelativeError=eos,
                  wallTemperatureMaximumErrorK=wall_temperature_error,
                  savedWallFluidTemperatureRangeK=[float(min(t[ids])), float(max(t[ids]))],
                  viscosityRelativeError=mu_error, heatGradientRelativeError=gradient_error,
                  heatFluxRangeWm2=[float(min(q)), float(max(q))],
                  temperatureRangeK=[float(min(t)), float(max(t))],
                  nativeHistoryHF=dict(value=native_heat, unit='W/m',
                                       meaning='planar dual-boundary sum in fixed SU2 source'),
                  nativePlanarIntegralRelativeError=float(planar_error),
                  revolvedPowerW=full_power, gaussPowerRelativeError=power_error,
                  revolvedPowerMethod='exact segment integral of linear q times linear r, 2*pi revolution',
                  yPlusPercentiles=dict(zip(['p50', 'p95', 'p99', 'max'],
                                           np.percentile(a[ids, columns.index('Y_Plus')], [50, 95, 99, 100]).tolist())),
                  finalResiduals={name: float(history[-1, labels.index(name)]) for name in
                                  ('rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]')},
                  conditions=c, elapsedSeconds=receipt['elapsedSeconds'],
                  sourceSHA256={str(p.relative_to(ROOT)): sha(p) for p in (
                      Path(__file__).resolve(), ROOT/'scripts/su2SSTInterface.py',
                      ROOT/'scripts/prepareSU2AxisymmetricBenchmark.py', folder/'engine_process.json')})
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    assert not output.exists()
    try:
        result = audit(args.folder)
    except Exception as error:
        save_json(output, dict(interfacePassed=False, errorType=type(error).__name__, error=str(error),
                               sourceSHA256=sha(Path(__file__))))
        raise
    save_json(output, result)
    print(json.dumps(result, ensure_ascii=False))
    raise SystemExit(0 if result['interfacePassed'] else 1)
