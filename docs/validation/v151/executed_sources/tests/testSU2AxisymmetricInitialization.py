"""独立检查二维相容初场、过零积分及公开算例验收的拒绝边界。"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareSU2AxisymmetricInitialState import initial_state, segment_projection
from su2AxisymmetricAnalysis import signed_segment_integral
from su2BenchmarkAcceptance import acceptance, stability


def verify():
    # 模块一，闭线段投影的内部点和两个端点均有解析距离。
    endpoints = np.array([[0., .1], [1., .1]])
    probe = np.array([[.5, .2], [-.2, .1], [1.2, .3]])
    distance, index, weight, feet = segment_projection(probe, endpoints, np.array([[0, 1]]))
    assert np.max(abs(distance-np.array([.1, .2, np.sqrt(.08)]))) < 1e-14
    assert np.array_equal(index, np.zeros(3))
    assert np.max(abs(weight-np.array([.5, 0, 1]))) < 1e-14
    assert np.max(abs(feet-np.array([[.5, .1], [0, .1], [1, .1]]))) < 1e-14
    points = np.array([[0., .1], [1., .1], [0., .2], [1., .2], [.5, .105]])
    markers = dict(wall=np.array([[0, 1]]), inlet=np.array([[1, 3]]), top=np.array([[2, 3]]))
    c = dict(gamma=1.4, gasConstant=287.058, speedMS=1274.902422634943,
             pressurePa=550., temperatureK=80., wallTemperatureK=311.,
             turbulenceIntensity=.00021, eddyToLaminarRatio=.001, unitReynoldsPerM=5.574e6)
    c['densityKgM3'] = c['pressurePa']/(c['gasConstant']*c['temperatureK'])
    data, detail = initial_state(points, markers, np.array([0, 1]), np.array([1e-5, 1e-5]), c)
    rho = data[:, 0]
    velocity = data[:, 1:3]/rho[:, None]
    thermal = data[:, 3]/rho-.5*np.sum(velocity**2, axis=1)-data[:, 4]
    reconstructed_t = thermal*(c['gamma']-1)/c['gasConstant']
    energy_error = float(max(abs(reconstructed_t/detail['temperature']-1)))
    assert energy_error < 2e-14
    assert np.all(velocity[:2] == 0) and np.all(data[:2, 4] == 0)
    assert np.allclose(reconstructed_t[:2], 311., rtol=0, atol=1e-12)
    assert np.all(velocity[[2, 3], 0] == -c['speedMS'])
    assert min(data[:, 5]) > 0
    cp = c['gamma']*c['gasConstant']/(c['gamma']-1)
    h0 = cp*reconstructed_t+.5*np.sum(velocity**2, axis=1)+data[:, 4]
    blend = detail['blend']
    expected = (1-blend)*cp*311+blend*(cp*80+c['speedMS']**2/2+detail['kInfinity'])
    assert max(abs(h0/expected-1)) < 2e-14

    # 模块二，r=1+t，f=-1+2t在t=0.5变号，分别使用解析原函数。
    flux = signed_segment_integral(np.array([[1., 2.]]), np.array([[-1., 1.]]))
    assert abs(flux['net']-np.pi/3) < 1e-14
    assert abs(flux['incoming']-7*np.pi/12) < 1e-14
    assert abs(flux['outgoing']-11*np.pi/12) < 1e-14

    # 模块三，真实跨度/缺帧及严格门槛都有拒绝断言。
    frames = [dict(iteration=99+100*i, p=np.array([500.]), q=np.array([10000.]), power=100.)
              for i in range(12)]
    steady = stability(frames, 1199, 1.)
    assert steady['available'] and steady['spanUpdates'] == 1000
    assert not stability(frames[-10:], 1199, 1.)['available']
    assert not stability(frames, 1200, 1.)['available']
    assert not stability(frames[:5]+frames[6:], 1199, 1.)['available']
    current = dict(EOSRelativeError=0., heatGradientRelativeError=0., yPlus99=.5,
                   residuals=dict(normalizedRMS=dict(Density=0., Energy=0.)),
                   boundaryBalance=dict(massRelativeImbalance=0., energyRelativeImbalance=0.),
                   stations={key: dict(rmsRelativeError=0., maximumRelativeError=0.)
                             for key in ('pressure', 'heat')},
                   measuredStationHeatPeak=dict(relativeMagnitudeError=0., positionErrorM=0.),
                   profile=dict(velocityRatioRMSError=0., totalTemperatureRatioRMSError=0.),
                   localReference=dict(pressurePa=576., heatFluxWm2=9300.))
    residuals = {name: -9. for name in ('rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]')}
    good = acceptance(current, residuals, steady, True, True)
    assert good['passed'] and not good['pending'] and not good['failed']
    current['residuals']['normalizedRMS']['Density'] = 1e-5
    bad = acceptance(current, residuals, steady, True, True)
    assert bad['failed'] == ['normalized density residual']
    missing = acceptance(None, residuals, dict(available=False, reason='missing'), False, False)
    assert not missing['passed'] and len(missing['pending']) == 21
    return dict(passed=True, energyRelativeError=energy_error, signedIntegral=flux,
                gateCount=len(good['gates']), missingGateCount=len(missing['pending']),
                scope='algebra, geometry, exact signed integration and acceptance rules')


if __name__ == '__main__':
    print(json.dumps(verify()))
