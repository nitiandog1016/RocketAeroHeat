"""公开压缩拐角的逐项预设验收与跨检查点稳定性，缺失量保留原因。"""
import numpy as np


def stability(frames, final_iteration, wall_area):
    """最后至少1000更新的所有保存帧与末帧逐站比较，跨度不足保持待验。"""
    valid = sorted(frames, key=lambda item: item['iteration'])
    if not valid or valid[-1]['iteration'] != final_iteration:
        return dict(available=False, reason='last recorded state has no matching saved frame')
    before = [frame for frame in valid if frame['iteration'] <= final_iteration-1000]
    if not before:
        return dict(available=False, reason='saved frames span less than 1000 updates')
    first = before[-1]['iteration']
    selected = [frame for frame in valid if frame['iteration'] >= first]
    if any(b['iteration']-a['iteration'] > 100 for a, b in zip(selected, selected[1:])):
        return dict(available=False, reason='checkpoint gap exceeds the planned 100 updates')
    final = selected[-1]
    pscale = np.maximum(np.abs(final['p']), 576.)
    qscale = np.maximum(np.abs(final['q']), 9300.)
    power_scale = max(abs(final['power']), 9300.*wall_area)
    pressure = max(float(np.max(abs(np.asarray(f['p'])-final['p'])/pscale)) for f in selected)
    heat = max(float(np.max(abs(np.asarray(f['q'])-final['q'])/qscale)) for f in selected)
    power = max(abs(f['power']-final['power'])/power_scale for f in selected)
    return dict(available=True, firstIteration=first, lastIteration=final_iteration,
                spanUpdates=final_iteration-first, frames=len(selected),
                maximumPressureChange=pressure, maximumHeatChange=heat,
                maximumPowerChange=power)


def acceptance(current, residuals, steady, run_succeeded, field_is_final):
    """把事前计划的数值与实验门槛逐项展开，网格/外域资格另行保留。"""
    gates = []

    def gate(name, value, limit, strict=False, missing='current saved field unavailable'):
        available = value is not None and np.isfinite(value)
        passed = bool(value < limit if strict else value <= limit) if available else None
        gates.append(dict(name=name, value=value, limit=limit,
                          comparison='<' if strict else '<=', passed=passed,
                          reason=None if available else missing))

    gate('native successful terminal', 0. if run_succeeded else 1., 0.)
    gate('saved field matches final history', 0. if field_is_final else None, 0.)
    for name, value in residuals.items():
        gate(name, value, -8., True)
    values = [
        ('same-stage EOS', ('EOSRelativeError',), 2e-12, True),
        ('wall heat gradient', ('heatGradientRelativeError',), 2e-6, True),
        ('normalized density residual', ('residuals', 'normalizedRMS', 'Density'), 1e-5, True),
        ('normalized energy residual', ('residuals', 'normalizedRMS', 'Energy'), 1e-5, True),
        ('mass balance', ('boundaryBalance', 'massRelativeImbalance'), .01, True),
        ('energy balance', ('boundaryBalance', 'energyRelativeImbalance'), .01, True),
        ('y+ 99th percentile', ('yPlus99',), 1., False),
        ('pressure station RMS', ('stations', 'pressure', 'rmsRelativeError'), .1, False),
        ('pressure station maximum', ('stations', 'pressure', 'maximumRelativeError'), .2, False),
        ('heat station RMS', ('stations', 'heat', 'rmsRelativeError'), .2, False),
        ('heat station maximum', ('stations', 'heat', 'maximumRelativeError'), .3, False),
        ('heat peak magnitude', ('measuredStationHeatPeak', 'relativeMagnitudeError'), .2, False),
        ('heat peak position m', ('measuredStationHeatPeak', 'positionErrorM'), .01, False),
        ('upstream velocity profile', ('profile', 'velocityRatioRMSError'), .05, False),
        ('upstream total temperature profile', ('profile', 'totalTemperatureRatioRMSError'), .03, False)]
    for name, path, limit, strict in values:
        value = current
        for key in path:
            value = value.get(key) if isinstance(value, dict) else None
        gate(name, value, limit, strict)
    reference = current['localReference'] if current else {}
    for name, key, reference_value, limit in (
            ('upstream pressure reference', 'pressurePa', 576., .05),
            ('upstream heat reference', 'heatFluxWm2', 9300., .2)):
        value = abs(reference[key]/reference_value-1) if key in reference else None
        gate(name, value, limit)
    for name, key, limit in (
            ('pressure stability', 'maximumPressureChange', .01),
            ('heat stability', 'maximumHeatChange', .02),
            ('revolved power stability', 'maximumPowerChange', .02)):
        gate(name, steady.get(key), limit, True, steady.get('reason'))
    return dict(gates=gates, passed=all(item['passed'] is True for item in gates),
                pending=[item['name'] for item in gates if item['passed'] is None],
                failed=[item['name'] for item in gates if item['passed'] is False],
                furtherQualificationRequired=['grid sensitivity', 'outer-domain sensitivity',
                                              'model and material-property applicability'])
