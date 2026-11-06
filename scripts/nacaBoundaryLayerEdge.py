"""模块用途：从真实法向列提取受限的局部平台外缘；不向无定义区域填入来流。

本模块只实现零攻角附着层中可辨认平台的传统 δ99 诊断，不实现一般
三维分离、驻点极限或有显著熵层/曲率变化时的无黏外缘重构。
"""
from dataclasses import asdict, dataclass

import numpy as np


METHOD = 'local-normal-column-plateau-delta99-v1'
SCOPE = 'attached-laminar-near-two-dimensional-profile-with-resolved-local-plateau'
FIELDS = ('pressurePa', 'temperatureK', 'densityKgM3', 'enthalpyJkg',
          'specificHeatJkgK', 'viscosityPas', 'conductivityWmK')


@dataclass(frozen=True)
class EdgePolicy:
    """模块一：公开数值识别门槛；这些是诊断规则，不是实验误差或物理定理。"""
    minimum_points: int = 18
    minimum_plateau_points: int = 6
    plateau_distance_ratio: float = 1.5
    plateau_relative_spread: float = .005
    shock_pressure_log_change: float = .10
    shock_density_log_change: float = .05
    shock_stencil_edges: int = 3
    shock_guard_points: int = 2
    minimum_edge_speed_ms: float = 1.
    maximum_crossflow_ratio: float = .02
    maximum_normal_speed_ratio: float = .10
    reverse_flow_tolerance_ratio: float = 1e-5
    maximum_overshoot_ratio: float = .02
    wall_velocity_tolerance_ratio: float = 1e-5
    minimum_points_below_delta99: int = 8
    maximum_delta99_bracket_ratio: float = .20
    maximum_decimation_delta99_change: float = .05
    maximum_decimation_edge_change: float = .01


def _base(policy, provenance):
    return dict(method=METHOD, scope=SCOPE, policy=asdict(policy), provenance=dict(provenance),
                accepted=False, numericalExtractionAccepted=False, edge=None,
                physicalQualification=False, experimentalQualification=False,
                freestreamFallbackUsed=False)


def _reject(result, reason, **details):
    result.update(status=reason, **details)
    return result


def _relative_spread(values, scale=None):
    a = np.asarray(values, float)
    divisor = max(abs(float(np.median(a))), 1e-30) if scale is None else max(float(scale), 1e-30)
    return float(np.ptp(a)/divisor)


def _crossing(y, u, target, end):
    # 首次跨越必须持续至平台，禁止越过一处回落后仍冒称单一 δ99。
    crossings = np.flatnonzero((u[:end] < target) & (u[1:end+1] >= target))
    if len(crossings) != 1:
        return None
    lo = int(crossings[0]); hi = lo+1
    if np.any(u[hi:end+1] < target):
        return None
    value = y[lo]+(target-u[lo])*(y[hi]-y[lo])/(u[hi]-u[lo])
    return dict(distanceM=float(value), lowerIndex=lo, upperIndex=hi,
                bracketWidthM=float(y[hi]-y[lo]))


def extract_edge(profile, normal, tangent, station_s_m, provenance, policy=None,
                 evaluate_sensitivity=True):
    """模块二：检查、选平台、计算 δ99 并报告隔点采样敏感性。

    profile 必含 distanceM[N]、coordinatesM[N,3]、velocityMS[N,3] 和 FIELDS。
    第一项必须是真实壁点；normal/tangent 为该壁点外法向/流向单位向量。
    provenance 必含 stateKind（computed-solution 或 synthetic-test）、
    flowStateSource、normalColumnSource、thermodynamicModel。返回状态而不补零。
    edge 中五项可送给 heat_parameter；edgeSample 保留同一真实节点全部状态。
    δ99 的插值位置不同于 edgeSample 所选平台节点，二者绝不混称。
    """
    policy = EdgePolicy() if policy is None else policy
    result = _base(policy, provenance)
    if (provenance.get('stateKind') not in ('computed-solution', 'synthetic-test') or
            not all(provenance.get(k) for k in
                    ('flowStateSource', 'normalColumnSource', 'thermodynamicModel'))):
        return _reject(result, 'missing-local-solution-provenance')
    if any('freestream-substitut' in str(v).lower() for v in provenance.values()):
        return _reject(result, 'freestream-substitution-forbidden')
    if not np.isfinite(station_s_m) or station_s_m <= 0:
        return _reject(result, 'stagnation-limit-not-defined')
    try:
        y = np.asarray(profile['distanceM'], float)
        points = np.asarray(profile['coordinatesM'], float)
        velocity = np.asarray(profile['velocityMS'], float)
        fields = {k: np.asarray(profile[k], float) for k in FIELDS}
        n, t = np.asarray(normal, float), np.asarray(tangent, float)
    except (KeyError, TypeError, ValueError):
        return _reject(result, 'missing-or-malformed-profile')
    count = len(y) if y.ndim == 1 else 0
    if (count < policy.minimum_points or points.shape != (count, 3) or
            velocity.shape != (count, 3) or any(a.shape != (count,) for a in fields.values())):
        return _reject(result, 'insufficient-or-inconsistent-samples')
    if (n.shape != (3,) or t.shape != (3,) or
            not all(np.isfinite(a).all() for a in [y, points, velocity, n, t]+list(fields.values()))):
        return _reject(result, 'nonfinite-or-malformed-profile')
    if any(np.any(fields[k] <= 0) for k in FIELDS if k != 'enthalpyJkg'):
        return _reject(result, 'nonpositive-thermodynamic-or-transport-state')
    if (abs(np.linalg.norm(n)-1) > 1e-10 or abs(np.linalg.norm(t)-1) > 1e-10 or
            abs(np.dot(n,t)) > 1e-10):
        return _reject(result, 'invalid-local-basis')
    if y[0] != 0 or np.any(np.diff(y) <= 0):
        return _reject(result, 'missing-wall-or-unordered-normal-distance')
    geometric_error = float(np.max(np.linalg.norm(points-points[0]-y[:,None]*n, axis=1)))
    if geometric_error > max(1e-10, 1e-8*y[-1]):
        return _reject(result, 'not-a-straight-normal-column', maximumCoordinateErrorM=geometric_error)
    u = velocity @ t
    vn = velocity @ n
    cross = velocity-u[:,None]*t-vn[:,None]*n
    h0 = fields['enthalpyJkg']+.5*np.sum(velocity*velocity, axis=1)
    # 总焓差的尺度取 cp*T+动能，避免任意焓零点改变平台判定。
    hscale = fields['specificHeatJkgK']*fields['temperatureK']+.5*np.sum(velocity*velocity, axis=1)

    # 模块三：压力与密度同向显著变化视为激波/强外场变化屏障。
    # 该保守规则不证明激波身份；屏障后数据均不得作为当地边界层外缘。
    lp, lr = np.log(fields['pressurePa']), np.log(fields['densityKgM3'])
    barriers = []
    for i in range(1, count-1):
        j = min(count-1, i+policy.shock_stencil_edges)
        dp, dr = lp[j]-lp[i], lr[j]-lr[i]
        if (dp*dr > 0 and abs(dp) >= policy.shock_pressure_log_change and
                abs(dr) >= policy.shock_density_log_change):
            barriers.append(i)
    stop = max(1, barriers[0]-policy.shock_guard_points) if barriers else count
    result.update(normalSamples=count, stationS_m=float(station_s_m),
                  strongOuterChangeFirstIndex=int(barriers[0]) if barriers else None,
                  usableLastIndex=stop-1, maximumCoordinateErrorM=geometric_error)
    if stop < policy.minimum_points:
        return _reject(result, 'outer-change-overlaps-unresolved-boundary-layer')
    speed_scale = float(np.max(u[:stop]))
    if speed_scale < policy.minimum_edge_speed_ms:
        return _reject(result, 'stagnant-or-near-stagnation-column')
    if np.min(u[:stop]) < -policy.reverse_flow_tolerance_ratio*speed_scale:
        return _reject(result, 'reverse-flow-outside-method-scope')
    if np.linalg.norm(velocity[0]) > policy.wall_velocity_tolerance_ratio*speed_scale:
        return _reject(result, 'wall-no-slip-not-satisfied')
    if np.max(np.linalg.norm(cross[:stop], axis=1))/speed_scale > policy.maximum_crossflow_ratio:
        return _reject(result, 'crossflow-outside-method-scope')

    # 模块四：仅把有足够点数、足够法向跨度、全部状态稳定的区间作为平台。
    # 不用相邻一个网格层的微小梯度，也不用最外侧/最大速度直接指定 Ue。
    windows = []
    spread_fields = dict(tangentialVelocity=u, pressure=fields['pressurePa'],
                         temperature=fields['temperatureK'], density=fields['densityKgM3'],
                         viscosity=fields['viscosityPas'], conductivity=fields['conductivityWmK'])
    for i in range(2, stop-policy.minimum_plateau_points+1):
        j = max(i+policy.minimum_plateau_points-1,
                int(np.searchsorted(y, y[i]*policy.plateau_distance_ratio)))
        if j >= stop:
            continue
        span = slice(i,j+1)
        spreads = {k:_relative_spread(a[span]) for k,a in spread_fields.items()}
        spreads['totalEnthalpyThermalScaled'] = _relative_spread(h0[span], np.median(hscale[span]))
        if max(spreads.values()) <= policy.plateau_relative_spread:
            windows.append([i,j])
    if not windows:
        return _reject(result, 'no-resolved-local-platform')
    groups = []
    for lo, hi in windows:
        if groups and lo <= groups[-1][1]:
            groups[-1][1] = max(groups[-1][1], hi)
        else:
            groups.append([lo,hi])
    result['candidatePlatforms'] = [dict(lowerIndex=a,upperIndex=b,
        lowerDistanceM=float(y[a]),upperDistanceM=float(y[b])) for a,b in groups]
    if len(groups) != 1:
        return _reject(result, 'multiple-platforms-ambiguous')
    lo, hi = groups[0]
    span = slice(lo,hi+1)
    spreads = {k:_relative_spread(a[span]) for k,a in spread_fields.items()}
    spreads['totalEnthalpyThermalScaled'] = _relative_spread(h0[span], np.median(hscale[span]))
    result['platformRelativeSpread'] = spreads
    # 一串局部缓变窗口可能组成明显漂移的大区间；不能用传递性伪造平台。
    if max(spreads.values()) > policy.plateau_relative_spread:
        return _reject(result, 'joined-platform-has-systematic-drift')
    index = int(lo+np.argmin(abs(y[lo:hi+1]-np.sqrt(y[lo]*y[hi]))))
    ue = float(u[index])
    if ue < policy.minimum_edge_speed_ms:
        return _reject(result, 'stagnant-or-near-stagnation-platform')
    if np.max(abs(vn[:hi+1]))/ue > policy.maximum_normal_speed_ratio:
        return _reject(result, 'normal-velocity-outside-method-scope')
    if np.max(u[:hi+1])/ue > 1+policy.maximum_overshoot_ratio:
        return _reject(result, 'velocity-overshoot-outside-method-scope')
    crossings = {str(f):_crossing(y,u,f*ue,hi) for f in (.985,.99,.995)}
    result['velocityFractionCrossings'] = crossings
    if any(v is None for v in crossings.values()):
        return _reject(result, 'missing-or-nonunique-sustained-velocity-crossing')
    delta = crossings['0.99']
    if (delta['upperIndex'] < policy.minimum_points_below_delta99 or
            delta['bracketWidthM']/delta['distanceM'] > policy.maximum_delta99_bracket_ratio):
        return _reject(result, 'delta99-under-resolved')
    if crossings['0.995']['distanceM'] > y[lo]:
        return _reject(result, 'platform-not-outside-velocity-layer')

    # 模块五：取真实平台中点的同一状态，绝不分别平均 rho/T/mu/k 而破坏相容性。
    sample = {k:float(a[index]) for k,a in fields.items()}
    sample.update(index=index,distanceM=float(y[index]),coordinatesM=points[index].tolist(),
                  velocityMS=velocity[index].tolist(),tangentialSpeedMS=ue,
                  totalSpeedMS=float(np.linalg.norm(velocity[index])),
                  totalEnthalpyJkg=float(h0[index]))
    result.update(edgeSample=sample, delta99M=delta['distanceM'],
                  definitionSensitivityM=[crossings['0.985']['distanceM'],crossings['0.995']['distanceM']])
    edge = dict(temperatureK=sample['temperatureK'], densityKgM3=sample['densityKgM3'],
                speedMS=sample['totalSpeedMS'], viscosityPas=sample['viscosityPas'],
                conductivityWmK=sample['conductivityWmK'])
    if evaluate_sensitivity:
        keep = np.unique(np.r_[np.arange(0,count,2),count-1]).astype(int)
        reduced = {k:np.asarray(profile[k])[keep] for k in
                   ('distanceM','coordinatesM','velocityMS')+FIELDS}
        other = extract_edge(reduced,n,t,station_s_m,provenance,policy,False)
        sensitivity = dict(kind='every-second-sample-with-wall-and-outer-end-preserved',
                           accepted=other['accepted'], status=other['status'])
        if other['accepted']:
            dd = abs(other['delta99M']/result['delta99M']-1)
            de = max(abs(other['edge'][k]/edge[k]-1) for k in edge)
            sensitivity.update(delta99RelativeChange=float(dd),maximumEdgeRelativeChange=float(de),
                               edgeSampleDistanceM=other['edgeSample']['distanceM'])
            result['samplingSensitivity'] = sensitivity
            if (dd > policy.maximum_decimation_delta99_change or
                    de > policy.maximum_decimation_edge_change):
                return _reject(result,'normal-sampling-sensitivity-failed')
        else:
            result['samplingSensitivity'] = sensitivity
            return _reject(result,'normal-sampling-sensitivity-unresolved')
    else:
        result['samplingSensitivity'] = dict(kind='not-requested',accepted=None)
    result.update(status='accepted-local-platform-diagnostic', accepted=True,
                  numericalExtractionAccepted=True, edge=edge)
    return result
