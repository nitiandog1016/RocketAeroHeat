"""轴对称公开算例的控制体、物理通量、壁面测站及剖面后处理。"""
import csv
from pathlib import Path

import numpy as np

from su2SSTInterface import state, wall_geometry


def mesh_data(path):
    """读取原四边形和全部有向边界，边界右法向指向流体域外。"""
    with Path(path).open() as f:
        assert f.readline().strip() == 'NDIME= 2'
        count = int(f.readline().split('=')[1])
        rows = np.array([list(map(int, f.readline().split())) for _ in range(count)])
        assert np.all(rows[:, 0] == 9)
        n = int(f.readline().split('=')[1])
        p = np.array([list(map(float, f.readline().split())) for _ in range(n)])
        assert np.array_equal(p[:, 2], np.arange(n))
        markers = {}
        for _ in range(int(f.readline().split('=')[1])):
            tag = f.readline().split('=')[1].strip()
            size = int(f.readline().split('=')[1])
            e = np.array([list(map(int, f.readline().split())) for _ in range(size)])
            assert np.all(e[:, 0] == 3)
            markers[tag] = e[:, 1:]
    return p[:, :2], rows[:, 1:], markers


def median_volumes(points, cells):
    """原四边形均值中心/边中点构成中位子域，分别累计平面面积和旋转体积。"""
    p = points[cells]
    center = np.mean(p, axis=1)
    area = np.zeros(len(points))
    volume = np.zeros(len(points))
    for j in range(4):
        polygon = np.stack((p[:, j], (p[:, j]+p[:, (j+1)%4])/2,
                            center, (p[:, j]+p[:, (j-1)%4])/2), axis=1)
        nxt = np.roll(polygon, -1, axis=1)
        cross = polygon[:, :, 0]*nxt[:, :, 1]-nxt[:, :, 0]*polygon[:, :, 1]
        local_area = np.sum(cross, axis=1)/2
        local_volume = np.pi*np.sum((polygon[:, :, 1]+nxt[:, :, 1])*cross, axis=1)/3
        assert min(local_area) > 0 and min(local_volume) > 0
        np.add.at(area, cells[:, j], local_area)
        np.add.at(volume, cells[:, j], local_volume)
    return area, volume


def wls_gradient(points, cells, values):
    """按原单元棱边邻接做距离加权最小二乘，批量求多个场量的梯度。"""
    values = np.asarray(values)
    if values.ndim == 1:
        values = values[:, None]
    pairs = np.concatenate([cells[:, [j, (j+1)%4]] for j in range(4)])
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    pairs = np.concatenate((pairs, pairs[:, ::-1]))
    i, j = pairs[:, 0], pairs[:, 1]
    delta = points[j]-points[i]
    length = np.linalg.norm(delta, axis=1)
    unit = delta/length[:, None]
    matrix = np.zeros((len(points), 2, 2))
    rhs = np.zeros((len(points), 2, values.shape[1]))
    np.add.at(matrix, i, unit[:, :, None]*unit[:, None, :])
    np.add.at(rhs, i, unit[:, :, None]*(values[j]-values[i])[:, None, :]/length[:, None, None])
    assert np.min(np.linalg.eigvalsh(matrix)) > 1e-12
    return np.linalg.solve(matrix, rhs).transpose(0, 2, 1)


def integrate_vector(points, segments, nodal_flux):
    """对分段线性的物理通量做精确整周积分，同时保存流入/流出绝对量。"""
    p = points[segments]
    d = p[:, 1]-p[:, 0]
    outward = np.column_stack((d[:, 1], -d[:, 0]))
    f = nodal_flux[segments]
    total = incoming = outgoing = 0.
    for s in ((1-1/np.sqrt(3))/2, (1+1/np.sqrt(3))/2):
        radius = (1-s)*p[:, 0, 1]+s*p[:, 1, 1]
        normal_flux = np.sum(((1-s)*f[:, 0]+s*f[:, 1])*outward, axis=1)
        contribution = np.pi*radius*normal_flux
        total += float(np.sum(contribution))
        incoming += float(np.sum(np.maximum(-contribution, 0)))
        outgoing += float(np.sum(np.maximum(contribution, 0)))
    return dict(net=total, incoming=incoming, outgoing=outgoing)


def integrate_scalar(points, segments, values):
    """标量按线段弧长和半径积分，供壁面热功率与完整壁面积使用。"""
    p = points[segments]
    length = np.linalg.norm(p[:, 1]-p[:, 0], axis=1)
    q = values[segments]
    r = p[:, :, 1]
    return float(np.pi/3*np.dot(length, 2*q[:, 0]*r[:, 0]+q[:, 0]*r[:, 1]
                               +q[:, 1]*r[:, 0]+2*q[:, 1]*r[:, 1]))


def fields(path, spatial=True):
    """具名字段保持原全局编号，并核对流动阶段primitive与守恒速度。"""
    a, labels, nd, rho, u, thermal = state(path, True, True, spatial)
    assert nd == 2
    value = {name: a[:, i] for i, name in enumerate(labels)}
    value['velocity'] = u
    value['positiveThermalEnergy'] = thermal
    saved_u = np.column_stack((value['Velocity_x'], value['Velocity_y']))
    assert np.max(np.abs(saved_u-u)/np.maximum(np.abs(u), 1)) < 2e-12
    return value


def boundary_balance(points, cells, markers, f, conditions):
    """独立物理边界积分，包含轴对称应力、开放边界导热与静止壁面热流。"""
    rho = f['Density']
    u = f['velocity']
    temperature = f['Temperature']
    cp = conditions['gamma']*conditions['gasConstant']/(conditions['gamma']-1)
    grad = wls_gradient(points, cells, np.column_stack((u, temperature)))
    divergence = grad[:, 0, 0]+grad[:, 1, 1]+u[:, 1]/points[:, 1]
    effective_mu = f['Laminar_Viscosity']+f['Eddy_Viscosity']
    stress = effective_mu[:, None, None]*(grad[:, :2]+grad[:, :2].transpose(0, 2, 1))
    stress -= (2*effective_mu*divergence/3)[:, None, None]*np.eye(2)
    conductivity = cp*(f['Laminar_Viscosity']/.72+f['Eddy_Viscosity']/.9)
    mass = rho[:, None]*u
    energy = (f['Energy']+f['Pressure'])[:, None]*u
    energy -= np.einsum('nij,nj->ni', stress, u)+conductivity[:, None]*grad[:, 2]
    by_marker = {}
    for name, segments in markers.items():
        m = integrate_vector(points, segments, mass)
        e = integrate_vector(points, segments, energy)
        if name == 'wall':
            power = integrate_scalar(points, segments, f['Heat_Flux'])
            e = dict(net=power, incoming=max(-power, 0), outgoing=max(power, 0))
        by_marker[name] = dict(massKgS=m, energyW=e)
    mass_in = sum(row['massKgS']['incoming'] for row in by_marker.values())
    energy_in = sum(row['energyW']['incoming'] for row in by_marker.values())
    mass_net = sum(row['massKgS']['net'] for row in by_marker.values())
    energy_net = sum(row['energyW']['net'] for row in by_marker.values())
    return dict(markers=by_marker, massNetKgS=mass_net, energyNetW=energy_net,
                massRelativeImbalance=abs(mass_net)/max(mass_in, 1e-30),
                energyRelativeImbalance=abs(energy_net)/max(energy_in, 1e-30),
                definition='independent nodal physical flux quadrature with WLS viscous fluxes')


def residual_scales(points, cells, f, conditions):
    """原生二维积分RHS除以控制面积与来流尺度，RMS按旋转体积加权。"""
    area, volume = median_volumes(points, cells)
    c = conditions
    speed = c['speedMS']
    rho = c['densityKgM3']
    cp = c['gamma']*c['gasConstant']/(c['gamma']-1)
    k = 1.5*(speed*c['turbulenceIntensity'])**2
    h0 = cp*c['temperatureK']+.5*speed**2+k
    normalized = {}
    raw = {}
    for name, factor in (('Density', rho*speed), ('Energy', rho*speed*h0)):
        residual = f['Residual_'+name]
        scaled = residual/(area*factor)
        normalized[name] = float(np.sqrt(np.dot(volume, scaled**2)/sum(volume)))
        raw[name] = float(np.log10(np.sqrt(np.mean(residual**2))))
    return dict(referenceLengthM=1., fullRevolvedVolumeM3=float(sum(volume)),
                normalizedRMS=normalized, nativeRawLogRMS=raw)


def wall_values(points, markers, f):
    """沿下游顺序保存表面S、压力、正向入壁热流、壁流体温度和有符号剪切。"""
    ids = np.unique(markers['wall'])
    ids = ids[np.argsort(-points[ids, 0])]
    x = -points[ids, 0]
    s = np.where(x < 0, x, x/np.cos(np.deg2rad(20)))
    tangent = np.gradient(points[ids], axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    shear = np.column_stack((f['Skin_Friction_Coefficient_x'][ids],
                             f['Skin_Friction_Coefficient_y'][ids]))
    return dict(ids=ids, S=s, x=x, radius=points[ids, 1], p=f['Pressure'][ids],
                q=f['Heat_Flux'][ids], fluidT=f['Temperature'][ids],
                yplus=f['Y_Plus'][ids], cf=np.sum(shear*tangent, axis=1))


def sample_stations(wall, contract):
    """对原实验弧长测站插值，原重复压力行按合同排除。"""
    pressure = [r for r in contract['measurements']['pressure_rows'] if r['independent_comparison_point']]
    heat = contract['measurements']['heating_rows']
    records = {}
    for name, rows, field, key in (
            ('pressure', pressure, 'p', 'pressure_Pa'),
            ('heat', heat, 'q', 'reported_heating_magnitude_W_m2')):
        positions = np.array([r['S_m'] for r in rows])
        assert min(positions) >= min(wall['S']) and max(positions) <= max(wall['S'])
        predicted = np.interp(positions, wall['S'], wall[field])
        measured = np.array([r[key] for r in rows])
        relative = (predicted-measured)/measured
        records[name] = dict(S=positions, predicted=predicted, measured=measured,
                             signedRelativeError=relative,
                             rmsRelativeError=float(np.sqrt(np.mean(relative**2))),
                             maximumRelativeError=float(np.max(np.abs(relative))),
                             meanRelativeBias=float(np.mean(relative)))
    return records


def upstream_profile(points, f, conditions, contract):
    """在source_x=-0.06m沿真实原网格逐径向层插值，比较总温比和轴向速度比。"""
    first = np.arange(81*101).reshape(101, 81)
    cp = conditions['gamma']*conditions['gasConstant']/(conditions['gamma']-1)
    total_t = f['Temperature']+np.sum(f['velocity']**2, axis=1)/(2*cp)
    radius, velocity, total = [], [], []
    for row in first:
        order = row[np.argsort(points[row, 0])]
        assert min(points[order, 0]) <= .06 <= max(points[order, 0])
        radius.append(np.interp(.06, points[order, 0], points[order, 1]))
        velocity.append(np.interp(.06, points[order, 0], -f['velocity'][order, 0]))
        total.append(np.interp(.06, points[order, 0], total_t[order]))
    radial_distance = np.array(radius)-radius[0]
    data = contract['measurements']['upstream_profile']
    y = np.array([r['y_from_wall_m'] for r in data])
    reference = contract['experimental_local_reference']
    tref = reference['temperature_K']+reference['velocity_m_s']**2/(2*cp)
    predicted_u = np.interp(y, radial_distance, velocity)/reference['velocity_m_s']
    predicted_t = np.interp(y, radial_distance, total)/tref
    measured_u = np.array([r['velocity_ratio'] for r in data])
    measured_t = np.array([r['total_temperature_ratio'] for r in data])
    return dict(y=y, predictedU=predicted_u, measuredU=measured_u,
                predictedT0=predicted_t, measuredT0=measured_t, referenceT0K=tref,
                velocityRatioRMSError=float(np.sqrt(np.mean((predicted_u-measured_u)**2))),
                totalTemperatureRatioRMSError=float(np.sqrt(np.mean((predicted_t-measured_t)**2))))
