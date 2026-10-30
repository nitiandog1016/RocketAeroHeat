"""独立复算相容速度、解析Sod体平均、最后一步通量/几何/CFL和逐步五量账。"""

from pathlib import Path
import argparse
import json
import math

import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from scipy.special import beta

from auditShockFittedALEKernel import exact_at, geometry, pressure_curve, uniform_step
from auditThermalPressureSurfaceFlow import sha


def primitive(Q, gas):
    q = np.asarray(Q)
    rho = q[..., 0]
    v = q[..., 1:4] / rho[..., None]
    p = (gas['gamma'] - 1) * (q[..., 4] - rho * np.sum(v * v, axis=-1) / 2)
    return np.concatenate((rho[..., None], v, (p / (rho * gas['R']))[..., None]), axis=-1)


def shock_speed(left, normal, up, gas):
    right = np.r_[up['density'], up['velocityVector'], up['temperature']]
    _, pressure = exact_at(left, right, normal, 0., gas)
    sound = math.sqrt(gas['gamma'] * gas['R'] * up['temperature'])
    speed = float(np.asarray(up['velocityVector']) @ normal) + sound * math.sqrt(
        (gas['gamma'] + 1) / (2 * gas['gamma']) * pressure / up['pressure'] +
        (gas['gamma'] - 1) / (2 * gas['gamma']))
    return speed, pressure


def sod_reference(gamma):
    """模块一：Brent求星状态，直接在实验室相似坐标采样解析波系。"""
    rhoL, rhoR, pL, pR = 1., .125, 1e5, 1e4
    aL, aR = math.sqrt(gamma * pL / rhoL), math.sqrt(gamma * pR / rhoR)
    root = lambda p: pressure_curve(p, rhoL, pL, aL, gamma) + pressure_curve(p, rhoR, pR, aR, gamma)
    pressure = brentq(root, pR * 1e-8, pL, xtol=1e-10, rtol=1e-14)
    velocity = (pressure_curve(pressure, rhoR, pR, aR, gamma) -
                pressure_curve(pressure, rhoL, pL, aL, gamma)) / 2
    rhoSL = rhoL * (pressure / pL) ** (1 / gamma)
    z = (gamma - 1) / (gamma + 1)
    rhoSR = rhoR * (pressure / pR + z) / (z * pressure / pR + 1)
    head = -aL
    tail = velocity - aL * (pressure / pL) ** ((gamma - 1) / (2 * gamma))
    shock = aR * math.sqrt((gamma + 1) / (2 * gamma) * pressure / pR + (gamma - 1) / (2 * gamma))

    def point(xi):
        if xi < head:
            rho, u, p = rhoL, 0., pL
        elif xi < tail:
            a = 2 / (gamma + 1) * (aL - (gamma - 1) * xi / 2)
            u = 2 / (gamma + 1) * (aL + xi)
            rho, p = rhoL * (a / aL) ** (2 / (gamma - 1)), pL * (a / aL) ** (2 * gamma / (gamma - 1))
        elif xi < velocity:
            rho, u, p = rhoSL, velocity, pressure
        elif xi < shock:
            rho, u, p = rhoSR, velocity, pressure
        else:
            rho, u, p = rhoR, 0., pR
        return np.array([rho, rho * u, 0., 0., p / (gamma - 1) + rho * u * u / 2])

    return point, np.array([head, tail, velocity, shock]), dict(
        pressureStarPa=pressure, velocityStarMS=velocity, leftHeadMS=head, leftTailMS=tail,
        rightShockMS=shock, leftStarDensity=rhoSL, rightStarDensity=rhoSR)


def exact_cell_averages(study, gas):
    """模块二：按四面体X投影的Beta权重积分；波前处分段，不把重心值冒充体平均。"""
    point, waves, star = sod_reference(gas['gamma'])
    nodes = np.asarray(study['box']['nodes'])
    cells = np.asarray(study['box']['mesh']['tetrahedra'], int) - 1
    time = float(study['timeS'])
    breakpoints = .02 + time * waves
    gauss, weights = np.polynomial.legendre.leggauss(12)
    result = np.zeros((len(cells), 5))
    for i, ids in enumerate(cells):
        x = nodes[ids, 0]
        lo, hi = x.min(), x.max()
        right = int(np.sum(x == hi))
        assert hi > lo and right in (1, 2, 3) and np.all((x == lo) | (x == hi))
        edges = np.r_[lo, breakpoints[(breakpoints > lo) & (breakpoints < hi)], hi]
        for a, b in zip(edges[:-1], edges[1:]):
            sample = (a + b) / 2 + (b - a) / 2 * gauss
            z = (sample - lo) / (hi - lo)
            density = z ** (right - 1) * (1 - z) ** (3 - right) / (beta(right, 4 - right) * (hi - lo))
            values = np.array([point((v - .02) / time) for v in sample])
            result[i] += np.sum(values * (weights * density * (b - a) / 2)[:, None], axis=0)
    return result, star


def last_step(study, gas):
    """模块三：独立Gauss几何、实验室xi=D Riemann和五量更新，另核当前几何波扇CFL。"""
    mesh = study['box']['mesh']
    X0, X1 = np.asarray(study['lastNodes0']), np.asarray(study['nodes1'])
    dt = float(study['history'][-1, 1])
    geom = geometry(mesh, X0, X1, dt)
    area = np.linalg.norm(geom['Abar'], axis=1)
    normal = geom['Abar'] / area[:, None]
    moving = geom['swept'] / (dt * area)
    owner, next_cell = geom['owner'], geom['neighbour']
    bound = np.asarray(mesh['boundaryFaces'], int) - 1
    rows = {int(face): row for row, face in enumerate(bound)}
    initial = np.asarray(study['lastInputQ'])
    left_all = primitive(initial, gas)
    kind = np.atleast_1d(study['boundary']['kind'])
    exterior = np.asarray(study['boundary']['exteriorQ'])
    flux = np.zeros((len(owner), 5))
    radii = np.zeros(len(owner))
    pstars = np.zeros(len(owner))
    g = gas['gamma']
    boost = np.asarray(study['boost'])
    for face, cell in enumerate(owner):
        left = left_all[cell]
        right = left_all[next_cell[face]] if next_cell[face] >= 0 else left.copy()
        wall = False
        if next_cell[face] < 0:
            row = rows[face]
            if kind[row] == 1:
                right = primitive(exterior[row], gas)
            else:
                assert kind[row] == 2
                wall = True
                right[1:4] -= 2 * (left[1:4] @ normal[face] - moving[face]) * normal[face]
        value, star = exact_at(left, right, normal[face], moving[face], gas)
        flux[face] = np.r_[0., star * normal[face], star * moving[face]] if wall else value
        # CFL在步初几何/明确给定节点速度处定义，不能用差分节点的舍入速度替代输入。
        D = float(boost @ normal[face])
        pL, pR = left[0] * gas['R'] * left[4], right[0] * gas['R'] * right[4]
        aL, aR = math.sqrt(g * pL / left[0]), math.sqrt(g * pR / right[0])
        uL, uR = float(left[1:4] @ normal[face] - D), float(right[1:4] @ normal[face] - D)
        speedL = uL - aL * (math.sqrt((g + 1) / (2 * g) * star / pL + (g - 1) / (2 * g)) if star > pL else 1.)
        speedR = uR + aR * (math.sqrt((g + 1) / (2 * g) * star / pR + (g - 1) / (2 * g)) if star > pR else 1.)
        radii[face] = max(abs(uL) + aL, abs(uR) + aR, abs(speedL), abs(speedR))
        pstars[face] = star
    saved = study['lastStep']
    component_scale = np.maximum(1., np.max(abs(flux), axis=0))
    flux_error = float(np.max(abs(saved['faceALEFlux'] - flux) / component_scale))
    assert flux_error < 2e-11, flux_error
    divergence = np.zeros_like(initial)
    integrated = area[:, None] * flux
    np.add.at(divergence, owner, integrated)
    inside = next_cell >= 0
    np.add.at(divergence, next_cell[inside], -integrated[inside])
    expected = (geom['V0'][:, None] * initial - dt * divergence) / geom['V1'][:, None]
    state_scale = np.array([1., math.sqrt(g * 1e5), math.sqrt(g * 1e5), math.sqrt(g * 1e5), 1e5 / (g - 1)])
    state_error = float(np.max(abs(saved['Q'] - expected) / state_scale))
    assert state_error < 2e-11, state_error
    rates = np.zeros(len(initial))
    np.add.at(rates, owner, area * radii)
    np.add.at(rates, next_cell[inside], area[inside] * radii[inside])
    estimate = study['lastTimeStepEstimate']
    cfl_dt = float(estimate['CFL']) * np.min(geom['V0'] / rates)
    cfl_error = abs(cfl_dt / estimate['dt'] - 1)
    pressure_error = float(np.max(abs(pstars / estimate['riemannPressureStarPa'] - 1)))
    assert cfl_error < 2e-11 and pressure_error < 2e-11, (cfl_error, pressure_error)
    assert dt <= estimate['dt'] * (1 + 2e-14)
    history = np.asarray(study['history'])
    integrals = np.asarray(study['conservedIntegrals'])
    boundary = np.asarray(study['boundaryFluxRates'])
    physical_scale = geom['V0'].sum() * state_scale
    errors = abs(np.diff(integrals, axis=0) + history[:, 1, None] * boundary) / physical_scale
    balance = float(errors.max())
    assert balance < 2e-11 and history[:, 5].max() < 2e-11
    return dict(lastStepFluxRelativeError=flux_error, lastStepStateRelativeError=state_error,
                independentGCLRelative=geom['gcl'], cflRelativeError=cfl_error,
                cflPressureRootRelativeError=pressure_error, allStepsBalanceByPhysicalScale=balance)


def audit(evidence, process_file, destination):
    assert not destination.exists()
    process = json.loads(process_file.read_text(encoding='utf-8-sig'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
    data = loadmat(evidence, simplify_cells=True)
    assert isinstance(data['report'], dict) and not data['report']['productionQualified']
    gas, up = data['gas'], data['up']
    normals = np.asarray(data['normals'])
    known = primitive(data['known']['state']['Q'], gas)
    inverse = [shock_speed(left, n, up, gas)[0] for left, n in zip(known, normals)]
    inverse_error = float(np.max(abs(np.asarray(inverse) - data['result']['normalSpeedMS'])))
    assert inverse_error < 2e-8
    np.testing.assert_allclose(data['galilean']['normalSpeedMS'], np.asarray(inverse) + normals @ data['boost'], rtol=0, atol=2e-8)
    np.testing.assert_allclose(data['rotated']['normalSpeedMS'], inverse, rtol=0, atol=2e-8)
    tubes = []
    scales = np.array([1., math.sqrt(gas['gamma'] * 1e5), math.sqrt(gas['gamma'] * 1e5),
                       math.sqrt(gas['gamma'] * 1e5), 1e5 / (gas['gamma'] - 1)])
    for study in np.atleast_1d(data['tubes']):
        exact, star = exact_cell_averages(study, gas)
        nodes = study['box']['nodes']
        cells = np.asarray(study['box']['mesh']['tetrahedra'], int) - 1
        volume = np.linalg.det(nodes[cells[:, 1:]] - nodes[cells[:, 0, None]]) / 6
        errors = np.sum(volume[:, None] * abs(study['Q'] - exact), axis=0) / (volume.sum() * scales)
        velocities = primitive(study['Q'], gas)[:, 1:4]
        item = dict(nx=int(study['box']['counts'][0]), actualSteps=len(study['history']),
                    conservedCellAverageL1ByPhysicalScale=errors.tolist(),
                    transverseRmsVelocityMS=float(np.sqrt(np.sum(volume * np.sum(velocities[:, 1:] ** 2, axis=1)) / volume.sum())),
                    exactSod=star, **last_step(study, gas))
        tubes.append(item)
    errors = np.array([v['conservedCellAverageL1ByPhysicalScale'] for v in tubes])
    assert np.all(np.diff(errors[:, [0, 1, 4]], axis=0) < 0), errors
    shifted = last_step(data['shiftedTube'], gas)
    original = data['tubes'][0]['Q']
    rotated_momentum = original[:, 1:4] @ data['rotation'].T
    boost = data['boost']
    expected = np.column_stack((original[:, 0], rotated_momentum + original[:, 0, None] * boost,
                                original[:, 4] + rotated_momentum @ boost + original[:, 0] * (boost @ boost) / 2))
    frame_error = float(np.max(abs(data['shiftedTube']['Q'] - expected) / scales))
    assert frame_error < 2e-10
    cap = data['cap']
    moving = data['moving']
    projection = geometry(cap['shell']['mesh'], cap['shell']['nodes'],
                          cap['shell']['nodes'] + 2e-7 * moving['nodeVelocityMS'], 2e-7)
    np.testing.assert_array_equal(moving['nodeVelocityMS'][:len(cap['shell']['originalVertexIDs'])], 0.)
    box = data['validOutletBox']
    P = np.array([1., 600., 0., 0., 300.])
    Q = np.r_[P[0], P[0] * P[1:4], P[0] * gas['R'] * P[4] / (gas['gamma'] - 1) + P[0] * (P[1:4] @ P[1:4]) / 2]
    outlet = uniform_step(data['validOutlet'], box['mesh'], box['nodes'], box['nodes'],
                          min(1e-8, float(data['estimate']['dt'])), Q, gas)
    result = dict(status='passed-independent-Riemann-Sod-volume-average-ALE-CFL-and-conservation-audit',
                  evidenceSHA256=sha(evidence), processSHA256=sha(process_file), shockSpeedErrorMS=inverse_error,
                  tubes=tubes, transformedLastStep=shifted, frameErrorByPhysicalScales=frame_error,
                  nonuniformSharedNodeGCLRelative=projection['gcl'], supersonicUniformOutlet=outlet,
                  exactSodComparison='tetrahedron volume averages integrated with Beta(k,4-k) X-density; wave-front split Gauss quadrature',
                  productionQualified=False, newProductionSolvedAreaM2=0)
    destination.write_bytes((json.dumps(result, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('evidence', type=Path)
    parser.add_argument('process', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    audit(args.evidence, args.process, args.destination)
