"""核对曲面候选的真实三维通量、固定壁面、全局运动投影和五量账；不授予稳态资格。"""

from pathlib import Path
import argparse
import json
import math

import numpy as np
from scipy.io import loadmat

from auditShockFittedALEKernel import exact_at, geometry
from auditShockFittedDynamics import primitive, shock_speed
from auditThermalPressureSurfaceFlow import sha


def euler_ale(Q, normal, speed, gamma):
    rho = Q[0]
    v = Q[1:4] / rho
    pressure = (gamma - 1) * (Q[4] - rho * float(v @ v) / 2)
    un = float(v @ normal)
    return np.r_[rho * un, rho * v * un + pressure * normal, (Q[4] + pressure) * un] - speed * Q


def last_flux(study):
    """模块一：激波直接使用上游F-DQ；墙用独立镜像根；出口检查全部特征出流。"""
    gas, up = study['gas'], study['upstream']
    mesh = study['cap']['shell']['mesh']
    Q = study['lastInputQ']
    dt = float(study['history'][-1, 1])
    geo = geometry(mesh, study['lastNodes0'], study['nodes'], dt)
    area = np.linalg.norm(geo['Abar'], axis=1)
    normals = geo['Abar'] / area[:, None]
    D = geo['swept'] / (dt * area)
    owner, neighbour = geo['owner'], geo['neighbour']
    bound = mesh['boundaryFaces'].astype(int) - 1
    row_by_face = {int(face): i for i, face in enumerate(bound)}
    kind = study['boundary']['kind']
    PL = primitive(Q, gas)
    U = np.asarray(up['velocityVector'])
    Qup = np.r_[up['density'], up['density'] * U,
                up['pressure'] / (gas['gamma'] - 1) + up['density'] * float(U @ U) / 2]
    flux = np.zeros((len(owner), 5))
    wall_pressures = np.full(len(bound), np.nan)
    shock_pressures = np.full(len(bound), np.nan)
    minimum_outlet_mach = math.inf
    for face, cell in enumerate(owner):
        left, n, w = PL[cell], normals[face], D[face]
        if neighbour[face] >= 0:
            flux[face] = exact_at(left, PL[neighbour[face]], n, w, gas)[0]
            continue
        row = row_by_face[face]
        if kind[row] == 2:
            right = left.copy()
            right[1:4] -= 2 * (float(left[1:4] @ n) - w) * n
            _, p = exact_at(left, right, n, w, gas)
            flux[face] = np.r_[0., p * n, p * w]
            wall_pressures[row] = p
        elif kind[row] == 3:
            flux[face] = euler_ale(Qup, n, w, gas['gamma'])
            relative = w - float(U @ n)
            mach2 = relative ** 2 / (gas['gamma'] * gas['R'] * up['temperature'])
            assert mach2 > 1
            ratio = (gas['gamma'] + 1) * mach2 / ((gas['gamma'] - 1) * mach2 + 2)
            shock_pressures[row] = up['pressure'] + up['density'] * relative ** 2 * (1 - 1 / ratio)
        elif kind[row] == 4:
            normal_mach = (float(left[1:4] @ n) - w) / math.sqrt(gas['gamma'] * gas['R'] * left[4])
            assert normal_mach > 1
            minimum_outlet_mach = min(minimum_outlet_mach, normal_mach)
            flux[face] = euler_ale(Q[cell], n, w, gas['gamma'])
        else:
            raise AssertionError('未列明的球鼻边界类型。')
    saved = study['lastStep']
    scale = np.maximum(1., np.max(abs(flux), axis=0))
    flux_error = float(np.max(abs(saved['faceALEFlux'] - flux) / scale))
    assert flux_error < 2e-11, flux_error
    np.testing.assert_allclose(saved['wallPressurePa'], wall_pressures, rtol=2e-11, atol=0., equal_nan=True)
    np.testing.assert_allclose(saved['shockPressurePa'], shock_pressures, rtol=2e-11, atol=0., equal_nan=True)
    balance = np.zeros_like(Q)
    integrated = area[:, None] * flux
    np.add.at(balance, owner, integrated)
    inside = neighbour >= 0
    np.add.at(balance, neighbour[inside], -integrated[inside])
    expected = (geo['V0'][:, None] * Q - dt * balance) / geo['V1'][:, None]
    U0 = np.linalg.norm(U)
    qscale = np.r_[up['density'], np.full(3, up['density'] * U0), up['density'] * U0 ** 2]
    state_error = float(np.max(abs(expected - study['Q']) / qscale))
    assert state_error < 2e-11, state_error
    history = study['history']
    integrals, boundary = study['conservedIntegrals'], study['boundaryFluxRates']
    assert len(integrals) == len(history) + 1 and len(boundary) == len(history)
    global_error = float(np.max(abs(np.diff(integrals, axis=0) + history[:, 1, None] * boundary) /
                                (geo['V0'].sum() * qscale)))
    assert global_error < 2e-11
    return dict(fluxRelativeError=flux_error, stateRelativeError=state_error,
                independentGCLRelative=geo['gcl'], allStepsBalanceByFreestreamScale=global_error,
                minimumOutletNormalMach=minimum_outlet_mach)


def motion_audit(study):
    """模块二：独立SVD求真实三节点面速度投影，并用Brent相容根核对请求速度。"""
    shell = study['cap']['shell']
    mesh = shell['mesh']
    X = study['lastNodes0']
    N = len(shell['originalVertexIDs'])
    layer = int(shell['layers'])
    outer = (mesh['boundaryFaces'].astype(int) - 1)[np.asarray(shell['outerBoundaryMask'], bool)]
    faces = mesh['faces'][outer].astype(int) - 1
    xyz = X[faces]
    vectors = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]) / 2
    area = np.linalg.norm(vectors, axis=1)
    normals = vectors / area[:, None]
    left = primitive(study['lastInputQ'][mesh['owner'][outer].astype(int) - 1], study['gas'])
    requested = np.array([shock_speed(p, n, study['upstream'], study['gas'])[0] for p, n in zip(left, normals)])
    actual = study['lastMotion']
    np.testing.assert_allclose(actual['requestedFaceNormalSpeedMS'], requested, rtol=2e-11, atol=2e-8)
    rays = X[layer * N:layer * N + N] - X[:N]
    direction = rays / np.linalg.norm(rays, axis=1)[:, None]
    C = np.zeros((len(faces), N))
    ids = faces - layer * N
    for i in range(3):
        C[np.arange(len(faces)), ids[:, i]] = np.sum(direction[ids[:, i]] * normals, axis=1) / 3
    A, b = np.sqrt(area)[:, None] * C, np.sqrt(area) * requested
    solution, _, rank, singular = np.linalg.lstsq(A, b, rcond=1e-12)
    node_error = float(np.max(abs(solution - actual['radialNodeSpeedMS'])))
    assert node_error < 2e-7
    normal_error = float(np.max(abs(C @ solution - actual['faceNormalSpeedMS'])))
    assert normal_error < 2e-7
    np.testing.assert_array_equal(actual['nodeVelocityMS'][:N], 0.)
    np.testing.assert_array_equal(study['nodes'][:N], study['initialNodes'][:N])
    np.testing.assert_allclose(X + study['history'][-1, 1] * actual['nodeVelocityMS'], study['nodes'], rtol=0, atol=1e-17)
    residual = A @ solution - b
    stationarity = float(np.linalg.norm(A.T @ residual) / (np.linalg.norm(A) * np.linalg.norm(b)))
    assert stationarity < 2e-10
    return dict(maximumNodeSpeedDifferenceMS=node_error, maximumFaceSpeedDifferenceMS=normal_error,
                matrixRank=int(rank), radialUnknowns=N, faces=len(faces),
                singularValueRatio=float(singular[-1] / singular[0]),
                normalEquationStationarityRelative=stationarity,
                weightedResidualRmsMS=float(np.linalg.norm(residual) / np.sqrt(area.sum())),
                maximumUnrepresentedFaceSpeedMS=float(np.max(abs(C @ solution - requested))))


def audit(evidence, process_file, destination):
    assert not destination.exists()
    process = json.loads(process_file.read_text(encoding='utf-8-sig'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
    saved = loadmat(evidence, simplify_cells=True)
    s = saved['study']
    assert not s['productionQualified'] and not s['shockShapeSolved'] and not s['heatFluxSolved']
    assert not s['wallTemperatureSolved'] and not s['processPeakMemoryMeasured'] and s['bodyVerticesFixed']
    history = s['history']
    assert history.shape[1] == 16 and np.all(history[:, 1] > 0)
    assert history[:, 6].min() >= s['gas']['temperatureRangeK'][0]
    assert history[:, 7].max() <= s['gas']['temperatureRangeK'][1]
    assert np.all(history[:, 8] > 0) and np.all(history[:, 12] > 1)
    # 候选稳态与生产资格分开；不能仅因一次守恒通过就授予稳态。
    all_steady = bool(history[-1, 2] < 5e-4 and history[-1, 3] < 1e-4 and history[-1, 4] < 1e-4)
    assert bool(s['candidateDiscreteSteady']) == all_steady
    up, gas = s['upstream'], s['gas']
    M2 = np.dot(up['velocityVector'], up['velocityVector']) / (gas['gamma'] * gas['R'] * up['temperature'])
    g = gas['gamma']
    p2 = up['pressure'] * (2 * g * M2 - (g - 1)) / (g + 1)
    downstream_mach2 = (1 + (g - 1) * M2 / 2) / (g * M2 - (g - 1) / 2)
    pitot = p2 * (1 + (g - 1) * downstream_mach2 / 2) ** (g / (g - 1))
    wall = s['lastStep']['wallPressurePa']
    pressure = wall[np.isfinite(wall)]
    report = dict(status='passed-conservative-curved-transient-audit-not-steady-qualification',
                  evidenceSHA256=sha(evidence), processSHA256=sha(process_file),
                  acceptedSteps=len(history), timeS=float(s['timeS']),
                  bodyFacets=len(s['cap']['F']), tetrahedra=len(s['Q']),
                  lastStep=last_flux(s), sharedNodeMotion=motion_audit(s),
                  finalResidual=float(history[-1, 2]), finalIdealShockSpeedByFreestream=float(history[-1, 3]),
                  finalProjectionResidualByFreestream=float(history[-1, 4]),
                  finalShockPressureMismatch=float(history[-1, 15]),
                  minimumTemperatureK=float(history[:, 6].min()), maximumTemperatureK=float(history[:, 7].max()),
                  wallPressureRangePa=[float(pressure.min()), float(pressure.max())],
                  normalShockIsentropicPitotPressurePa=float(pitot),
                  maximumTransientWallPressureDifferenceFromPitot=float(pressure.max() / pitot - 1),
                  pitotComparisonRole='解析正常激波后等熵驻止极限参考；有限面元/未稳态最大壁压不是驻点验证通过',
                  maximumFacetCentroidSagM=float(s['cap']['maximumFacetCentroidSagM']),
                  declaredRegionArrayUpperBytes=float(s['maximumDeclaredRegionArrayBytes']),
                  processPeakMemoryMeasured=False, candidateDiscreteSteady=all_steady,
                  productionQualified=False, newProductionSolvedAreaM2=0, experimentalValidation=False)
    destination.write_bytes((json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('evidence', type=Path)
    parser.add_argument('process', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    audit(args.evidence, args.process, args.destination)
