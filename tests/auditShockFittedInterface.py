"""独立核对激波通量、节点重构和实际候选推进；不授予工程区域资格。"""

from pathlib import Path
import argparse
import json
import math

import numpy as np
from scipy.io import loadmat
from scipy.linalg import lstsq
from scipy.optimize import brentq

from auditShockFittedALEKernel import exact_at, geometry
from auditShockFittedCap import euler_ale, motion_audit
from auditShockFittedDynamics import primitive, shock_speed
from auditThermalPressureSurfaceFlow import sha


def nodal_motion(study):
    """模块一：用带列主元QR和Brent根，复核共享节点仿射状态及局部J+相容。"""
    shell, gas = study['cap']['shell'], study['gas']
    mesh, X, Q = shell['mesh'], study['lastNodes0'], study['lastInputQ']
    actual = study['lastMotion']
    count, layer = len(shell['originalVertexIDs']), int(shell['layers'])
    node_ids = layer * count + np.arange(count)
    tetra = np.asarray(mesh['tetrahedra'], int) - 1
    centers = X[tetra].mean(axis=1)
    owners, neighbours = np.asarray(mesh['owner'], int) - 1, np.asarray(mesh['neighbour'], int) - 1
    adjacency = [set() for _ in range(len(Q))]
    for a, b in zip(owners, neighbours):
        if b >= 0:
            adjacency[a].add(b)
            adjacency[b].add(a)
    values = primitive(Q, gas)
    values[:, 4] *= values[:, 0] * gas['R']
    recovered = np.empty((count, 5))
    factors = np.ones(count)
    maximum_fit_error = 0.
    saved = actual['nodalTrace']
    scales = np.r_[study['upstream']['density'], np.full(3, np.linalg.norm(study['upstream']['velocityVector'])),
                   study['upstream']['density'] * np.linalg.norm(study['upstream']['velocityVector']) ** 2]
    for j, node in enumerate(node_ids):
        visited = set(np.flatnonzero(np.any(tetra == node, axis=1)).tolist())
        front = set(visited)
        for ring in range(4):
            ids = np.array(sorted(visited), int)
            if ring >= 1 and len(ids) >= 6:
                dx = centers[ids] - X[node]
                distance = np.linalg.norm(dx, axis=1)
                h = distance.max()
                A = np.column_stack((np.ones(len(ids)), dx / h)) / distance[:, None]
                if np.linalg.svd(A, compute_uv=False)[-1] > 1e-10 * np.linalg.svd(A, compute_uv=False)[0]:
                    anchor = values[ids].mean(axis=0)
                    if np.all(values[ids] == values[ids[0]]):
                        anchor = values[ids[0]].copy()
                    fit, _, rank, _ = lstsq(A, (values[ids] - anchor) / distance[:, None],
                                          cond=1e-12, lapack_driver='gelsy')
                    assert rank == 4
                    break
            next_ids = set().union(*(adjacency[k] for k in front)).difference(visited)
            visited.update(next_ids)
            front = next_ids
        else:
            raise AssertionError('独立流体侧节点邻域秩不足。')
        np.testing.assert_array_equal(ids, np.asarray(saved['stencilCellIDs'][j], int) - 1)
        assert ring == int(saved['stencilRings'][j])
        unlimited = anchor + fit[0]
        lower, upper = gas['R'] * np.asarray(gas['temperatureRangeK'])
        base = np.array([anchor[0], anchor[4], anchor[4] - lower * anchor[0], upper * anchor[0] - anchor[4]])
        trial = np.array([unlimited[0], unlimited[4], unlimited[4] - lower * unlimited[0], upper * unlimited[0] - unlimited[4]])
        bad = trial < 0
        bad[:2] = trial[:2] <= 0
        if np.any(bad):
            factors[j] = max(0., min(1., np.min(base[bad] / (base[bad] - trial[bad])))) * (1 - 64 * np.finfo(float).eps)
        recovered[j] = anchor + factors[j] * fit[0]
        maximum_fit_error = max(maximum_fit_error, float(np.max(abs(recovered[j] - saved['primitiveRP'][j]) / scales)))
    assert maximum_fit_error < 2e-10
    np.testing.assert_allclose(factors, saved['slopeFactor'], rtol=2e-11, atol=2e-12)
    boundary = np.asarray(mesh['boundaryFaces'], int) - 1
    outside = boundary[np.asarray(shell['outerBoundaryMask'], bool)]
    faces = np.asarray(mesh['faces'], int)[outside] - 1
    xyz = X[faces]
    vectors = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]) / 2
    area = np.linalg.norm(vectors, axis=1)
    normals = vectors / area[:, None]
    node_vectors = np.zeros((count, 3))
    for column in range(3):
        np.add.at(node_vectors, faces[:, column] - layer * count, vectors)
    node_normals = node_vectors / np.linalg.norm(node_vectors, axis=1)[:, None]
    np.testing.assert_allclose(node_normals, actual['nodalNormal'], rtol=2e-11, atol=2e-13)
    gamma, up = gas['gamma'], study['upstream']
    sound_up = math.sqrt(gamma * gas['R'] * up['temperature'])
    speed = np.zeros(count)
    for j in range(count):
        un = float(np.asarray(up['velocityVector']) @ node_normals[j])
        sound_down = math.sqrt(gamma * recovered[j, 4] / recovered[j, 0])
        target = (recovered[j, 1:4] @ node_normals[j] + 2 * sound_down / (gamma - 1) - un) / sound_up
        def residual(M):
            density_ratio = (gamma + 1) * M * M / ((gamma - 1) * M * M + 2)
            pressure_ratio = (2 * gamma * M * M - gamma + 1) / (gamma + 1)
            return M * (1 - 1 / density_ratio) + 2 / (gamma - 1) * math.sqrt(pressure_ratio / density_ratio) - target
        high = 2.
        while residual(high) < 0:
            high *= 2
        speed[j] = un + sound_up * brentq(residual, 1., high, xtol=1e-13, rtol=1e-14)
    np.testing.assert_allclose(speed, actual['nodalNormalSpeedMS'], rtol=2e-11, atol=2e-8)
    rays = X[node_ids] - X[:count]
    rays /= np.linalg.norm(rays, axis=1)[:, None]
    radial = speed / np.sum(rays * node_normals, axis=1)
    outer_velocity = radial[:, None] * rays
    velocity = np.zeros_like(X)
    for k in range(1, layer + 1):
        velocity[k * count + np.arange(count)] = k / layer * outer_velocity
    np.testing.assert_allclose(velocity, actual['nodeVelocityMS'], rtol=2e-11, atol=2e-8)
    requested = np.array([shock_speed(row, n, up, gas)[0]
                          for row, n in zip(primitive(Q, gas)[owners[outside]], normals)])
    face_speed = velocity[faces].mean(axis=1)
    face_speed = np.sum(face_speed * normals, axis=1)
    residual = face_speed - requested
    np.testing.assert_allclose(residual, actual['projectionResidualMS'], rtol=2e-11, atol=2e-8)
    dt = study['history'][-1, 1]
    np.testing.assert_allclose(X + dt * velocity, study['nodes'], rtol=0, atol=2e-17)
    np.testing.assert_array_equal(study['nodes'][:count], study['initialNodes'][:count])
    return dict(method='independent-pivoted-QR-affine-trace-and-Brent-Jplus-root',
                maximumTraceErrorByPhysicalScale=maximum_fit_error,
                maximumNodeSpeedDifferenceMS=float(np.max(abs(velocity - actual['nodeVelocityMS']))),
                maximumUnrepresentedFaceSpeedMS=float(np.max(abs(residual))),
                limitedNodes=int(np.count_nonzero(factors < 1)),
                faceLeastSquaresOptimalityClaimed=False)


def flux_audit(study):
    """模块二：在真实网格速度处独立采样Riemann；保留RH五量缺陷及原单元压力差。"""
    mesh = study['cap']['shell']['mesh']
    gas, up, Q = study['gas'], study['upstream'], study['lastInputQ']
    dt = float(study['history'][-1, 1])
    geo = geometry(mesh, study['lastNodes0'], study['nodes'], dt)
    area = np.linalg.norm(geo['Abar'], axis=1)
    normals, D = geo['Abar'] / area[:, None], geo['swept'] / (dt * area)
    bound = np.asarray(mesh['boundaryFaces'], int) - 1
    rows = {int(face): j for j, face in enumerate(bound)}
    PL, right = primitive(Q, gas), np.r_[up['density'], up['velocityVector'], up['temperature']]
    U = np.asarray(up['velocityVector'])
    Qup = np.r_[up['density'], up['density'] * U,
                up['pressure'] / (gas['gamma'] - 1) + up['density'] * (U @ U) / 2]
    speed_scale = max(np.linalg.norm(U), math.sqrt(gas['gamma'] * gas['R'] * up['temperature']))
    physical_scale = up['density'] * np.r_[speed_scale, np.full(3, speed_scale ** 2), speed_scale ** 3]
    flux = np.zeros((len(geo['owner']), 5))
    pressures = np.full(len(bound), np.nan)
    wall = np.full(len(bound), np.nan)
    defects = np.full((len(bound), 5), np.nan)
    star = np.full(len(bound), np.nan)
    min_outlet = math.inf
    for face, cell in enumerate(geo['owner']):
        left, normal, speed = PL[cell], normals[face], D[face]
        neighbour = geo['neighbour'][face]
        if neighbour >= 0:
            flux[face] = exact_at(left, PL[neighbour], normal, speed, gas)[0]
            continue
        row, kind = rows[face], study['boundary']['kind'][rows[face]]
        if kind == 2:
            mirror = left.copy()
            mirror[1:4] -= 2 * (left[1:4] @ normal - speed) * normal
            wall[row] = exact_at(left, mirror, normal, speed, gas)[1]
            flux[face] = np.r_[0., wall[row] * normal, wall[row] * speed]
        elif kind in (3, 5):
            reference = euler_ale(Qup, normal, speed, gas['gamma'])
            relative = speed - U @ normal
            M2 = relative ** 2 / (gas['gamma'] * gas['R'] * up['temperature'])
            assert M2 > 1
            pressures[row] = up['pressure'] * (2 * gas['gamma'] * M2 - gas['gamma'] + 1) / (gas['gamma'] + 1)
            if kind == 5:
                flux[face], star[row] = exact_at(left, right, normal, speed, gas)
                assert star[row] > up['pressure']
            else:
                flux[face] = reference
            defects[row] = (flux[face] - reference) / physical_scale
        elif kind == 4:
            mach = (left[1:4] @ normal - speed) / math.sqrt(gas['gamma'] * gas['R'] * left[4])
            assert mach > 1
            min_outlet = min(min_outlet, mach)
            flux[face] = euler_ale(Q[cell], normal, speed, gas['gamma'])
        else:
            raise AssertionError('球鼻边界身份不符。')
    step = study['lastStep']
    scale = np.maximum(1., np.max(abs(flux), axis=0))
    error = float(np.max(abs(flux - step['faceALEFlux']) / scale))
    assert error < 2e-11
    np.testing.assert_allclose(step['wallPressurePa'], wall, rtol=2e-11, atol=0, equal_nan=True)
    np.testing.assert_allclose(step['shockPressurePa'], pressures, rtol=2e-11, atol=0, equal_nan=True)
    cell_pressure = PL[np.asarray(mesh['owner'], int)[bound] - 1, 0] * gas['R'] * PL[np.asarray(mesh['owner'], int)[bound] - 1, 4]
    mismatch = cell_pressure / pressures - 1
    np.testing.assert_allclose(step['shockInteriorPressureRelativeMismatch'], mismatch, rtol=2e-11, atol=2e-12, equal_nan=True)
    if 'shockBoundaryFluxRelativeMismatch' in step:
        np.testing.assert_allclose(step['shockBoundaryFluxRelativeMismatch'], defects, rtol=2e-10, atol=2e-12, equal_nan=True)
        np.testing.assert_allclose(step['shockRiemannStarPressurePa'], star, rtol=2e-11, atol=0, equal_nan=True)
    balance = np.zeros_like(Q)
    integrated = area[:, None] * flux
    np.add.at(balance, geo['owner'], integrated)
    inside = geo['neighbour'] >= 0
    np.add.at(balance, geo['neighbour'][inside], -integrated[inside])
    expected = (geo['V0'][:, None] * Q - dt * balance) / geo['V1'][:, None]
    qscale = up['density'] * np.r_[1., np.full(3, np.linalg.norm(U)), np.linalg.norm(U) ** 2]
    state_error = float(np.max(abs(expected - study['Q']) / qscale))
    assert state_error < 2e-11
    history = study['history']
    global_error = float(np.max(abs(np.diff(study['conservedIntegrals'], axis=0) +
                                  history[:, 1, None] * study['boundaryFluxRates']) /
                               (geo['V0'].sum() * qscale)))
    assert global_error < 2e-11
    return dict(fluxRelativeError=error, stateRelativeError=state_error,
                independentGCLRelative=geo['gcl'], allStepsBalanceByFreestreamScale=global_error,
                minimumOutletNormalMach=float(min_outlet),
                maximumIndependentRHFluxDefect=float(np.nanmax(abs(defects))),
                maximumIndependentCellPressureMismatch=float(np.nanmax(abs(mismatch))))


def current_wall_pressure(study, Q=None, nodes=None):
    """模块三：用当前保存帧状态派生当前壁压；避免把旧步通量采样值标成终态。"""
    Q = study['Q'] if Q is None else Q
    nodes = study['nodes'] if nodes is None else nodes
    shell, gas = study['cap']['shell'], study['gas']
    mesh = shell['mesh']
    boundary = np.asarray(mesh['boundaryFaces'], int) - 1
    walls = np.asarray(shell['wallBoundaryMask'], bool)
    faces = np.asarray(mesh['faces'], int)[boundary[walls]] - 1
    xyz = nodes[faces]
    vectors = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    normals = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    left = primitive(Q, gas)[np.asarray(mesh['owner'], int)[boundary[walls]] - 1]
    result = np.full(len(boundary), np.nan)
    for j, row in enumerate(np.flatnonzero(walls)):
        right = left[j].copy()
        right[1:4] -= 2 * (left[j, 1:4] @ normals[j]) * normals[j]
        result[row] = exact_at(left[j], right, normals[j], 0., gas)[1]
    return result


def audit(evidence, process_file, target):
    """模块四：绑定已终态实际任务，完整检查方法分支、温域、边界及未合格身份。"""
    assert not target.exists()
    process = json.loads(process_file.read_text(encoding='utf-8-sig'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
    saved = loadmat(evidence, simplify_cells=True)
    study = saved['study']
    for key in ('productionQualified', 'heatFluxSolved', 'wallTemperatureSolved', 'shockShapeSolved'):
        assert not study[key]
    assert study['bodyVerticesFixed'] and not study['processPeakMemoryMeasured']
    history = np.asarray(study['history'])
    assert history.shape[1] == 16 and np.all(history[:, 1] > 0)
    assert np.all(history[:, 6] >= study['gas']['temperatureRangeK'][0])
    assert np.all(history[:, 7] <= study['gas']['temperatureRangeK'][1])
    assert np.all(history[:, 8] > 0) and np.all(history[:, 12] > 1)
    assert study['maximumDeclaredRegionArrayBytes'] <= 128 * 1024 ** 2
    method = study.get('motionMethod', 'face-riemann-L2')
    motion = nodal_motion(study) if method == 'nodal-characteristic' else motion_audit(study)
    flux = flux_audit(study)
    steady = history[-1, 2] < 5e-4 and history[-1, 3] < 1e-4 and history[-1, 4] < 1e-4
    if method == 'nodal-characteristic':
        steady = steady and study['shockTracePressureMismatchHistory'][-1] < 1e-4 and np.max(abs(study['lastMotion']['nodalNormalSpeedMS'])) / np.linalg.norm(study['upstream']['velocityVector']) < 1e-4
    if study.get('shockFluxMethod', 'RH') == 'upstream-riemann':
        steady = steady and study['shockBoundaryFluxMismatchHistory'][-1] < 1e-4
    assert bool(study['candidateDiscreteSteady']) == bool(steady)
    current_wall = current_wall_pressure(study)
    old_wall = study['lastStep']['wallPressurePa']
    walls = np.isfinite(current_wall)
    result = dict(status='passed-independent-interface-transient-and-current-frame-wall-audit-not-region-qualification',
                  evidenceSHA256=sha(evidence), processSHA256=sha(process_file), motionMethod=method,
                  shockFluxMethod=study.get('shockFluxMethod', 'RH'),
                  acceptedSteps=len(history), timeS=float(study['timeS']),
                  lastStepInputTimeS=float(history[-1, 0] - history[-1, 1]),
                  bodyFacets=len(study['cap']['F']), tetrahedra=len(study['Q']),
                  radialLayers=int(study['cap']['shell']['layers']),
                  lastStep=flux, sharedNodeMotion=motion, finalResidual=float(history[-1, 2]),
                  finalShockPressureMismatch=float(history[-1, 15]),
                  currentWallPressureRangePa=[float(np.nanmin(current_wall)), float(np.nanmax(current_wall))],
                  lastStepWallVsCurrentMaximumRelativeDifference=float(np.max(abs(old_wall[walls] / current_wall[walls] - 1))),
                  currentWallMethod='当前保存五量/当前固定壁面法向/独立镜像Riemann压力根；没有重算流场',
                  maximumFacetCentroidSagM=float(study['cap']['maximumFacetCentroidSagM']),
                  declaredRegionArrayUpperBytes=float(study['maximumDeclaredRegionArrayBytes']),
                  processPeakMemoryMeasured=False, candidateDiscreteSteady=bool(steady),
                  productionQualified=False, newProductionSolvedAreaM2=0, experimentalValidation=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes((json.dumps(result, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('evidence', 'process', 'target'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    audit(args.evidence, args.process, args.target)
