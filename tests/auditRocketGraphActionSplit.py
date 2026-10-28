# -*- coding: utf-8 -*-
"""独立重建分项WLS/两点图作用、装配矩阵与误差相消；不求新流场。"""
from pathlib import Path
import csv, hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
case, kind = sys.argv[1:3]
assert case.replace('_', '').isalnum() and kind in ('unit', 'real')
folder = VAL / case
target = VAL / (case + '_independent_graph_split.json')
assert not target.exists()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


bindings = {}


def read(path, variables=None):
    path = Path(path)
    bindings[str(path.relative_to(ROOT))] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def rows(x):
    return [x] if isinstance(x, dict) else list(np.atleast_1d(x))


def rel(a, b):
    return float(np.linalg.norm(a - b) / max(1., np.linalg.norm(b)))


if kind == 'unit':
    context = read(folder / 'graph_split_unit_evidence.mat')
    items = rows(context['evidence'])
    operator = context['operator']
    source = read(Path(context['sourceFile']), ['unitA'])
    np.testing.assert_array_equal(source['unitA'], context['unitA'])
else:
    context = read(folder / 'graph_split_context.mat')
    source = read(Path(context['sourceFile']), ['result'])['result']
    parent = read(Path(context['parentFile']), ['result'])['result']
    np.testing.assert_array_equal(source['Q'], context['Q'])
    np.testing.assert_array_equal(source['Q'], parent['Q'])
    assert source['acceptedNonlinearIterations'] == parent['acceptedNonlinearIterations'] == 0
    operator = source['firstOperator']
    last = rows(source['linearLog'])[-1]['steadyDirectionEvidence']
    previous = rows(parent['linearLog'])[-1]['steadyDirectionEvidence']
    expected_directions = [last['directionU'], last['directionU'] - previous['directionU'], context['rhs'] + last['value']]
    expected_directions[2] /= np.max(abs(expected_directions[2]))
    np.testing.assert_array_equal(last['value'], context['savedOriginalNSAction'])
    del source, parent, last, previous
    items = None

Q, cut, gas, scale, base = [context[k] for k in ('Q', 'cut', 'gas', 'scale', 'baseInfo')]
f, w = cut['openFaces'], cut['walls']
v, centers = cut['cells']['volume'], cut['cells']['centroid']
n = len(v)
m = f['cellMinus'].astype(int) - 1
p = f['cellPlus'].astype(int) - 1
owner = w['ownerCell'].astype(int) - 1
inside, outside = np.flatnonzero(p >= 0), np.flatnonzero(p < 0)
assert np.all(m >= 0) and np.all(owner >= 0) and np.all(f['area'] > 0)
gamma, R, mu, pr = [gas[k] for k in ('gamma', 'R', 'mu', 'Pr')]
gm = gamma - 1
cv, kappa = R / gm, mu * gamma * R / gm / pr
rho = Q[:, 0]
u = Q[:, 1:4] / rho[:, None]
T = (Q[:, 4] - .5 * rho * np.sum(u * u, axis=1)) / (rho * cv)
P = np.column_stack((rho, u, T))
assert rel(P, base['primitive']) < 1e-14
nf, nl = int(operator['sourceFaces']), int(operator['layers'])
assert nf * nl == n
names = ['interiorNormalEuler', 'interiorSideEuler', 'outerEuler', 'solidWallEuler',
         'momentumStress', 'heatConduction', 'workGradient', 'workVelocity']


def resolve(P, b, normal):
    """采用直接幂关系重建给定储库的特征状态，不调用MATLAB缓存。"""
    density, velocity, temperature = P[0], P[1:4], P[4]
    un = float(velocity @ normal)
    acoustic = np.sqrt(gamma * R * temperature)
    if b['type'] == 'supersonicOutflow':
        return P.copy(), False
    if b['type'] == 'supersonicInflow':
        return np.asarray(b['primitive']), True
    assert b['type'] == 'characteristic'
    reservoir = np.asarray(b['reservoirPrimitive'])
    if un > acoustic:
        return P.copy(), False
    if un < -acoustic:
        return reservoir.copy(), True
    ar = np.sqrt(gamma * R * reservoir[4])
    ur = float(reservoir[1:4] @ normal)
    K, Kr = R * temperature / density**gm, R * reservoir[4] / reservoir[0]**gm
    jp, jm = un + 2 * acoustic / gm, ur - 2 * ar / gm
    ub, ab = (jp + jm) / 2, gm * (jp - jm) / 4
    incoming = ub < 0
    entropy, tangent = (Kr, reservoir[1:4] - ur * normal) if incoming else (K, velocity - un * normal)
    acoustic_state = np.r_[(ab * ab / (gamma * entropy))**(1 / gm), tangent + ub * normal, ab * ab / (gamma * R)]
    if b['mode'] == 'riemannFarfield':
        return acoustic_state, incoming
    pressure = float(b['pressure'])
    rb = (pressure / K)**(1 / gamma)
    tb = pressure / (rb * R)
    aa = np.sqrt(gamma * R * tb)
    uu = jp - 2 * aa / gm
    pressure_state = np.r_[rb, velocity - un * normal + uu * normal, tb]
    return (acoustic_state, True) if (un < 0 or uu < 0) and incoming else (pressure_state, False)


state = np.full((len(m), 5), np.nan)
incoming = np.zeros(len(m), dtype=bool)
boundary_names = [['xmin', 'xmax'], ['ymin', 'ymax'], ['zmin', 'zmax']]
for face in outside:
    b = context['boundary'][boundary_names[int(f['axis'][face]) - 1][1]]
    state[face], incoming[face] = resolve(P[m[face]], b, f['normal'][face])


def divergence(face, wall):
    face = np.asarray(face).reshape((len(m), -1))
    wall = np.asarray(wall).reshape((len(owner), -1))
    out = np.zeros((n, face.shape[1]))
    np.add.at(out, m, f['area'][:, None] / v[m, None] * face)
    np.add.at(out, p[inside], -f['area'][inside, None] / v[p[inside], None] * face[inside])
    np.add.at(out, owner, w['area'][:, None] / v[owner, None] * wall)
    return out


def graph_matrix_action(direction):
    out = np.einsum('abn,nb->na', operator['diagonal'], direction)
    out[nf:] += np.einsum('abn,nb->na', operator['lower'][:, :, nf:], direction[:-nf])
    out[:-nf] += np.einsum('abn,nb->na', operator['upper'][:, :, :-nf], direction[nf:])
    mi = operator['sideMinus'].astype(int) - 1
    pi = operator['sidePlus'].astype(int) - 1
    for first in range(0, len(mi), 8192):
        ids = slice(first, first + 8192)
        np.add.at(out, mi[ids], np.einsum('abn,nb->na', operator['minusCoupling'][:, :, ids], direction[pi[ids]]))
        np.add.at(out, pi[ids], np.einsum('abn,nb->na', operator['plusCoupling'][:, :, ids], direction[mi[ids]]))
    return direction / operator['pseudoTimeStepSeconds'] - out


def force(G, normal):
    velocity_gradient = G[:, :3, :]
    stress = mu * (velocity_gradient + velocity_gradient.swapaxes(1, 2) -
                   (2 / 3) * np.trace(velocity_gradient, axis1=1, axis2=2)[:, None, None] * np.eye(3))
    return np.einsum('nij,nj->ni', stress, normal)


# WLS几何正规矩阵只装配一次，每个方向另行装配线性右端。
normal_matrix = np.zeros((n, 3, 3))
for cells, delta, area in [(m[inside], centers[p[inside]] - centers[m[inside]], f['area'][inside]),
                           (p[inside], centers[p[inside]] - centers[m[inside]], f['area'][inside]),
                           (m[outside[incoming[outside]]], f['centroid'][outside[incoming[outside]]] - centers[m[outside[incoming[outside]]]], f['area'][outside[incoming[outside]]]),
                           (owner, w['centroid'] - centers[owner], w['area'])]:
    weight = area / np.sum(delta * delta, axis=1)
    np.add.at(normal_matrix, cells, weight[:, None, None] * delta[:, :, None] * delta[:, None, :])


def audit(item, ordinal):
    direction = item['directionU']
    if kind == 'real':
        np.testing.assert_array_equal(direction, expected_directions[ordinal])
    detail, full_saved, graph_saved = item['detail'], item['full'], item['graph']
    dq = direction * scale
    du = (dq[:, 1:4] - u * dq[:, 0, None]) / rho[:, None]
    dT = (dq[:, 4] - np.sum(u * dq[:, 1:4], axis=1) +
          (.5 * np.sum(u * u, axis=1) - cv * T) * dq[:, 0]) / (rho * cv)
    dphi = np.column_stack((du, dT))
    dGiven = detail['boundaryDirection']['primitiveDerivative']
    np.testing.assert_array_equal(incoming, detail['boundaryDirection']['incoming'].astype(bool))
    assert not np.any(detail['boundaryDirection']['wall'])
    B = np.zeros((n, 3, 4))
    in_given = outside[incoming[outside]]
    for cells, delta, difference, area in [
            (m[inside], centers[p[inside]] - centers[m[inside]], dphi[p[inside]] - dphi[m[inside]], f['area'][inside]),
            (p[inside], centers[p[inside]] - centers[m[inside]], dphi[p[inside]] - dphi[m[inside]], f['area'][inside]),
            (m[in_given], f['centroid'][in_given] - centers[m[in_given]], dGiven[in_given, 1:5] - dphi[m[in_given]], f['area'][in_given]),
            (owner, w['centroid'] - centers[owner], -dphi[owner], w['area'])]:
        weight = area / np.sum(delta * delta, axis=1)
        np.add.at(B, cells, weight[:, None, None] * delta[:, :, None] * difference[:, None, :])
    gradient = np.linalg.solve(normal_matrix, B).swapaxes(1, 2)
    vd, ed = detail['viscousDetail'], detail['eulerDetail']
    gradient_error = rel(gradient, np.asarray(vd['gradientDerivative']).transpose(2, 0, 1))
    assert gradient_error < 2e-9
    stress = np.zeros((len(m), 3))
    heat, wg, wv = [np.zeros(len(m)) for _ in range(3)]
    graph_stress = np.zeros_like(stress)
    graph_heat, graph_work = np.zeros(len(m)), np.zeros(len(m))
    for first in range(0, len(m), 8192):
        ids = np.arange(first, min(first + 8192, len(m)))
        a, b = m[ids], p[ids]
        inmask = b >= 0
        G = gradient[a].copy()
        velocity, dvelocity = u[a].copy(), du[a].copy()
        inside_ids, outside_ids = np.flatnonzero(inmask), np.flatnonzero(~inmask)
        mi, pi = a[inside_ids], b[inside_ids]
        delta = centers[pi] - centers[mi]
        G[inside_ids] = .5 * (gradient[mi] + gradient[pi])
        G[inside_ids] += (dphi[pi] - dphi[mi] - np.einsum('nij,nj->ni', G[inside_ids], delta))[:, :, None] * delta[:, None, :] / np.sum(delta * delta, axis=1)[:, None, None]
        velocity[inside_ids], dvelocity[inside_ids] = .5 * (u[mi] + u[pi]), .5 * (du[mi] + du[pi])
        ff = ids[outside_ids]
        ii = outside_ids[incoming[ff]]
        faces, cells = ids[ii], a[ii]
        delta = f['centroid'][faces] - centers[cells]
        G[ii] += (dGiven[faces, 1:5] - dphi[cells] - np.einsum('nij,nj->ni', G[ii], delta))[:, :, None] * delta[:, None, :] / np.sum(delta * delta, axis=1)[:, None, None]
        ii = outside_ids[~incoming[ff]]
        normal = f['normal'][ids[ii]]
        G[ii] -= np.einsum('nij,nj->ni', G[ii], normal)[:, :, None] * normal[:, None, :]
        velocity[outside_ids], dvelocity[outside_ids] = state[ff, 1:4], dGiven[ff, 1:4]
        normal = f['normal'][ids]
        stress[ids] = force(G, normal)
        heat[ids] = kappa * np.sum(G[:, 3, :] * normal, axis=1)
        wg[ids] = np.sum(velocity * stress[ids], axis=1)
        wv[ids] = np.sum(dvelocity * base['viscousFaceFlux'][ids, 1:4], axis=1)
        dv, dt, frozen_velocity = -du[a], -dT[a], u[a].copy()
        delta = f['centroid'][ids] - centers[a]
        delta[inmask] = centers[b[inmask]] - centers[a[inmask]]
        dn = np.sum(delta * normal, axis=1)
        assert np.all(dn > 0)
        dv[inmask], dt[inmask] = du[b[inmask]] - du[a[inmask]], dT[b[inmask]] - dT[a[inmask]]
        frozen_velocity[inmask] = .5 * (u[a[inmask]] + u[b[inmask]])
        graph_stress[ids] = mu / dn[:, None] * (dv + np.sum(dv * normal, axis=1)[:, None] * normal / 3)
        graph_heat[ids], graph_work[ids] = kappa * dt / dn, np.sum(frozen_velocity * graph_stress[ids], axis=1)
    dn = np.sum((w['centroid'] - centers[owner]) * w['normal'], axis=1)
    wall_stress = -mu / dn[:, None] * (du[owner] + np.sum(du[owner] * w['normal'], axis=1)[:, None] * w['normal'] / 3)
    wall_heat = -kappa * dT[owner] / dn
    full, graph = {}, {}
    flux, wall = ed['faceEulerFluxDerivative'], ed['wallEulerFluxDerivative']
    for name, mask in [('interiorNormalEuler', (p >= 0) & (f['kind'] == 1)),
                       ('interiorSideEuler', (p >= 0) & (f['kind'] == 3)), ('outerEuler', p < 0)]:
        selected = np.zeros_like(flux)
        selected[mask] = flux[mask]
        full[name] = -divergence(selected, np.zeros_like(wall)) / scale
        graph[name] = full[name].copy()
    full['solidWallEuler'] = -divergence(np.zeros_like(flux), wall) / scale
    graph['solidWallEuler'] = full['solidWallEuler'].copy()
    wall_pressure = gm * (dq[owner, 4] - np.sum(u[owner] * dq[owner, 1:4], axis=1) + .5 * np.sum(u[owner]**2, axis=1) * dq[owner, 0])
    assert rel(wall[:, 1:4], wall_pressure[:, None] * w['normal']) < 2e-13
    # 独立Euler链式项重建外边界完整通量，以及冻结储库的图通量。
    def euler_action(primitive, dprimitive, normal):
        density, velocity, temperature = primitive[:, 0], primitive[:, 1:4], primitive[:, 4]
        dr, dv, dt = dprimitive[:, 0], dprimitive[:, 1:4], dprimitive[:, 4]
        pressure = density * R * temperature
        dp = R * (temperature * dr + density * dt)
        energy = density * (cv * temperature + .5 * np.sum(velocity * velocity, axis=1))
        de = dr * (cv * temperature + .5 * np.sum(velocity * velocity, axis=1)) + density * (cv * dt + np.sum(velocity * dv, axis=1))
        un, dun = np.sum(velocity * normal, axis=1), np.sum(dv * normal, axis=1)
        momentum = density[:, None] * velocity
        dm = dr[:, None] * velocity + density[:, None] * dv
        return np.column_stack((dr * un + density * dun,
                                dm * un[:, None] + momentum * dun[:, None] + dp[:, None] * normal,
                                (de + dp) * un + (energy + pressure) * dun))
    outer_error = rel(euler_action(state[outside], dGiven[outside], f['normal'][outside]), flux[outside])
    assert outer_error < 2e-10
    po = P[m[outside]]
    do = np.column_stack((dq[m[outside], 0], du[m[outside]], dT[m[outside]]))
    wave = abs(np.sum(po[:, 1:4] * f['normal'][outside], axis=1)) + np.sqrt(gamma * R * po[:, 4])
    selected = np.zeros_like(flux)
    selected[outside] = .5 * euler_action(po, do, f['normal'][outside]) + .5 * wave[:, None] * dq[m[outside]]
    graph['outerEuler'] = -divergence(selected, np.zeros_like(wall)) / scale
    for name, ff, fw, gf, gw, columns in [
            ('momentumStress', stress, wall_stress, graph_stress, wall_stress, slice(1, 4)),
            ('heatConduction', heat, wall_heat, graph_heat, wall_heat, slice(4, 5)),
            ('workGradient', wg, np.zeros(len(owner)), graph_work, np.zeros(len(owner)), slice(4, 5)),
            ('workVelocity', wv, np.zeros(len(owner)), np.zeros(len(m)), np.zeros(len(owner)), slice(4, 5))]:
        full[name], graph[name] = np.zeros_like(Q), np.zeros_like(Q)
        full[name][:, columns] = divergence(ff, fw)
        graph[name][:, columns] = divergence(gf, gw)
        full[name] /= scale
        graph[name] /= scale
    full_error = max(rel(full[name], full_saved[name]) for name in names)
    graph_error = max(rel(graph[name], graph_saved[name]) for name in names)
    assert max(full_error, graph_error) < 2e-9
    full_total, graph_total = sum(full.values()), sum(graph.values())
    assert max(rel(full_total, detail['fullTotal']), rel(graph_total, detail['graphTotal'])) < 2e-9
    assembled = graph_matrix_action(direction)
    matrix_error = rel(graph_total, assembled)
    assert matrix_error < 2e-10
    difference = graph_total - full_total
    errors = np.column_stack([(graph[name] - full[name]).ravel() for name in names])
    error2 = np.sum(difference**2)
    gram = errors.T @ errors / error2
    signed = gram.sum(axis=1)
    summary = item.get('summary')
    if summary is not None:
        np.testing.assert_allclose(summary['componentDifferenceGram'], gram, rtol=3e-9, atol=2e-11)
        values = np.array([[np.linalg.norm(full[name]), np.linalg.norm(graph[name]), np.linalg.norm(graph[name] - full[name]),
                            signed[j], np.linalg.norm(difference - (graph[name] - full[name])) / np.linalg.norm(full_total)]
                           for j, name in enumerate(names)])
        np.testing.assert_allclose(summary['componentRows'], values, rtol=2e-9, atol=2e-10)
        layer = np.sum(difference.reshape((nl, nf, 5))**2, axis=(1, 2)) / error2
        np.testing.assert_allclose(summary['layerDifferenceSquaredNormFraction'], layer, rtol=2e-9, atol=1e-12)
        order = np.argsort(np.sum(difference**2, axis=1))[::-1][:20]
        np.testing.assert_array_equal(summary['largestDifferenceCellIDs'], order + 1)
        np.testing.assert_array_equal(summary['largestDifferenceOriginalFaceIDs'], order % nf + 1)
    else:
        measured = (item['plus'] - item['minus']) / (2 * context['h']) / scale
        assert rel(full_total, measured) < 2e-6
        assert rel(full_total, (-context['unitA'] @ direction.ravel(order='F')).reshape(Q.shape, order='F')) < 2e-13
    return dict(name=summary['name'] if summary else 'unit' + str(ordinal + 1),
                independentWlsGradientRelativeDifference=gradient_error, independentFullComponentRelativeDifference=full_error,
                independentGraphComponentRelativeDifference=graph_error, independentAssembledGraphRelativeDifference=matrix_error,
                independentOuterEulerRelativeDifference=outer_error,
                graphActionRelativeDifference=float(np.linalg.norm(difference) / np.linalg.norm(full_total)),
                componentNames=names, componentSignedDifferenceContributions=signed.tolist(),
                componentDifferenceNormRelativeToFull=[float(np.linalg.norm(graph[name] - full[name]) / np.linalg.norm(full_total)) for name in names],
                layerDifferenceSquaredNormFraction=(np.sum(difference.reshape((nl, nf, 5))**2, axis=(1, 2)) / error2).tolist())


checks = []
if kind == 'unit':
    for ordinal, item in enumerate(items):
        checks.append(audit(item, ordinal))
else:
    for ordinal, name in enumerate(('totalNewton', 'lastRecycledIncrement', 'normalizedRemaining')):
        item = read(folder / ('graph_split_' + name + '.mat'))
        checks.append(audit(item, ordinal))
        del item
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure = folder / 'graph_split_components.png'
    assert not figure.exists()
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), layout='constrained')
    picked = [2, 4, 5, 6, 7]
    labels = ['Outer Euler', 'WLS stress', 'Heat conduction', 'Work: U d(stress)', 'Work: dU stress']
    colors = ['#426a8c', '#ad6f3b', '#658569']
    for j, check in enumerate(checks):
        axes[0].bar(np.arange(len(picked)) + (j - 1) * .24,
                    np.asarray(check['componentDifferenceNormRelativeToFull'])[picked], width=.24, color=colors[j], label=check['name'])
        axes[1].plot(np.arange(1, nl + 1), 100 * np.asarray(check['layerDifferenceSquaredNormFraction']), 'o-', color=colors[j], ms=4, label=check['name'])
    axes[0].set(yscale='log', xticks=np.arange(len(picked)), xticklabels=labels,
                ylabel='Component difference norm / full Jd norm', title='Component norms include cancellation; they are not additive area shares')
    axes[0].tick_params(axis='x', labelsize=9)
    axes[0].legend(fontsize=9, ncol=3)
    axes[1].set(xlabel='Wall-normal layer', ylabel='Share of total squared action difference (%)', xticks=np.arange(1, nl + 1))
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
    fig.suptitle('M5 / 30 km / alpha 20 deg / UNCONVERGED RESEARCH\nCurrent J versus approximate Jgraph on three directions; not pressure or heat-flux errors', fontsize=11)
    fig.savefig(figure, dpi=160)
    plt.close(fig)
    with (folder / 'graph_split_components.csv').open('w', newline='', encoding='utf8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['direction', 'component', 'difference_norm_divided_by_full_Jd_norm', 'signed_contribution_to_total_squared_difference'])
        for check in checks:
            for name, norm, signed in zip(names, check['componentDifferenceNormRelativeToFull'], check['componentSignedDifferenceContributions']):
                writer.writerow([check['name'], name, norm, signed])

assert all(sha(ROOT / path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-current-J-graph-component-WLS-two-point-and-assembled-action-audit', case=case, kind=kind,
           processID=os.getpid(), directions=checks, newActualNSEvaluations=0, newFlowAcceptances=0,
           originalNSEvaluationsInBoundMatlabUnitTest=5 if kind == 'unit' else 0,
           independentFullNSSolver=False, fullJacobianConditioningClaimed=False, physicalSolution=False, newProductionSolvedAreaM2=0,
           inputSHA256=bindings, scriptSHA256=sha(__file__))
if kind == 'real':
    out['figure'] = dict(path=str(figure.relative_to(ROOT)), sha256=sha(figure))
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], case=case, processID=os.getpid(),
                      relativeActionDifferences=[x['graphActionRelativeDifference'] for x in checks])))
