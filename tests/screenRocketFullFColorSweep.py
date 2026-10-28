# -*- coding: utf-8 -*-
"""诊断完整 F 的跨组连接，并筛查固定颜色顺序的五量块乘性扫描；不推进流场。"""
from pathlib import Path
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import bsr_matrix, coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
FOLDER = VAL / 'fullFcolorscreen01'
assert not FOLDER.exists()
FOLDER.mkdir()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(part)
    return h.hexdigest()


def read(path, variables):
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def buckets(factor):
    return [factor['buckets']] if isinstance(factor['buckets'], dict) else list(np.atleast_1d(factor['buckets']).ravel())


def arrays(bucket):
    ng = np.size(bucket['groupIDs'])
    k = np.asarray(bucket['members']).size // ng
    nodes = np.asarray(bucket['members']).reshape(k, ng, order='F').astype(np.int64).T - 1
    inverse = np.asarray(bucket['inverse']).reshape(5*k, 5*k, ng, order='F').transpose(2, 0, 1)
    return nodes, inverse


def matrix(graph):
    n = int(graph['cellCount'])
    order = np.atleast_1d(graph['order']).astype(np.int64) - 1
    np.testing.assert_array_equal(order, np.arange(n))
    ptr = np.atleast_1d(graph['rowPointer']).astype(np.int64) - 1
    col = np.atleast_1d(graph['column']).astype(np.int64) - 1
    values = np.asarray(graph['values']).reshape(5, 5, len(col), order='F').transpose(2, 0, 1)
    return bsr_matrix((values, col, ptr), shape=(5*n, 5*n)), ptr, col


def partition(factor):
    n = int(factor['cellCount'])
    group = np.full(n, -1, np.int64)
    ordered_nodes, ordered_inverses = [], []
    for bucket in buckets(factor):
        if not np.size(bucket['groupIDs']):
            continue
        nodes, inverse = arrays(bucket)
        for page, ids in enumerate(nodes):
            assert np.all(group[ids] == -1)
            group[ids] = len(ordered_nodes)
            ordered_nodes.append(ids)
            ordered_inverses.append(inverse[page])
    assert np.all(group >= 0)
    return group, ordered_nodes, ordered_inverses


def color_groups(group, ptr, col):
    owner = np.repeat(np.arange(len(group)), np.diff(ptr))
    row_group, col_group = group[owner], group[col]
    cross = row_group != col_group
    ng = int(group.max()) + 1
    directed = coo_matrix((np.ones(int(np.count_nonzero(cross)), np.int8),
        (row_group[cross], col_group[cross])), shape=(ng, ng)).tocsr()
    directed.data[:] = 1
    undirected = directed.maximum(directed.T)
    degree = np.diff(undirected.indptr)
    order = np.lexsort((np.arange(ng), -degree))
    color = np.full(ng, -1, np.int64)
    for g in order:
        neighbours = undirected.indices[undirected.indptr[g]:undirected.indptr[g+1]]
        used = set(color[neighbours].tolist())
        c = 0
        while c in used:
            c += 1
        color[g] = c
    assert np.all(color[row_group[cross]] != color[col_group[cross]])
    nscc, component = connected_components(directed, directed=True, connection='strong')
    largest = int(np.bincount(component).max())
    return color, cross, dict(groups=ng, crossGroupOriginalBlocks=int(cross.sum()),
        colorCount=int(color.max())+1, maximumGroupNeighbourDegree=int(degree.max()),
        strongComponents=int(nscc), largestStrongComponentGroups=largest,
        allOriginalCrossGroupEdgesConnectDifferentColors=True,
        groupOrder='descending-exact-group-degree-then-frozen-factor-group-index')


def schedule(factor, group, color):
    per_color = [[] for _ in range(int(color.max())+1)]
    for bucket in buckets(factor):
        if not np.size(bucket['groupIDs']):
            continue
        nodes, inverse = arrays(bucket)
        colors = color[group[nodes[:, 0]]]
        assert np.all(group[nodes] == group[nodes[:, :1]])
        for c in np.unique(colors):
            pages = np.flatnonzero(colors == c)
            per_color[c].append((nodes[pages], inverse[pages]))
    return per_color


def sweep(A, per_color, rhs, reverse=False, symmetric=False):
    """按颜色解原组局部五量块，颜色间立即更新完整原 F 剩余。"""
    x = np.zeros_like(rhs)
    history = []
    order = list(range(len(per_color)))
    if reverse:
        order.reverse()
    if symmetric:
        order += order[::-1]
    denominator = max(np.finfo(float).tiny, np.linalg.norm(rhs))
    for c in order:
        r = rhs - (A @ x.ravel()).reshape(rhs.shape)
        for nodes, inverse in per_color[c]:
            ng, k = nodes.shape
            update = np.einsum('pij,pj->pi', inverse, r[nodes].reshape(ng, 5*k)).reshape(ng, k, 5)
            x[nodes] += update
        assert np.all(np.isfinite(x))
        history.append(float(np.linalg.norm(x)))
    image = (A @ x.ravel()).reshape(rhs.shape)
    return x, image, dict(residualRecomputations=len(order),
        finalImageEvaluations=1, maximumPartialDirectionNorm=max(history),
        relativeResidual=float(np.linalg.norm(rhs-image)/denominator))


def optimal_step(rhs, image):
    bn, sn = np.linalg.norm(rhs), np.linalg.norm(image)
    alpha = float(np.sum((rhs/bn)*(image/sn))) * bn/sn
    return dict(coefficient=alpha, relativeResidual=float(np.linalg.norm(rhs-alpha*image)/bn))


small = read(VAL / 'overlapFunit01/restricted_overlap_unit_evidence.mat', ['cases'])
small_reports = []
for case in small['cases']:
    factor, graph = case['plainFactor'], case['graph']
    A, ptr, col = matrix(graph)
    group, nodes, inverses = partition(factor)
    colors, cross, stats = color_groups(group, ptr, col)
    plan = schedule(factor, group, colors)
    b = case['rhs']
    dense = A.toarray()
    flat_color = np.repeat(colors[group], 5)
    # 独立稠密块三角参考：同颜色组之间没有原连接，故同色可并行。
    lower = np.where(flat_color[:, None] >= flat_color[None, :], dense, 0.)
    upper = np.where(flat_color[:, None] <= flat_color[None, :], dense, 0.)
    forward = np.linalg.solve(lower, b.ravel())
    reference = forward + np.linalg.solve(upper, b.ravel() - dense @ forward)
    x, image, detail = sweep(A, plan, b, symmetric=True)
    error = float(np.linalg.norm(x.ravel()-reference)/max(1., np.linalg.norm(reference)))
    assert error < 1e-12
    other = np.sin(np.arange(b.size)).reshape(b.shape)
    combined = sweep(A, plan, .7*b-.3*other, symmetric=True)[0]
    linear_error = float(np.linalg.norm(combined-(.7*x-.3*sweep(A, plan, other, symmetric=True)[0]))/max(1., np.linalg.norm(combined)))
    assert linear_error < 1e-12
    # 零右端直接沿相同固定映射验证，不改变扫描次数。
    zero, zero_image, zero_detail = sweep(A, plan, np.zeros_like(b), symmetric=True)
    assert np.all(zero == 0) and np.all(zero_image == 0) and zero_detail['relativeResidual'] == 0
    small_reports.append(dict(denseTriangularRelativeDifference=error, linearityRelativeDifference=linear_error,
        colors=stats['colorCount'], fullFRelativeResidual=detail['relativeResidual'], zeroRHSExactlyZero=True))
del small, case, A, dense, lower, upper, reference, factor, graph, plan
gc.collect()
print('SMALL_MULTICOLOR_TRIANGULAR_REFERENCES_PASSED', flush=True)

base_file = VAL / 'fullFreal01/full_F_real_evidence.mat'
base = read(base_file, ['assemblyFile', 'rhsF', 'Fimage', 'sourceBinding'])
assembly_file = Path(base['assemblyFile'])
audit = json.loads((VAL / 'spatialreductionreal01_independent_spatial_reduction.json').read_text())
assert sha(assembly_file) == audit['inputSHA256'][str(assembly_file.relative_to(ROOT))]
raw = read(assembly_file, ['reduction'])
factor = raw['reduction']['Ffactor']
del raw
gc.collect()
graph_file = VAL / 'fullFreal01/full_F_graph_evidence.mat'
graph_audit = json.loads((VAL / 'fullFreal01_independent_full_F.json').read_text())
assert sha(graph_file) == graph_audit['inputSHA256'][str(graph_file.relative_to(ROOT))]
graph = read(graph_file, ['graph'])['graph']
A, ptr, col = matrix(graph)
assert A.data.nbytes + A.indices.nbytes + A.indptr.nbytes < 384*1024**2
group, nodes, inverses = partition(factor)
colors, cross, stats = color_groups(group, ptr, col)
assert stats['groups'] == 28902 and stats['crossGroupOriginalBlocks'] == 682056
plan = schedule(factor, group, colors)
physical = np.atleast_1d(graph['originalPhysicalCVs']).astype(np.int64) - 1
remaining = base['rhsF'] - base['Fimage']
energy = np.sum(remaining**2, axis=1)
ranking = np.argsort(-energy, kind='stable')
top = ranking[:64]
layer = physical // 19752
statistics = dict(bestSavedFullFResidual=float(np.linalg.norm(remaining)/np.linalg.norm(base['rhsF'])),
    remainingLayerSquaredNormFractions=(np.bincount(layer, weights=energy, minlength=11)/energy.sum()).tolist(),
    remainingEquationSquaredNormFractions=(np.sum(remaining**2, axis=0)/energy.sum()).tolist(),
    top64OriginalCVsOneBased=(physical[top]+1).tolist(),
    top64OriginalFacesOneBased=(physical[top] % 19752+1).tolist(),
    top64LayersOneBased=(layer[top]+1).tolist(),
    top64SquaredNormFraction=float(energy[top].sum()/energy.sum()),
    top64DistinctOriginalGroups=int(len(np.unique(group[top]))))
topology_file = FOLDER / 'full_F_color_topology.mat'
savemat(topology_file, dict(groupIndexOneBased=group+1, groupColorOneBased=colors+1,
    graphEvidenceSHA256=sha(graph_file), originalFactorAssemblySHA256=sha(assembly_file),
    method='fixed-forward-and-reverse-multiplicative-original-F-block-sweep'), do_compression=True)
print(json.dumps(dict(groupGraph=stats, remaining=statistics)), flush=True)
comparisons = []
for label, rhs in [('original-F-RHS', base['rhsF']), ('best-old-32-step-F-residual', remaining)]:
    x, image, detail = sweep(A, plan, rhs, symmetric=True)
    detail['rhs'] = label
    detail['optimalSingleStep'] = optimal_step(rhs, image)
    comparisons.append(detail)
    print(json.dumps(detail), flush=True)
    if label == 'original-F-RHS':
        savemat(FOLDER / 'full_F_color_screen_vectors.mat', dict(rhsF=rhs,
            symmetricDirection=x, symmetricImage=image, sourceBinding=base['sourceBinding']), do_compression=True)
assert all(sha(ROOT / name) == digest for name, digest in bindings.items())
out = dict(status='tested-independent-dense-multicolor-F-symmetric-block-sweep-and-real-action-screening',
    processID=os.getpid(), recordedUTC=datetime.now(timezone.utc).isoformat(), smallReferences=small_reports,
    groupGraph=stats, remaining=statistics, actionComparisons=comparisons,
    topologyFile=topology_file.relative_to(ROOT).as_posix(), topologySHA256=sha(topology_file),
    fixedMap=True, allFiveConservativeFieldsRetained=True, newLocalInverseConstructions=0,
    newKrylovSolves=0, originalCOrSchurIdentityClaimed=False, physicalStateChanged=False,
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    inputSHA256=bindings, processPeakMemoryClaimed=False, scriptSHA256=sha(__file__))
(FOLDER / 'full_F_color_screen_summary.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], summary=str(FOLDER / 'full_F_color_screen_summary.json'))))
