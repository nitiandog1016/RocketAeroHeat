# -*- coding: utf-8 -*-
"""复用已认证扩展逆，比较核心输出与成员归一输出；只诊断原 F 方程，不更新流场。"""
from pathlib import Path
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
TARGET = VAL / 'overlap_output_diagnostic01.json'
assert not TARGET.exists()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read_mat(path, variables):
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def buckets(factor):
    value = factor['buckets']
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def bucket_arrays(bucket):
    count = np.size(bucket['groupIDs'])
    k, m = int(bucket['overlapSize']), int(bucket['coreSize'])
    members = np.asarray(bucket['members']).reshape(k, count, order='F').astype(np.int64) - 1
    core = np.asarray(bucket['coreMembers']).reshape(m, count, order='F').astype(np.int64) - 1
    inverse = np.asarray(bucket['inverse']).reshape(5*k, 5*k, count, order='F').transpose(2, 0, 1)
    return count, k, m, members, core, inverse


def maps(factor, rhs):
    """模块一：同一完整局部解分别散回核心和全部成员，再按成员次数归一。"""
    core_value, weighted = np.zeros_like(rhs), np.zeros_like(rhs)
    membership = np.zeros(len(rhs), np.int64)
    for bucket in buckets(factor):
        if not np.size(bucket['groupIDs']):
            continue
        count, k, m, members, core, inverse = bucket_arrays(bucket)
        solved = np.einsum('pij,pj->pi', inverse, rhs[members.T].reshape(count, 5*k))
        core_value[core.T] = solved[:, :5*m].reshape(count, m, 5)
        np.add.at(weighted, members.T.ravel(), solved.reshape(count*k, 5))
        np.add.at(membership, members.T.ravel(), 1)
    assert np.all(membership > 0)
    weighted /= membership[:, None]
    assert np.all(np.isfinite(core_value)) and np.all(np.isfinite(weighted))
    return core_value, weighted, membership


def norm(value):
    scale = float(np.max(abs(value)))
    return float(np.linalg.norm(value / scale) * scale) if scale else 0.


def rel(actual, reference):
    return norm(actual - reference) / max(1., norm(reference))


def dense_weighted(factor):
    """模块二：逐组稠密矩阵散加参考，独立于批量右端映射。"""
    n = int(factor['cellCount'])
    matrix, counts = np.zeros((5*n, 5*n)), np.zeros(n, np.int64)
    for bucket in buckets(factor):
        if not np.size(bucket['groupIDs']):
            continue
        count, k, m, members, core, inverse = bucket_arrays(bucket)
        for page in range(count):
            nodes = members[:, page]
            dofs = (5*nodes[:, None] + np.arange(5)[None, :]).ravel()
            matrix[np.ix_(dofs, dofs)] += inverse[page]
            counts[nodes] += 1
    matrix /= np.repeat(counts, 5)[:, None]
    return matrix, counts


def action(graph, value):
    """模块三：完整原 A_FF 五量作用，保留全部跨组连接。"""
    n = int(graph['cellCount'])
    ptr = np.atleast_1d(graph['rowPointer']).astype(np.int64) - 1
    col = np.atleast_1d(graph['column']).astype(np.int64) - 1
    order = np.atleast_1d(graph['order']).astype(np.int64) - 1
    owner = np.repeat(np.arange(n), np.diff(ptr))
    blocks = np.asarray(graph['values']).reshape(5, 5, len(col), order='F').transpose(2, 0, 1)
    ordered, output = value[order], np.zeros((n, 5))
    for first in range(0, len(col), 4096):
        ii = slice(first, first + 4096)
        np.add.at(output, owner[ii], np.einsum('kab,kb->ka', blocks[ii], ordered[col[ii]]))
    answer = np.empty_like(output)
    answer[order] = output
    return answer


def one_step(rhs, image):
    # 只计算一维最小残差系数，不把它作为新的物理步长或接受更新。
    ss, bb = norm(image), norm(rhs)
    assert ss > 0 and bb > 0
    alpha = float(np.sum((rhs / bb) * (image / ss))) * bb / ss
    return dict(unitStepRelativeResidual=norm(rhs - image) / bb,
        optimalSingleStepCoefficient=alpha,
        optimalSingleStepRelativeResidual=norm(rhs - alpha * image) / bb)


small = read_mat(VAL / 'overlapFunit01/restricted_overlap_unit_evidence.mat', ['cases'])
small_results = []
for case in small['cases']:
    f, b = case['factor'], case['rhs']
    core, weighted, counts = maps(f, b)
    dense, reference_counts = dense_weighted(f)
    np.testing.assert_array_equal(counts, reference_counts)
    dense_error = rel(weighted.ravel(), dense @ b.ravel())
    core_error = rel(core, case['single'])
    assert max(dense_error, core_error) < 1e-12
    other = np.cos(np.arange(b.size)).reshape(b.shape)
    combined = maps(f, .7*b - .3*other)[1]
    linear_error = rel(combined, .7*weighted - .3*maps(f, other)[1])
    assert linear_error < 1e-12
    plain = maps(case['plainFactor'], b)
    assert rel(plain[0], plain[1]) < 1e-12 and np.all(plain[2] == 1)
    assert norm(maps(f, np.zeros_like(b))[1]) == 0
    small_results.append(dict(denseMapRelativeDifference=dense_error,
        oldCoreMapRelativeDifference=core_error, linearityRelativeDifference=linear_error,
        zeroOverlapEqualsCoreMap=True, zeroRHSExactlyZero=True,
        # 存档 A 使用 MATLAB rhs(:) 的按列排列；上方独立新映射稠密参考使用逐CV排列。
        completeFImageDifference=rel(action(case['graph'], weighted).ravel(order='F'), case['A'] @ weighted.ravel(order='F'))))
    assert small_results[-1]['completeFImageDifference'] < 1e-12
del small, case, f, dense, combined, core, weighted, plain
gc.collect()
print('SMALL_DENSE_MAPS_PASSED', flush=True)

audit = json.loads((VAL / 'overlapFreal01_independent_restricted_overlap.json').read_text())
assert audit['status'].startswith('passed-')
factor_path = VAL / 'overlapFreal01/restricted_overlap_factor_evidence.mat'
graph_path = VAL / 'fullFreal01/full_F_graph_evidence.mat'
for file in (factor_path, graph_path):
    assert sha(file) == audit['inputSHA256'][str(file.relative_to(ROOT))]
factor = read_mat(factor_path, ['factor'])['factor']
graph = read_mat(graph_path, ['graph'])['graph']
base = read_mat(VAL / 'fullFreal01/full_F_real_evidence.mat', ['rhsF', 'Fimage'])
ras = read_mat(VAL / 'overlapFreal01/restricted_overlap_real_evidence.mat', ['rhsF', 'Fimage', 'single'])
np.testing.assert_array_equal(base['rhsF'], ras['rhsF'])
results = []
for label, b in [('original-F-RHS', base['rhsF']),
        ('best-old-32-step-F-residual', base['rhsF'] - base['Fimage']),
        ('restricted-overlap-32-step-F-residual', ras['rhsF'] - ras['Fimage'])]:
    core, weighted, counts = maps(factor, b)
    if label == 'original-F-RHS':
        assert rel(core, ras['single']) < 1e-12
    core_image, weighted_image = action(graph, core), action(graph, weighted)
    item = dict(rhs=label, coreOutput=one_step(b, core_image),
        membershipNormalizedOutput=one_step(b, weighted_image),
        membershipRange=[int(counts.min()), int(counts.max())],
        mappedVectorRelativeDifference=rel(weighted, core),
        normalizedOutputBetterOptimalSingleStep=False)
    item['normalizedOutputBetterOptimalSingleStep'] = item['membershipNormalizedOutput']['optimalSingleStepRelativeResidual'] < item['coreOutput']['optimalSingleStepRelativeResidual']
    results.append(item)
    print(json.dumps(item), flush=True)

assert all(sha(ROOT / file) == digest for file, digest in bindings.items())
out = dict(status='tested-read-only-core-versus-membership-normalized-overlap-action-screening',
    processID=os.getpid(), recordedUTC=datetime.now(timezone.utc).isoformat(),
    smallDenseReferences=small_results, realActionComparisons=results,
    method='D^-1 sum_g R_g^T A_g^-1 R_g; D is original-F membership count',
    certifiedExtendedInversesReusedByExactFileIdentity=True, newFactorizations=0,
    newKrylovSolves=0, fullFActionEvaluations=6, originalCOrSchurIdentityClaimed=False,
    physicalFieldAveraging=False, physicalStateChanged=False, physicalSolution=False,
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0,
    newProductionSolvedAreaM2=0, originalNSOrRegionalQualificationProven=False,
    screeningDoesNotProveKrylovConvergenceOrGlobalSpectrum=True,
    processPeakMemoryClaimed=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
TARGET.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], output=str(TARGET), sha256=sha(TARGET))))
