"""读取完整或紧凑QR坐标基，以原物理排列独立计算投影。"""
import numpy as np
from fractions import Fraction
from scipy.sparse import lil_matrix


def conditioned_rows_reference(H, condition):
    """用二进制浮点数的精确有理表示验证可逆行变换，不以QR秩选行。"""
    q = H.shape[0]
    changed = np.atleast_1d(condition['transformedRowIndices']).astype(int)-1
    references = np.atleast_1d(condition['referenceRowIndices']).astype(int)-1
    powers = np.atleast_1d(condition['powerTwoExponents']).astype(int)
    assert len(changed) == len(references) == len(powers) and len(set(changed)) == len(changed)
    assert condition['inputDistinctRowCount'] == condition['outputDistinctConditionCount'] == q
    assert not condition['originalHChanged'] and not condition['rankBasedRowSelection'] and not condition['rankTruncationUsed']
    assert condition['allOriginalConditionsRetained'] and condition['exactDifferenceBySterbenzVerified']
    W = H.tolil(copy=True); checked = 0
    for j, k, power in zip(changed, references, powers):
        assert 0 <= k < j < q and k not in changed
        a = H.getrow(j); b = H.getrow(k)
        np.testing.assert_array_equal(a.indices, b.indices)
        assert np.all(np.signbit(a.data) == np.signbit(b.data))
        assert np.all(abs(a.data) >= .5*abs(b.data)) and np.all(abs(a.data) <= 2*abs(b.data))
        exact = [(Fraction(float(x))-Fraction(float(y)))*Fraction(2)**int(power) for x, y in zip(a.data, b.data)]
        values = np.array([float(x) for x in exact])
        for value, truth, x, y in zip(values, exact, a.data, b.data):
            assert Fraction(float(value)) == truth
            assert Fraction(float(y))+Fraction(float(value))*Fraction(2)**(-int(power)) == Fraction(float(x))
            checked += 1
        assert np.any(values)
        W[j, a.indices] = values
    return W.tocsr(), dict(transformedRows=len(changed), originalConditions=q,
                            exactRationalCoefficientIdentities=checked, rankTruncationUsed=False)


def project_coordinates(projector, coordinates):
    value = np.asarray(coordinates).ravel(order='F').copy()
    if 'basisCoordinateIndices' in projector:
        ids = np.atleast_1d(projector['basisCoordinateIndices']).astype(int)-1
        assert len(np.unique(ids)) == len(ids) and np.all((ids >= 0) & (ids < len(value)))
    else:
        ids = np.arange(len(value))
    U = np.asarray(projector['basis']).reshape(len(ids), -1, order='F')
    np.testing.assert_allclose(U.T@U, np.eye(U.shape[1]), atol=1e-12, rtol=0)
    value[ids] -= U@(U.T@value[ids])
    return value.reshape(np.shape(coordinates), order='F')


def certify_projector_span(projector, H, transpose_right_map):
    """逐个独立原行计算完整P转置，不把紧凑支撑正确性作为输入假设。"""
    reps = np.atleast_1d(projector['exactRepresentativeRowIndices']).astype(int)-1
    mapping = np.atleast_1d(projector['originalToUniqueRow']).astype(int)-1
    assert (H-H[reps[mapping]]).nnz == 0
    ids = (np.atleast_1d(projector['basisCoordinateIndices']).astype(int)-1
           if 'basisCoordinateIndices' in projector else np.arange(H.shape[1]))
    U = np.asarray(projector['basis']).reshape(len(ids), len(reps), order='F')
    np.testing.assert_allclose(U.T@U, np.eye(len(reps)), atol=1e-12, rtol=0)
    outside = np.ones(H.shape[1], dtype=bool); outside[ids] = False
    errors = []; norms = []
    for row in reps:
        column = transpose_right_map(H.getrow(row).toarray().ravel())
        assert np.all(column[outside] == 0), '完整P转置H转置在声称支撑外不为零'
        norm = np.linalg.norm(column); assert np.isfinite(norm) and norm > 0
        error = float(np.linalg.norm(column[ids]-U@(U.T@column[ids]))/norm)
        assert error < 1e-10
        errors.append(error); norms.append(norm)
    np.testing.assert_allclose(np.asarray(norms)[mapping], np.atleast_1d(projector['constraintColumnNorms']), rtol=1e-12)
    conditioning = None
    if 'constraintRowPreconditioning' in projector:
        working, conditioning = conditioned_rows_reference(H[reps], projector['constraintRowPreconditioning'])
        images = []
        for j in range(working.shape[0]):
            column = transpose_right_map(working.getrow(j).toarray().ravel())
            assert np.all(column[outside] == 0)
            norm = np.linalg.norm(column)
            image = column[ids]/norm
            assert np.linalg.norm(image-U@(U.T@image)) < 1e-10
            images.append(image)
            np.testing.assert_allclose(norm, np.atleast_1d(projector['conditionedConstraintColumnNorms'])[j], rtol=1e-12)
        _, R = np.linalg.qr(np.column_stack(images), mode='reduced')
        pivots = abs(np.diag(R)); assert np.all(pivots > 1e-10)
        conditioning['minimumIndependentNormalizedQRDiagonal'] = float(pivots.min())
    if 'basisCoordinateIndices' in projector:
        assert U.nbytes == projector['basisArrayBytes']
        assert projector['additionalConstructionArrayBoundBytes'] <= projector['originalConstructionBudgetBytes'] == 448*1024**2
    result = dict(originalRows=H.shape[0], distinctConditions=len(reps), storedCoordinateCount=len(ids),
                compactStorage='basisCoordinateIndices' in projector,
                maximumCompleteRightMapSpanDifference=max(errors), allOmittedCoordinatesExactlyZero=True)
    if conditioning is not None: result['exactRowPreconditioning'] = conditioning
    return result
