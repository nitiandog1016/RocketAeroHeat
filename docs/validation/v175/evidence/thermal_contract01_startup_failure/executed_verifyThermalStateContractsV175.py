"""用途：独立有理数反例区分状态正性/物理温域，并区分两种成熟远场重构；不执行CFD。"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GAMMA = F(7, 5)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')


def square_root(value):
    """本检查有意选有理声速；不使用生产物性或边界程序生成期望结果。"""
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    assert n * n == value.numerator and d * d == value.denominator
    return F(n, d)


def conservative(rho, velocity, pressure):
    assert rho > 0 and pressure > 0
    return [rho, rho * velocity, pressure / (GAMMA - 1) + rho * velocity**2 / 2]


def primitive(values, gas_constant=F(1)):
    rho, momentum, energy = values
    assert rho > 0
    velocity = momentum / rho
    rho_e = energy - momentum**2 / (2 * rho)
    assert rho_e > 0
    pressure = (GAMMA - 1) * rho_e
    return dict(rho=rho, velocity=velocity, pressure=pressure, temperature=pressure / (rho * gas_constant))


def as_numbers(values):
    return {key: float(value) for key, value in values.items()}


def boundary_states(inside, outside):
    """法向一维嵌入：局部守恒特征投影与非线性R±/熵重构；不是原生通量执行。"""
    pi, pe = primitive(inside), primitive(outside)
    rho, u, p = pi['rho'], pi['velocity'], pi['pressure']
    re, ue, pext = pe['rho'], pe['velocity'], pe['pressure']
    a, ae = square_root(GAMMA * p / rho), square_root(GAMMA * pext / re)
    enthalpy = (inside[2] + p) / rho
    columns = [[F(1), u-a, enthalpy-u*a], [F(1), u, u*u/2], [F(1), u+a, enthalpy+u*a]]
    dr, dm, dE = [right-left for left, right in zip(inside, outside)]
    dp = (GAMMA-1)*(dE-u*dm+u*u*dr/2)
    momentum_jump = dm-u*dr
    alphas = [(dp-a*momentum_jump)/(2*a*a), dr-dp/(a*a), (dp+a*momentum_jump)/(2*a*a)]
    eigenvalues = [u-a, u, u+a]
    projected = [inside[i] + sum(columns[j][i]*alphas[j] for j in range(3) if eigenvalues[j] < 0) for i in range(3)]
    # 精确守恒跳跃分解，以及独立矩阵求逆核对解析投影系数。
    assert [sum(columns[j][i]*alphas[j] for j in range(3)) for i in range(3)] == [dr, dm, dE]
    matrix = np.asarray([[float(columns[j][i]) for j in range(3)] for i in range(3)])
    independent = np.linalg.solve(matrix, np.asarray([float(value) for value in (dr, dm, dE)]))
    error = float(np.max(abs(independent-np.asarray([float(value) for value in alphas]))))
    assert error < 2e-15
    rplus = u + 2*a/(GAMMA-1) if ue > -ae else ue + 2*ae/(GAMMA-1)
    rminus = u - 2*a/(GAMMA-1) if ue > ae else ue - 2*ae/(GAMMA-1)
    ub, ab = (rplus+rminus)/2, (rplus-rminus)*(GAMMA-1)/4
    assert ab > 0
    source_rho, source_a = (rho, a) if ue > 0 else (re, ae)
    assert 2/(GAMMA-1) == 5
    rb = source_rho*(ab/source_a)**5
    pb = rb*ab*ab/GAMMA
    farfield = conservative(rb, ub, pb)
    return projected, farfield, error


def verify(out, exact_source, expected_sha):
    """冻结计划、准确源与历史小证据；计算仅为独立代数，不读取原生大场。"""
    start = time.monotonic()
    out, exact_source = Path(out).resolve(), Path(exact_source).resolve()
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks = {}
    report = dict(passed=False, newCFDTasks=0, builds=0, largeFieldsRead=False,
        physicalQualified=False, experimentalQualification=False, v174FailureReclassified=False,
        oldTemperatureGuardChanged=False, originalCauseEstablished=False)

    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        assert expected is None or digest == expected, str(path)
        locks[str(path)] = digest

    try:
        source = Path(__file__).resolve()
        bind(source)
        (out / 'executed_verifyThermalStateContractsV175.py').write_bytes(source.read_bytes())
        bind(ROOT / 'docs/validation/v175/PLAN.md')
        bind(exact_source, expected_sha)
        # 只绑定旧配置/回执/小审核，不读restart或体场，也不重新执行旧验证。
        for name in ('engines/warm/case.cfg', 'engines/warm/engine_process.json',
                'checks/plate_warm01/warm_start_manifest.json', 'checks/plate_warm01/physical_audit.json',
                'source_versions/prepareSU2FlatPlateWarmStart.py.4865173ac3120306.raw',
                'source_versions/auditSU2WarmConvergence.py.8f30d49d1137354b.raw'):
            bind(ROOT / 'docs/validation/v108' / name)
        text = exact_source.read_text(encoding='utf8')
        assert 'void CEulerSolver::BC_Far_Field' in text and 'void CEulerSolver::BC_Riemann' in text
        assert 'if (Lambda_i[jVar] < 0)' in text and 'Qn_Infty > -SoundSpeed_Infty' in text
        assert 'conv_numerics->SetPrimitive(V_domain, V_infty)' in text
        R = F(831451, 100000) / F(2896546, 100000000)
        cv = R/(GAMMA-1)
        # 合法正内能可超出物理验证温域：仅是数学反例，不重新接纳v174真实任务。
        high = conservative(F(1), F(0), R*450)
        high_state = primitive(high, R)
        assert high_state['temperature'] == 450 and high_state['pressure'] > 0
        assert not 150 <= high_state['temperature'] <= 360
        rejected = []
        for name, values in [('zero-density', [F(0), F(0), F(1)]),
                ('negative-density', [F(-1), F(0), F(1)]), ('zero-internal-energy', [F(1), F(2), F(2)]),
                ('negative-internal-energy', [F(1), F(2), F(1)])]:
            try:
                primitive(values, R)
            except AssertionError:
                rejected.append(name)
            else:
                raise AssertionError('非正状态未拒绝')
        # 守恒凸映射保正，但速度方差变成内能，不能推导端点温度最大值原理。
        left, right = conservative(F(1), F(0), R*300), conservative(F(1), F(1000), R*300)
        mixed = [(a+b)/2 for a, b in zip(left, right)]
        mixed_state = primitive(mixed, R)
        expected_temperature = F(300) + F(125000)/cv
        assert mixed_state['temperature'] == expected_temperature and expected_temperature > 360
        # 给定Tw与壁流体T身份不同；保持k不变只是展示热壁离散式里的独立变量。
        conductivity, neighbour, tw, dn = F(1, 50), F(320), F(340), F(1, 10000)
        q_into_wall = conductivity*(neighbour-tw)/dn
        assert q_into_wall == -4000
        assert conductivity*(neighbour-F(450))/dn == -26000
        # 有限子声速正状态反例：两种边界重构不同；并非旧任务失效原因证明。
        inside = conservative(F(1), F(1, 5), F(5, 7))
        outside = conservative(F(1), F(1, 10), F(81, 140))
        riemann, farfield, error = boundary_states(inside, outside)
        assert riemann != farfield
        r_state, f_state = primitive(riemann), primitive(farfield)
        assert f_state['velocity'] == F(2, 5) and f_state['rho'] == F(24, 25)**5
        same_r, same_f, same_error = boundary_states(inside, inside)
        assert same_r == same_f == inside
        incoming_i = conservative(F(1), F(-2), F(5, 7))
        incoming_e = conservative(F(1), F(-3), F(81, 140))
        incoming_r, incoming_f, incoming_error = boundary_states(incoming_i, incoming_e)
        assert incoming_r == incoming_f == incoming_e
        report.update(passed=True, algebraicPositiveOutOfDevelopmentDomain=as_numbers(high_state),
            independentPositiveStateRefusals=rejected,
            conservativeConvexInterpolation=dict(endpointTemperaturesK=[300., 300.],
                endpointVelocitiesMS=[0., 1000.], mixture=as_numbers(mixed_state),
                exactExpectedTemperatureK=float(expected_temperature), positivityPreserved=True,
                endpointTemperatureMaximumPrinciple=False, actualNACAInterpolationClaim=False),
            weakWallIdentity=dict(givenTwK=float(tw), neighbourTemperatureK=float(neighbour),
                conductivityFixedForAlgebraOnly=float(conductivity), normalDistanceM=float(dn),
                qIntoWallWm2=float(q_into_wall), replacingGivenTwWithFluid450ProducesWrongFluxWm2=-26000.,
                fullConductivityMayDependOnFluidState=True, productionBoundaryVerified=False),
            characteristicCounterexample=dict(dimensionless=True, outwardNormal=1,
                interior=as_numbers(primitive(inside)), external=as_numbers(primitive(outside)),
                localConservativeProjection=as_numbers(r_state), nonlinearFarfieldReconstruction=as_numbers(f_state),
                differentPositiveReconstructedStates=True, independentMatrixError=error,
                sameStateIdentityPassed=True, supersonicIncomingIdentityPassed=True,
                maximumIndependentMatrixError=max(error, same_error, incoming_error),
                nativeSolverExecuted=False, finalUpwindFluxCompared=False, provesNACACause=False),
            evidenceScope='independent-exact-rational-algebra-and-bounded-read-only-source-binding')
    except Exception as error:
        report.update(passed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {name: Path(name).is_file() and sha(Path(name)) == digest for name, digest in locks.items()}
        if not all(stable.values()):
            report['passed'] = False
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic()-start)
        save(out / 'report.json', report)
    assert report['passed'] and all(stable.values())
    print(dict(passed=True, refusals=len(rejected), boundaryReconstructionDifferent=True, elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--exact-euler-source', required=True)
    parser.add_argument('--exact-euler-sha256', required=True)
    args = parser.parse_args()
    verify(args.output, args.exact_euler_source, args.exact_euler_sha256)
