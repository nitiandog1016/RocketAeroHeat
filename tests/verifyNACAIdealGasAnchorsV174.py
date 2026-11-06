"""模块用途：以50位Decimal独立重建NACA定比热来流和C1初猜；不启动CFD。"""
import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]


def check(out):
    """模块一：仅读取实验锚点，旧Lemmon派生来流保留为比较项。"""
    start = time.monotonic()
    out = Path(out).resolve()
    out.relative_to(ROOT / 'outputs')
    assert not out.exists(), '不能覆盖既有检查'
    source = ROOT / 'tests/data/naca_tn3344/contract.json'
    raw = source.read_bytes()
    original = json.loads(raw.decode('utf-8-sig'))
    old = original['freeStream']
    with localcontext() as ctx:
        ctx.prec = 50
        gamma = Decimal('1.4')
        gas = Decimal('8.31451') / Decimal('.02896546')
        diameter = Decimal(str(original['geometry']['diameter_m']))
        mach = Decimal(str(old['Mach']))
        reynolds = Decimal(str(old['ReynoldsBasedOnDiameter']))
        total = Decimal(str(old['totalTemperature_R'])) * Decimal(5) / Decimal(9)
        temp = total / (1 + (gamma - 1) * mach * mach / 2)
        sound = (gamma * gas * temp).sqrt()
        speed = mach * sound
        viscosity = (Decimal('1.716e-5') * (temp / Decimal('273.15')) ** 1
                     * (temp / Decimal('273.15')).sqrt()
                     * Decimal('383.55') / (temp + Decimal('110.4')))
        density = reynolds * viscosity / (speed * diameter)
        pressure = density * gas * temp
        cp = gamma * gas / (gamma - 1)
        cv = gas / (gamma - 1)
        exact = dict(gamma=gamma, R=gas, cp=cp, cv=cv, diameterM=diameter,
                     Mach=mach, ReD=reynolds, totalTemperatureK=total,
                     staticTemperatureK=temp, soundSpeedMS=sound, speedMS=speed,
                     viscosityPas=viscosity, densityKgM3=density, pressurePa=pressure,
                     totalEnthalpyCpT0Jkg=cp * total)
    values = {key: float(value) for key, value in exact.items()}

    # 模块二：用独立浮点守恒量反算T、p与全焓，检查12个解析初态。
    maximum_error = 0.
    rows = []
    for wall in original['wall']['TwRange_K']:
        assert wall > values['totalTemperatureK'] > values['staticTemperatureK'] > 0
        for blend in (0., .001, .1, .5, .9, 1.):
            target_h = (1 - blend) * values['cp'] * wall + blend * values['totalEnthalpyCpT0Jkg']
            velocity = -blend * values['speedMS']
            temperature = ((1 - blend) * wall + blend * values['totalTemperatureK']
                - blend * blend * (values['totalTemperatureK'] - values['staticTemperatureK']))
            rho = values['pressurePa'] / (values['R'] * temperature)
            momentum = rho * velocity
            rho_energy = rho * (values['cv'] * temperature + .5 * velocity * velocity)
            restored_t = (rho_energy / rho - .5 * (momentum / rho) ** 2) / values['cv']
            restored_p = (values['gamma'] - 1) * (rho_energy - .5 * momentum * momentum / rho)
            restored_h = (rho_energy + restored_p) / rho
            errors = [abs(restored_t / temperature - 1),
                      abs(restored_p / values['pressurePa'] - 1), abs(restored_h / target_h - 1)]
            maximum_error = max(maximum_error, *errors)
            assert all(math.isfinite(v) and v > 0 for v in (rho, temperature, rho_energy))
            assert values['staticTemperatureK'] - 1e-12 <= temperature <= wall + 1e-12
            rows.append(dict(TwK=wall, blend=blend, temperatureK=temperature,
                             velocityXMS=velocity, density=rho, rhoEnergy=rho_energy))
    assert maximum_error < 2e-13
    anchor_errors = {
        'Mach': abs(values['speedMS'] / values['soundSpeedMS'] / values['Mach'] - 1),
        'ReD': abs(values['densityKgM3'] * values['speedMS'] * values['diameterM']
                   / values['viscosityPas'] / values['ReD'] - 1),
        'T0': abs((values['cp'] * values['staticTemperatureK'] + .5 * values['speedMS'] ** 2)
                  / values['totalEnthalpyCpT0Jkg'] - 1),
        'EOS': abs(values['densityKgM3'] * values['R'] * values['staticTemperatureK']
                   / values['pressurePa'] - 1),
    }
    assert max(anchor_errors.values()) < 2e-13

    # 模块三：两种典型错误有明确辨别力；不是新生产入口拒绝测试。
    wrong_mach_vector_speed = values['soundSpeedMS']
    wrong_oth_temperature = (old['totalEnthalpy_OTH_Jkg'] - .5 * values['speedMS'] ** 2) / values['cp']
    assert abs(wrong_mach_vector_speed / values['speedMS'] - 1) > .4
    assert abs(wrong_oth_temperature / values['staticTemperatureK'] - 1) > .1
    locks = {str(source): hashlib.sha256(raw).hexdigest(),
             str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    stable = {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest for path, digest in locks.items()}
    assert all(stable.values())
    result = dict(passed=True, scope='independent-50-digit-anchored-algebra-no-production-validator',
        values=values, decimalTruth={key: str(value) for key, value in exact.items()},
        anchorMaximumRelativeErrors=anchor_errors, initialStateMaximumRelativeRoundtripError=maximum_error,
        initialStates=rows, monotonicTemperatureProof='dT/db=T0-Tw-2*b*(T0-Tinf)<0 for Tw>T0>Tinf',
        oldDerivedStateUsedForInitialization=False, priorSuccessfulUpdates=0,
        wrongUnitMachVectorSpeedMS=wrong_mach_vector_speed,
        wrongOTHTotalEnthalpyGivesStaticTemperatureK=wrong_oth_temperature,
        positiveAlgebraCases=12, discriminatingWrongInputExamples=2, productionContractRefusalsTested=0,
        modelComparabilityTested=False, physicalQualified=False, experimentalQualification=False,
        newQualifiedArea_m2=0, newCFDTasks=0, builds=0,
        inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - start)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    (out / 'report.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf8')
    print(json.dumps(dict(passed=True, values=values, maximumError=maximum_error,
                         elapsedSeconds=result['elapsedSeconds'])))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    check(args.output)
