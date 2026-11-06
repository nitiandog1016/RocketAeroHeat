"""模块用途：定比热低焓参考与既有IdealAir核的独立物性/能量/局部归约筛选；不启动CFD。"""
import argparse
import csv
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import nacaConstantCpModel as model
from idealAirProperties import IdealAir
from naca3344Model import DATA, source_data, heat_parameter, recovery_temperature
from naca3344Normalization import published_recovery
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

PLAN = ROOT/'docs/validation/v174/PLAN.md'
SOURCE = ROOT/'tools/runtime/su2-thermal-build/source'
PARENT_TREE = ROOT/'outputs/v171_validation/build01/source_tree.json'
LIMITS = dict(cp=.02, conductivity=.02, deltaH=.02, sound=.01, localReduction=.03)
OFFICIAL_FILES = ('SU2_CFD/include/fluid/CIdealGas.hpp', 'SU2_CFD/src/fluid/CIdealGas.cpp',
    'SU2_CFD/include/fluid/CSutherland.hpp', 'SU2_CFD/include/fluid/CConstantPrandtl.hpp')


def screen(errors):
    """固定发展额度，非有限结果或缺项拒绝；超额度正常返回不入选，不能改阈值。"""
    if set(errors) != set(LIMITS) or not all(np.isfinite(x) and x >= 0 for x in errors.values()):
        raise ValueError('发展筛选指标缺失或非法')
    decisions = {name: dict(maximumRelativeDifference=float(errors[name]), limit=limit,
                            passed=bool(errors[name] <= limit)) for name, limit in LIMITS.items()}
    return dict(eligible=all(row['passed'] for row in decisions.values()), criteria=decisions,
                physicalQualified=False, experimentalQualified=False)


def reference(core, pressure, temperature):
    """旧模型来自已绑定DLL；cv/gamma由真实cp与同R得到，输运沿原明确Sutherland/Pr。"""
    old = core.state('PT', pressure, temperature)
    t, cp = old['T'], old['cp']
    cv = cp-model.R
    mu = 1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
    old.update(cv=cv, gamma=cp/cv, mu=mu, conductivity=mu*cp/.72)
    assert np.max(abs(old['sound']/np.sqrt(cp/cv*model.R*t)-1)) < 2e-13
    return old


def write_csv(path, columns):
    names = list(columns)
    size = len(columns[names[0]])
    assert all(len(columns[name]) == size for name in names)
    with path.open('x', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(names)
        writer.writerows(zip(*(columns[name] for name in names)))


def temperature_comparison(core, temperatures, tref, pressure, output, name):
    """共同Tref焓增量；零增量及1e-9K数值邻域明确采用cp导数极限。"""
    t = np.unique(np.asarray(temperatures, float))
    old = reference(core, pressure, t); new = model.state('PT', pressure, t)
    old_ref = float(core.state('PT', pressure, tref)['h'][0])
    old_cp_ref = float(core.state('PT', pressure, tref)['cp'][0])
    dh_old = old['h']-old_ref; dh_new = model.CP*(t-tref)
    zero = np.abs(t-tref) <= 1e-9
    dh_ratio = np.empty_like(t)
    dh_ratio[zero] = model.CP/old_cp_ref-1
    dh_ratio[~zero] = dh_new[~zero]/dh_old[~zero]-1
    columns = dict(T_K=t, pressure_Pa=np.full(t.shape, pressure),
        deltaH_reference_Jkg=dh_old, deltaH_constant_Jkg=dh_new,
        deltaH_relativeDifference=dh_ratio, deltaH_derivativeLimitUsed=zero.astype(int))
    summary = {}
    for key in ('rho', 'cp', 'cv', 'gamma', 'sound', 'mu', 'conductivity'):
        difference = new[key]/old[key]-1
        columns[key+'_reference'] = old[key]
        columns[key+'_constant'] = new[key]
        columns[key+'_relativeDifference'] = difference
        index = int(np.argmax(abs(difference)))
        summary[key] = dict(maximumAbsoluteRelativeDifference=float(abs(difference[index])),
            signedRelativeDifferenceAtMaximum=float(difference[index]), temperatureAtMaximumK=float(t[index]),
            referenceRange=[float(old[key].min()), float(old[key].max())],
            constantRange=[float(new[key].min()), float(new[key].max())])
    index = int(np.argmax(abs(dh_ratio)))
    summary['deltaH'] = dict(maximumAbsoluteRelativeDifference=float(abs(dh_ratio[index])),
        signedRelativeDifferenceAtMaximum=float(dh_ratio[index]), temperatureAtMaximumK=float(t[index]),
        commonReferenceTemperatureK=tref, derivativeLimitRows=int(zero.sum()),
        zeroIncrementPolicy='abs(T-Tref)<=1e-9 K: cp_constant/cp_reference(Tref)-1; no 0/0',
        rawAbsoluteEnthalpyCompared=False)
    write_csv(output/(name+'_properties.csv'), columns)
    return dict(sampleCount=len(t), temperatureRangeK=[float(t.min()), float(t.max())],
                samePressurePa=pressure, metrics=summary, continuousIntervalErrorBoundProven=False)


def local_reduction(core, contract, wall, heat, edge_temperatures, output):
    """相同且独立指定的edge p/T/U与合成q：只比较模型因子，绝非真实外缘或实验值。"""
    diameter = contract['geometry']['diameter_m']
    s = heat['s_over_D']*diameter
    assert s.min() >= wall['s_over_D'].min()*diameter and s.max() <= wall['s_over_D'].max()*diameter
    # 两个压力、两个速度独立指定；不取新旧来流代填实际边界层外缘。
    ss, te, pe, ue = np.meshgrid(s, edge_temperatures, [20000., 60000.], [100., 300.], indexing='ij')
    ss, te, pe, ue = (array.ravel() for array in (ss, te, pe, ue))
    tw = np.interp(ss/diameter, wall['s_over_D'], wall['Tw_K'])
    cr = published_recovery(ss)
    tr = recovery_temperature(te, contract['freeStream']['totalTemperature_K'], cr,
                              'published-Cr-curve-applied-to-common-synthetic-edge-for-model-comparison')
    old = reference(core, pe, te); new = model.state('PT', pe, te)
    results = []
    for values in (old, new):
        edge = dict(temperatureK=te, densityKgM3=values['rho'], speedMS=ue,
                    viscosityPas=values['mu'], conductivityWmK=values['conductivity'])
        results.append(heat_parameter(ss, np.full(te.shape, -1000.), tw, edge, tr,
            dict(edge='independently-prescribed-synthetic-common-p-T-U-not-a-flow-solution',
                 recovery='same-published-Cr-and-experiment-T0-for-both-models')))
    old_reduction, new_reduction = results
    ratio = new_reduction['parameter']/old_reduction['parameter']
    nusselt_ratio = new_reduction['Nusselt']/old_reduction['Nusselt']
    expected = old['conductivity']/new['conductivity']
    maximum_identity_error = float(np.max(abs(ratio/expected-1)))
    assert maximum_identity_error < 2e-13
    assert np.max(abs(new_reduction['Reynolds']/old_reduction['Reynolds']-1)) < 2e-13
    assert np.max(abs(new_reduction['heatTransferCoefficient']/old_reduction['heatTransferCoefficient']-1)) < 2e-13
    write_csv(output/'synthetic_common_edge_reduction.csv', dict(s_m=ss, edgeT_K=te, edgeP_Pa=pe,
        edgeU_m_s=ue, Tw_K=tw, published_Cr=cr, commonTr_K=tr,
        syntheticQIntoWall_Wm2=np.full(te.shape, -1000.),
        Reynolds_reference=old_reduction['Reynolds'], Reynolds_constant=new_reduction['Reynolds'],
        Nu_reference=old_reduction['Nusselt'], Nu_constant=new_reduction['Nusselt'],
        parameter_reference=old_reduction['parameter'], parameter_constant=new_reduction['parameter'],
        parameter_relativeDifference=ratio-1))
    return dict(cases=len(te), maximumRelativeDifference=float(np.max(abs(ratio-1))),
        maximumNusseltRelativeDifference=float(np.max(abs(nusselt_ratio-1))),
        analyticSameStateFactorIdentityMaximumError=maximum_identity_error,
        exactFactor='new/old = k_old/k_new = cp_old/cp_constant at common p,T,U,s,q,Tw,Tr',
        edgeTemperatureRangeK=[float(te.min()), float(te.max())], edgePressurePa=[20000., 60000.],
        edgeSpeedMS=[100., 300.], syntheticQIntoWallWm2=-1000.,
        localEdgeIsActual=False, qIsMeasured=False, measuredParameterUsedAsHeatFlux=False,
        sourceTransportTableEquivalent=False, measurementOperatorQualified=False,
        wholeMainPropertyRangeIsValidEdgeRange=False, physicalQualified=False)


def rejection_checks():
    """非法输入/模型身份与筛选失败分开；拒绝用例不改变已登记模型或原数据。"""
    cases = [
        ('unsupported-OTH-pair', lambda: model.state('Ph', 101325., 1.)),
        ('zero-pressure', lambda: model.state('PT', 0., 300.)),
        ('negative-temperature', lambda: model.state('PT', 101325., -1.)),
        ('zero-density', lambda: model.state('rhoe', 0., 1.)),
        ('negative-static-energy', lambda: model.state('rhoe', 1., -1.)),
        ('nan-temperature', lambda: model.state('PT', 101325., np.nan)),
        ('infinite-pressure', lambda: model.state('PT', np.inf, 300.)),
        ('two-dimensional-array', lambda: model.state('PT', [[1., 2.]], 300.)),
        ('incompatible-arrays', lambda: model.state('PT', [1., 2.], [200., 250., 300.])),
        ('empty-array', lambda: model.state('PT', [], [])),
        ('state-arithmetic-overflow', lambda: model.state('PT', 1e308, 1e308)),
        ('zero-Mach', lambda: model.freestream(0., 2.28e6, 570*5/9, .1016)),
        ('negative-Reynolds', lambda: model.freestream(1.97, -1., 570*5/9, .1016)),
        ('zero-diameter', lambda: model.freestream(1.97, 2.28e6, 570*5/9, 0.)),
        ('missing-screening-metric', lambda: screen({'cp': 0.})),
        ('nonfinite-screening-metric', lambda: screen(dict.fromkeys(LIMITS, np.nan))),
    ]
    answers = []
    for name, operation in cases:
        try:
            operation()
        except ValueError as error:
            answers.append(dict(case=name, rejected=True, reason=str(error)))
        else:
            raise AssertionError('应拒绝而未拒绝：'+name)
    over = dict.fromkeys(LIMITS, 0.); over['cp'] = .020001
    assert screen(over)['eligible'] is False
    assert model.in_development_range([150., 360.]).all()
    assert not model.in_development_range([149., 361.]).any()
    assert model.state('PT', 101325., 149.)['T'][0] == 149.
    return dict(passed=True, rejectedCases=answers, overAllocationEligibilityRejected=True,
                developmentRangeEligibilityRejectsOutside=True,
                positiveOutOfDevelopmentStateIsNotClipped=True, thresholdRelaxation=False)


def main(output):
    """模块二：一次≤40秒小报告；准确源、旧输入及所有拒绝/未入选终态均留存。"""
    started = time.monotonic()
    output.relative_to(ROOT/'outputs'); output.mkdir(parents=True, exist_ok=False)
    (output/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    report = dict(passed=False, reviewCompleted=False, CFDLaunches=0, builds=0, largeFieldsRead=0,
                  physicalQualified=False, experimentalQualified=False, maximumSeconds=40)
    locks = {}
    try:
        plan_text = PLAN.read_text(encoding='utf8')
        assert '物性/能量/归约≤40秒' in plan_text and 'SU2原生IDEAL_GAS' in plan_text
        core = IdealAir()
        contract, wall, heat = source_data()
        calculated = model.conditions()
        source_tree = read(PARENT_TREE)
        official = {str(SOURCE/name): source_tree[str(SOURCE/name)] for name in OFFICIAL_FILES}
        assert all(sha(Path(name)) == digest for name, digest in official.items())
        paths = [Path(__file__).resolve(), Path(model.__file__).resolve(), PLAN, PARENT_TREE,
            ROOT/'scripts/naca3344Model.py', ROOT/'scripts/naca3344Normalization.py',
            ROOT/'scripts/su2ThermalRuntime.py', ROOT/'scripts/validateIdealAir.py',
            DATA/'contract.json', DATA/'source_manifest.json', DATA/'wall_temperature.csv',
            DATA/'heat_parameter.csv', DATA/'normalization/manifest.json',
            DATA/'normalization/recovery_factor.csv']+[Path(name) for name in official]
        paths += [DATA/'normalization'/name for name in read(DATA/'normalization/manifest.json')['files']]
        locks = {str(path): sha(path) for path in paths}; locks.update(core.locks)
        save(output/'input_manifest.json', locks)
        for name in locks:
            path = Path(name)
            if path.suffix.lower() in ('.py', '.cpp', '.hpp', '.json', '.csv', '.md'):
                target = output/'exact_inputs'/path.relative_to(ROOT)
                target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(path.read_bytes())
        c = contract['freeStream']; t0 = c['totalTemperature_K']; tref = c['staticTemperature_K']
        h0 = float(core.state('PT', 101325., t0)['h'][0])
        def energy_residual(temperature):
            state = core.state('PT', 101325., temperature)
            return float(state['h'][0]+.5*(c['Mach']*state['sound'][0])**2-h0)
        old_t = brentq(energy_residual, 60., t0, xtol=2e-12, rtol=1e-14)
        old = reference(core, 101325., old_t)
        old_u = float(c['Mach']*old['sound'][0])
        old_rho = float(c['ReynoldsBasedOnDiameter']*old['mu'][0]/(old_u*contract['geometry']['diameter_m']))
        rebuilt = dict(staticTemperature_K=old_t, density_kg_m3=old_rho,
                       pressure_Pa=old_rho*model.R*old_t, speed_m_s=old_u,
                       viscosity_Pas=float(old['mu'][0]), totalEnthalpy_OTH_Jkg=h0)
        old_input_differences = {key: float(rebuilt[key]/c[key]-1) for key in rebuilt}
        assert max(abs(value) for value in old_input_differences.values()) < 2e-10
        common_keys = ('staticTemperature_K', 'density_kg_m3', 'pressure_Pa', 'speed_m_s', 'viscosity_Pas')
        new_input_differences = {key: float(calculated[key]/rebuilt[key]-1) for key in common_keys}
        new_state = model.state('PT', calculated['pressure_Pa'], calculated['staticTemperature_K'])
        conservative_roundtrip = model.state('rhoe', new_state['rho'], new_state['e'])
        roundtrip_error = max(float(np.max(abs(conservative_roundtrip[key]/new_state[key]-1))) for key in new_state)
        anchors_error = dict(Mach=abs(calculated['speed_m_s']/calculated['soundSpeed_m_s']/c['Mach']-1),
            Reynolds=abs(calculated['density_kg_m3']*calculated['speed_m_s']*calculated['diameter_m']/calculated['viscosity_Pas']/c['ReynoldsBasedOnDiameter']-1),
            totalEnthalpy=abs((calculated['staticEnthalpy_Jkg']+.5*calculated['speed_m_s']**2)/calculated['totalEnthalpy_Jkg']-1))
        assert max(anchors_error.values()) < 2e-13 and roundtrip_error < 2e-13
        lower = min(old_t, calculated['staticTemperature_K'], tref); upper = float(wall['Tw_K'].max())
        main_t = np.unique(np.r_[np.linspace(lower, upper, 1025), wall['Tw_K'], tref, old_t,
                                 calculated['staticTemperature_K'], t0])
        extended_t = np.unique(np.r_[np.linspace(150., 360., 841), main_t, tref])
        primary = temperature_comparison(core, main_t, tref, 30000., output, 'primary')
        diagnostic = temperature_comparison(core, extended_t, tref, 30000., output, 'development')
        edge_t = np.unique(np.r_[np.linspace(lower, t0, 65), calculated['staticTemperature_K']])
        reduction = local_reduction(core, contract, wall, heat, edge_t, output)
        # 扩展区间的纯模型归约因子使用同p/T/U及同正Tw-Tr的形式因子；不伪造Te>T0的恢复温度。
        factor_old = reference(core, 30000., extended_t)['conductivity']
        factor_new = model.state('PT', 30000., extended_t)['conductivity']
        extended_factor = float(np.max(abs(factor_old/factor_new-1)))
        def metrics(comparison, normalization):
            return dict(cp=comparison['metrics']['cp']['maximumAbsoluteRelativeDifference'],
                conductivity=comparison['metrics']['conductivity']['maximumAbsoluteRelativeDifference'],
                deltaH=comparison['metrics']['deltaH']['maximumAbsoluteRelativeDifference'],
                sound=comparison['metrics']['sound']['maximumAbsoluteRelativeDifference'],
                localReduction=normalization)
        primary_screen = screen(metrics(primary, reduction['maximumRelativeDifference']))
        development_screen = screen(metrics(diagnostic, extended_factor))
        integral_checks = []
        for endpoint in (150., 250., t0, upper, 360.):
            integrated = quad(lambda temperature: float(core.state('PT', 30000., temperature)['cp'][0]),
                              tref, endpoint, epsabs=2e-7, epsrel=2e-12)[0]
            exact = float(core.state('PT', 30000., endpoint)['h'][0]-core.state('PT', 30000., tref)['h'][0])
            error = abs(integrated-exact)/max(abs(exact), 1.)
            assert error < 2e-10
            integral_checks.append(dict(endpointK=endpoint, integratedCpDeltaH_Jkg=integrated,
                                        kernelDeltaH_Jkg=exact, scaledDifference=error))
        rejections = rejection_checks()
        definitions = dict(
            cp='定压比热 J/(kg K)；原生IDEAL_GAS为gamma*R/(gamma-1)，参考值由既有IdealAir DLL独立返回。',
            cv='定容比热cp-R；gamma=cp/cv，两者分别报告但不私增发展误差阈值。',
            sound='a=sqrt(gamma*R*T)，单位m/s；参考DLL直接声速同时核对热力恒等式。',
            viscosity='mu=mu_ref*(T/T_ref)^1.5*(T_ref+S)/(T+S)，Pa s；两模型同T时采用完全相同明确输运。',
            conductivity='k=mu*cp/Pr，W/(m K)；同T时k差来自cp，未宣称等于原论文参考26输运表。',
            energy='新模型e=cv*T、h=cp*T，无OTH偏移。旧模型保留OTH；绝对h/e不可作为模型误差直接相减。',
            deltaH='同一Tref=旧输入T∞比较h(T)-h(Tref)，零增量用cp导数极限，表内逐行标记。',
            localReduction='Nu_s=h_transfer*s/k_e；Re_s=rho_e*U_e*s/mu_e；h_transfer=-q_into/(Tw-Tr)。同edge p/T/U/s/q/Tw/Tr时新旧Nu和Nu/sqrt(Re)比等于k_old/k_new。',
            syntheticEdge='p=20000/60000Pa、U=100/300m/s均独立指定，仅为模型比较；不来自来流替代、不是真实外缘，不保证每组满足全局流动。',
            developmentFactor='扩展[150,360]K只算共态输运归约因子k_old/k_new；Te>T0部分不代入原恢复温度式、不称可用外缘。',
            uncertainty='原Nu约15%、Re约1%为probable error，分别列出，不当标准差合成；模型小额度不授实验资格。')
        report.update(passed=True, reviewCompleted=True, model=calculated,
            newModelEvaluation='Python SI transcription of bound official CIdealGas/Sutherland/Pr formulas; no C++ build or execution',
            referenceModelEvaluation='existing independently bound IdealAir batch DLL; original OTH reference retained',
            originalAnchors=dict(Mach=c['Mach'], Re_D=c['ReynoldsBasedOnDiameter'], T0_K=t0,
                                 D_m=contract['geometry']['diameter_m'], wallPoints=22, heatParameters=20),
            rebuiltOriginalIdealAir=rebuilt, originalInputRelativeReconstructionDifferences=old_input_differences,
            newToOriginalDerivedInputRelativeDifferences=new_input_differences,
            newAnchorClosure=anchors_error, PT_rhoe_roundtripMaximumRelativeError=roundtrip_error,
            comparisonTemperatureRangeK=[lower, upper], developmentTemperatureRangeK=[150., 360.],
            primaryComparison=primary, developmentComparison=diagnostic,
            localReduction=reduction, developmentPureTransportReductionFactorMaximumDifference=extended_factor,
            commonReferenceEnthalpyIntegralChecks=integral_checks, rejectionChecks=rejections,
            primaryScreen=primary_screen, developmentScreen=development_screen,
            developmentEntryEligible=bool(primary_screen['eligible'] and development_screen['eligible']),
            officialSourceSHA256=official, definitions=definitions,
            originalExperimentalProbableErrors=dict(Nusselt=.15, Reynolds=.01, combinedAsStandardDeviation=False),
            originalTransportTableEquivalent=False, actualEdgeQualified=False,
            latterSixCenterIdentityClosed=False, originalContractModified=False,
            absoluteEnthalpyDifferenceUsedForEligibility=False)
        lines = ['v174 定比热低焓模型比较（仅发展入选筛选，非实验通过）',
            '新模型：gamma=1.4，R=8.31451/0.02896546，原Sutherland参数和Pr=0.72。',
            '原M/Re_D/T0/D及22点Tw/20个热参数保持，重建静态条件，未复用OTH守恒能量。',
            '比较温区 %.12g—%.12g K；开发检查150—360 K不是已解全场包络。'%(lower, upper),
            '发展入口入选：%s；物理资格：False。'%report['developmentEntryEligible'],
            '新T∞=%.12g K，rho∞=%.12g kg/m3，p∞=%.12g Pa，U∞=%.12g m/s。'%
            tuple(calculated[key] for key in ('staticTemperature_K','density_kg_m3','pressure_Pa','speed_m_s'))]
        lines += [key+'：'+value for key, value in definitions.items()]
        (output/'explanation.txt').write_text('\n\n'.join(lines)+'\n', encoding='utf8')
    except Exception as error:
        report['error'] = type(error).__name__+': '+str(error)
        raise
    finally:
        stable = {name: sha(Path(name)) == digest for name, digest in locks.items()}
        elapsed = time.monotonic()-started
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=elapsed,
                      budgetPassed=elapsed <= 40)
        report['passed'] = bool(report['passed'] and all(stable.values()) and elapsed <= 40)
        save(output/'report.json', report)
        assert all(stable.values()) and elapsed <= 40
    print(dict(passed=report['passed'], developmentEntryEligible=report['developmentEntryEligible'],
               elapsedSeconds=report['elapsedSeconds']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    main(Path(args.output).resolve())
