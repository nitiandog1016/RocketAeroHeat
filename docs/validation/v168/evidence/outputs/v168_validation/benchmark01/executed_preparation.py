"""模块用途：整理NACA实测半球柱热交换基准，图读值、单位和物性匹配可追溯。"""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.optimize import brentq

from idealAirProperties import IdealAir
from su2ThermalRuntime import ROOT, read
from validateIdealAir import save

DATA = ROOT/'tests/data/naca_tn3344'
SOURCE = ROOT/'outputs/v168_validation/sources'
OUT = ROOT/'outputs/v168_validation/benchmark01'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digitize(spec):
    """模块一：从已人工选择的点和三点标定换算，第四网格点独立核验。"""
    anchors = np.array(spec['pixelCalibration'], float)
    matrix = np.linalg.solve(np.column_stack((anchors[:, :2], np.ones(3))), anchors[:, 2:])
    px = np.array(spec['pixels'], float)
    values = np.column_stack((px, np.ones(len(px))))@matrix
    check = np.array(spec['independentAxisCheck'], float)
    error = np.abs(np.r_[check[:2], 1.]@matrix-check[2:])
    assert error[0] < .004 and error[1] < .0015, error
    assert np.all(np.diff(values[:, 0]) > 0)
    return values, error


def write_table(path, names, rows):
    with path.open('w', newline='', encoding='utf8') as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        writer.writerows(rows)


def main():
    """模块二：先形成可审查条件与测量合同，不生成网格或运行求解器。"""
    start = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_preparation.py').write_bytes(Path(__file__).read_bytes())
    input_files = [Path(__file__).resolve(), DATA/'digitization.json', SOURCE/'NACA_TN_3344.pdf',
                   SOURCE/'NACA_TN_3344_source.json', SOURCE/'tn3344-41.png', SOURCE/'tn3344-43.png',
                   SOURCE/'tn3344-46.png', SOURCE/'tn3344-49.png']
    locks = {str(p): sha(p) for p in input_files}
    source = read(SOURCE/'NACA_TN_3344_source.json')
    assert sha(SOURCE/'NACA_TN_3344.pdf') == source['sha256']
    raw = DATA/'raw'
    raw.mkdir(exist_ok=False)
    for p, name in [('tn3344-41.png', 'table1.png'), ('tn3344-43.png', 'figure1c.png'),
                    ('tn3344-46.png', 'figure3a.png'), ('tn3344-49.png', 'figure4b.png')]:
        shutil.copyfile(str(SOURCE/p), str(raw/name))
    spec = read(DATA/'digitization.json')
    heat, heat_check = digitize(spec['heat'])
    wall, wall_check = digitize(spec['wallTemperature'])
    raw_nose_s = float(wall[0, 0])
    assert abs(raw_nose_s) <= .006
    wall[0, 0] = 0.
    assert min(wall[:, 1]) > 1. and min(heat[:, 1]) > 0
    D, M, ReD, T0 = .1016, 1.97, 2.28e6, 570.*5/9
    model = IdealAir()
    locks.update(model.locks)
    h0 = float(model.state('PT', 101325, T0)['h'][0])
    def total_h(T):
        value = model.state('PT', 101325, T)
        return float(value['h'][0]+.5*M*M*value['sound'][0]**2)
    T = brentq(lambda t: total_h(t)-h0, 100., T0, xtol=1e-11)
    thermo = model.state('PT', 101325, T)
    U = M*float(thermo['sound'][0])
    mu = 1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4)
    rho = ReD*mu/(U*D)
    R = 287.0491267875601
    p = rho*R*T
    matched = model.state('PT', p, T)
    mach_error = abs(U/float(matched['sound'][0])-M)
    re_error = abs(float(matched['rho'][0])*U*D/mu/ReD-1)
    h_error = abs((float(matched['h'][0])+.5*U*U-h0)/max(1., abs(h0)))
    assert max(mach_error, re_error, h_error) < 2e-12
    # 独立常比热解析值只作模型差异对照，不同时强迫两套状态量相等。
    Tperfect = T0/(1.+.2*M*M)
    write_table(DATA/'heat_parameter.csv', ['s_over_D', 'Nu_s_over_sqrt_Re_s',
        'reading_uncertainty_s_over_D', 'reading_uncertainty_parameter', 'source_pixel_x', 'source_pixel_y'],
        [[*v, .004, .005, *px] for v, px in zip(heat, spec['heat']['pixels'])])
    write_table(DATA/'wall_temperature.csv', ['s_over_D', 'Tw_over_T0', 'Tw_K',
        'reading_uncertainty_ratio', 'source_pixel_x', 'source_pixel_y'],
        [[v[0], v[1], v[1]*T0, .0015, *px] for v, px in zip(wall, spec['wallTemperature']['pixels'])])
    manifest = dict(source, document='NACA TN 3344, Stine and Wanlass, December 1954',
        sourcePage='https://ntrs.nasa.gov/citations/19930084172',
        downloadedPDFKeptLocal=True, pdfPrefixBytes=436,
        pdfRead='pypdf recovered xref; Poppler successfully rendered and figures visually inspected',
        files=[dict(path=x.relative_to(DATA).as_posix(), bytes=x.stat().st_size, sha256=sha(x))
               for x in sorted(raw.glob('*.png'))],
        digitizationSHA256=sha(DATA/'digitization.json'),
        measurements=[dict(path=x.name, sha256=sha(x)) for x in (DATA/'heat_parameter.csv', DATA/'wall_temperature.csv')])
    save(DATA/'source_manifest.json', manifest)
    contract = dict(schema='RocketAeroHeat.naca-tn3344-heat.v1',
        caseID='NACA-TN3344-M1.97-ReD2.28e6-heated-hemisphere-cylinder',
        status='selected-primary-experiment-with-digitized-conditions; solver-and-normalization-adapter-pending',
        geometry=dict(diameter_m=D, noseRadius_m=D/2, tipToCylinderEnd_m=.127,
            endRingThickness_m=.003175, incidence_deg=0., incidenceTolerance_deg=.1,
            physicalShape='hemisphere plus cylinder, stainless steel electrically heated shell; Figure 1(c)',
            coordinate='s is meridional arclength from tip; project X=0 at tip, nose +X, flow -X',
            hemisphereMapping='theta=2*s/D; X=-R*(1-cos(theta)), r=R*sin(theta)',
            cylinderMapping='s>=pi*R/2; X=-R-(s-pi*R/2), r=R',
            experimentalStingAndElectricalHardware='not axisymmetric; downstream computational truncation must have an independence check',
            measuredHeatSurfaceOnly=True, aftBaseAndStingHeatMeasured=False),
        freeStream=dict(Mach=M, ReynoldsBasedOnDiameter=ReD, totalTemperature_R=570., totalTemperature_K=T0,
            source='Table I heat-transfer distribution, wind tunnel 1, final row',
            independentInputs=['Mach', 'Re_D', 'T0', 'D', 'gas/viscosity model'],
            uniformUpstreamFlow=True, prescribedIncomingBoundaryLayer=False,
            gas='dry air, laminar; do not enable SST to match this experiment',
            model='Lemmon2000 ideal Air OTH, Sutherland 1.716e-5/273.15/110.4, Pr=0.72 proposed',
            staticTemperature_K=T, pressure_Pa=p, density_kg_m3=rho, speed_m_s=U,
            viscosity_Pas=mu, totalEnthalpy_OTH_Jkg=h0, perfectGasGamma14Temperature_K=Tperfect,
            derivedStateIsMeasured=False, propertyChoiceRequiresDeclaredComparison=True),
        wall=dict(material='stainless steel with electrical resistance heating; given thermal boundary, no solid solve',
            measuredProfile='wall_temperature.csv', TwRange_K=[float(min(wall[:, 1])*T0), float(max(wall[:, 1])*T0)],
            isExactlyUniform=False, source='Figure 3(a), upper triangles, same Re_D=2.28e6',
            temperatureDataBeyondLastStation=None,
            extensionRule='none yet: prepare explicit downstream continuation and independence check before solving',
            sign='body is hotter than recovery temperature: heat flows OUT of wall; project positive-into-wall q is negative'),
        measurements=dict(file='heat_parameter.csv', points=len(heat), source='Figure 4(b), circles with dot',
            observable='Nu_s/sqrt(Re_s)', definitions='Nu_s=h*s/k_e; Re_s=rho_e*U_e*s/mu_e; h=q_out/(Tw-Tr)',
            localProperties='edge of boundary layer; s is arclength, NOT axial X; NOT freestream Re_D',
            sourceReduction='published heat parameter uses measured electrical input, local outer state and recovery correction (equations 10-12)',
            rawDimensionalHeatFluxAvailable=False,
            forbidden='cannot label these values as directly measured W/m2, use Re_D in local normalization, or use theoretical curves as measurements',
            publishedNusseltProbableError=0.15, publishedReynoldsProbableError=0.01,
            sourceReductionUncertainty='includes small Tw nonuniformity/conduction; measured unheated surface is not true recovery T near nose s/D<0.35',
            acceptance='not yet qualified: finalize normalization and experiment/reading/mesh uncertainties before any candidate'),
        nextImplementation=dict(engine='existing embedded mature SU2, dedicated laminar axisymmetric entry',
            scope='source data/units and measured-Tw boundary adapter, single geometry, internal external-flow solution',
            prerequisites=['freeze heat-parameter comparison including edge/recovery definitions',
                           'bind measured Tw(s) and downstream extension without default 300 K',
                           'resolve source laminar no-turbulence 2D contract and axis boundary',
                           'preregister numerical convergence, heat stability, grid/domain uncertainty and one budgeted candidate'],
            noNewKernelPatchWithoutSmallEvidence=True),
        experimentalQualification=False, newQualifiedArea_m2=0., originalMach7GatesUnchanged=True,
        sourceManifestSHA256=sha(DATA/'source_manifest.json'))
    save(DATA/'contract.json', contract)
    stable = {path: sha(Path(path)) == value for path, value in locks.items()}
    assert all(stable.values())
    save(OUT/'report.json', dict(passed=True, heatPoints=len(heat), wallTemperaturePoints=len(wall),
        heatIndependentAxisError=heat_check.tolist(), wallIndependentAxisError=wall_check.tolist(),
        rawNoseS_over_D=raw_nose_s, MachAbsoluteError=mach_error, ReDRelativeError=re_error,
        totalEnthalpyRelativeError=h_error, newNativeTasks=0, newQualifiedArea_m2=0,
        inputsSHA256=locks, inputsStable=stable, contractSHA256=sha(DATA/'contract.json'),
        elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True, heatPoints=len(heat), wallPoints=len(wall),
        staticTemperatureK=T, pressurePa=p, wallTemperatureRangeK=contract['wall']['TwRange_K'])))


if __name__ == '__main__':
    main()
