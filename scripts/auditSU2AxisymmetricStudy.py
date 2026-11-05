"""审核公开轴对称任务，失败历史、源场诊断和新检查点严格分别记录。"""
import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json
from su2AxisymmetricAnalysis import (mesh_data, fields, boundary_balance, residual_scales,
                                     wall_values, sample_stations, upstream_profile, integrate_scalar)
from su2SSTInterface import wall_geometry
from su2BenchmarkAcceptance import acceptance, stability


def assessment(path, mesh_path, conditions, contract, spatial):
    """核对任意明确身份的保存场，返回量值和误差，资格由调用方判定。"""
    p, cells, markers = mesh_data(mesh_path)
    f = fields(path, spatial)
    rho = f['Density']
    t = f['Temperature']
    eos = float(np.max(abs(f['Pressure']/(rho*conditions['gasConstant']*t)-1)))
    assert eos < 2e-12 and min(t) > 0 and min(f['Pressure']) > 0
    g = wall_geometry(mesh_path)
    ids = g['wallIDs']
    cp = conditions['gamma']*conditions['gasConstant']/(conditions['gamma']-1)
    tw = conditions['wallTemperatureK']
    q = f['Laminar_Viscosity'][ids]*cp/.72*(t[g['normalNeighbor']]-tw)/g['thermalDistance']
    heat_error = float(np.max(abs(q-f['Heat_Flux'][ids])/np.maximum(abs(q), 1)))
    assert heat_error < 2e-6
    wall = wall_values(p, markers, f)
    stations = sample_stations(wall, contract)
    profile = upstream_profile(p, f, conditions, contract)
    balance = boundary_balance(p, cells, markers, f, conditions)
    flux = integrate_scalar(p, markers['wall'], f['Heat_Flux'])
    reference = dict(pressurePa=float(np.interp(-.06, wall['S'], wall['p'])),
                     heatFluxWm2=float(np.interp(-.06, wall['S'], wall['q'])))
    heat = stations['heat']
    peak_pred = int(np.argmax(heat['predicted']))
    peak_measured = int(np.argmax(heat['measured']))
    crossing = []
    for i in np.flatnonzero(wall['cf'][:-1]*wall['cf'][1:] < 0):
        s = wall['S'][i]-wall['cf'][i]*(wall['S'][i+1]-wall['S'][i])/(wall['cf'][i+1]-wall['cf'][i])
        crossing.append(dict(S_m=float(s), kind='separation' if wall['cf'][i] > 0 else 'reattachment'))
    return dict(EOSRelativeError=eos, heatGradientRelativeError=heat_error,
                savedWallFluidTemperatureRangeK=[float(min(t[ids])), float(max(t[ids]))],
                maximumWallFluidTemperatureDifferenceK=float(max(abs(t[ids]-tw))),
                wall=wall, stations=stations, profile=profile, boundaryBalance=balance,
                residuals=residual_scales(p, cells, f, conditions) if spatial else None,
                revolvedHeatPowerW=flux, localReference=reference,
                yPlus99=float(np.percentile(wall['yplus'], 99)), shearZeroCrossings=crossing,
                measuredStationHeatPeak=dict(predictedWm2=float(heat['predicted'][peak_pred]),
                    measuredWm2=float(heat['measured'][peak_measured]),
                    relativeMagnitudeError=float(abs(heat['predicted'][peak_pred]/heat['measured'][peak_measured]-1)),
                    positionErrorM=float(abs(heat['S'][peak_pred]-heat['S'][peak_measured]))))


def serial(value):
    """仅把数值数组转换成明确JSON数组，不用缺失值补零。"""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serial(v) for v in value]
    return value


def plot_history(history, labels, report, target):
    """科研图显示实际原生历史，末次历史步与保存场身份同时标出。"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    iteration = history[:, labels.index('Inner_Iter')]
    for ax, names, title in (
            (axes[0, 0], ['rms[Rho]', 'rms[RhoE]'], 'Flow residuals'),
            (axes[0, 1], ['rms[k]', 'rms[w]'], 'SST-Vm residuals')):
        for name in names:
            ax.plot(iteration, history[:, labels.index(name)], label=name, linewidth=1.2)
        ax.axhline(-8, color='.5', linestyle=':', label='Original threshold')
        ax.set_ylabel('Native log10 RMS')
        ax.set_title(title)
        ax.legend(fontsize=8)
    axes[1, 0].plot(iteration, history[:, labels.index('HF')]/1000, color='#36688d')
    axes[1, 0].set_ylabel('Native planar heat integral (kW/m)')
    axes[1, 0].set_title('Native HF history')
    for name in ('Min DT', 'Max DT'):
        axes[1, 1].semilogy(iteration, history[:, labels.index(name)], label=name)
    axes[1, 1].set_ylabel('Local pseudo-time step (s)')
    axes[1, 1].set_title('Recorded time steps')
    axes[1, 1].legend(fontsize=8)
    for ax in axes.flat:
        ax.set_xlabel('Local recorded iteration')
        ax.grid(alpha=.2)
    status = 'terminal exit %d' % report['exitCode']
    frame = 'new field unavailable' if report['savedLocalIteration'] is None else 'saved field %d' % report['savedLocalIteration']
    fig.suptitle('NASA ACC-SBLI M7.11 | SST-Vm | '+status+' | '+frame, fontsize=12)
    fig.savefig(str(target), dpi=150)
    plt.close(fig)


def audit(folder, target):
    """先审终态与输出身份，再决定当前结果是否有可分析的新场。"""
    folder = Path(folder).resolve()
    target = Path(target).resolve()
    assert not target.exists()
    target.mkdir(parents=True)
    receipt = json.loads((folder/'engine_process.json').read_text(encoding='utf-8'))
    manifest = json.loads((folder/'study_manifest.json').read_text(encoding='utf-8'))
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    for path, expected in receipt['sourceSHA256'].items():
        assert sha(path) == expected, path
    for name, expected in receipt['outputSHA256'].items():
        assert sha(folder/name) == expected, name
    for path, expected in manifest['sourceSHA256'].items():
        assert sha(path) == expected, path
    c = manifest['benchmark']['conditions']
    contract_path = ROOT/'tests/data/tmr_mach7/contract.json'
    contract = json.loads(contract_path.read_text(encoding='utf-8'))
    with (folder/'history.csv').open() as f:
        rows = list(csv.reader(f, skipinitialspace=True))
    labels = [s.strip() for s in rows[0]]
    history = np.array(rows[1:], dtype=float)
    assert np.isfinite(history).all()
    assert np.array_equal(history[:, labels.index('Inner_Iter')], np.arange(len(history)))
    recorded = int(history[-1, labels.index('Inner_Iter')])
    files = sorted(folder.glob('restart_flow_*.csv'))
    current = None
    step = None
    current_error = None
    if files:
        checkpoint = files[-1]
        step = int(re.fullmatch(r'restart_flow_(\d+)\.csv', checkpoint.name)[1])
        assert checkpoint.name in receipt['outputSHA256'] and step <= recorded
        try:
            current = assessment(checkpoint, folder/'mesh.su2', c, contract, True)
        except (AssertionError, ValueError, KeyError) as error:
            current_error = type(error).__name__+': '+str(error)
    # 源场单独命名。原任务早期失败时，源场帮助复核处理链并保持其真实时刻。
    independent = manifest.get('independentInitialization')
    source = None if independent else assessment(folder/'initial_flow.csv', folder/'mesh.su2', c, contract, False)
    residual_names = ['rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]']
    residuals = {n: float(history[-1, labels.index(n)]) for n in residual_names}
    frames = []
    frame_errors = []
    p, _, markers = mesh_data(folder/'mesh.su2')
    for saved in files:
        iteration = int(re.fullmatch(r'restart_flow_(\d+)\.csv', saved.name)[1])
        if iteration < recorded-1100:
            continue
        try:
            f = fields(saved, True)
            stations = sample_stations(wall_values(p, markers, f), contract)
            frames.append(dict(iteration=iteration, p=stations['pressure']['predicted'],
                               q=stations['heat']['predicted'],
                               power=integrate_scalar(p, markers['wall'], f['Heat_Flux'])))
        except (AssertionError, ValueError, KeyError) as error:
            frame_errors.append(dict(iteration=iteration, error=type(error).__name__+': '+str(error)))
    steady = stability(frames, recorded, manifest['benchmark']['geometry']['revolvedWallAreaM2'])
    succeeded = receipt['exitCode'] == 0 and receipt['stopReason'] is None
    decision = acceptance(current, residuals, steady, succeeded, step == recorded and current is not None)
    report = dict(schema='RocketAeroHeat.axisymmetric-study-review.v1', auditCompleted=True,
                  reviewScope='preregistered numerical and station acceptance; domain qualification separate',
                  convergenceAcceptanceImplemented=True,
                  directory=str(folder), model=manifest['benchmark']['model'], conditions=c,
                  exitCode=receipt['exitCode'], stopReason=receipt['stopReason'],
                  elapsedSeconds=receipt['elapsedSeconds'], lastRecordedIteration=recorded,
                  sourceLocalIteration=manifest['sourceLocalIteration'], savedLocalIteration=step,
                  runSucceeded=succeeded,
                  hasSavedNewField=bool(files), validationPassed=decision['passed'],
                  physicsQualified=False, newQualifiedArea=0,
                  finalRecordedResiduals=residuals,
                  originalResidualThresholdPassed=all(v < -8 for v in residuals.values()),
                  currentFieldAssessment=current,
                  currentFieldError=current_error, checkpointFrameErrors=frame_errors,
                  acceptance=decision, temporalStability=steady,
                  sourceInitialFieldAssessment=source,
                  independentInitialization=independent,
                  sourceOriginalStartupPassed=manifest['sourceOriginalStartupPassed'],
                  missingCurrentChecks=decision['pending']+decision['furtherQualificationRequired'],
                  numericalFailure='NaN detected' if 'NaN detected' in (folder/'engine.log').read_text() else None,
                  inputsStable=True,
                  sourceSHA256={str(p): sha(p) for p in (Path(__file__).resolve(),
                      ROOT/'scripts/su2AxisymmetricAnalysis.py', ROOT/'scripts/su2SSTInterface.py',
                      ROOT/'scripts/su2BenchmarkAcceptance.py',
                      folder/'engine_process.json', folder/'study_manifest.json', contract_path)})
    plot_history(history, labels, report, target/'convergence.png')
    report['historyPlot'] = dict(path=str(target/'convergence.png'), sha256=sha(target/'convergence.png'))
    save_json(target/'report.json', serial(report))
    print(json.dumps({k: report[k] for k in ('auditCompleted', 'exitCode', 'lastRecordedIteration',
                     'savedLocalIteration', 'hasSavedNewField', 'originalResidualThresholdPassed',
                     'numericalFailure', 'elapsedSeconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    audit(args.folder, args.output)
