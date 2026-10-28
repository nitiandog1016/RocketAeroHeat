"""显示真实连续步的残差、方向门槛和温域激活；不推定未保存状态或物理收敛。"""
import json
import os, sys
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha

OUT = ROOT / 'outputs/v87_validation'
if len(sys.argv) > 1:
    supplied = Path(sys.argv[1])
    OUT = (supplied if supplied.is_absolute() else ROOT / supplied).resolve()
    OUT.relative_to(ROOT.resolve())
study_file = OUT / 'continuous01/active_continuation.json'


def main():
    study = json.loads(study_file.read_text('utf8'))
    steps = study['steps'] if isinstance(study['steps'], list) else [study['steps']]
    actual = [x for x in steps if x['main'] and x['main']['acceptedNonlinearIterations']]
    if actual:
        initial = actual[0]['main']; initial_rms = initial['initialRmsMerit']; initial_peak = initial['initialMaximumResidual']
    else:
        initial = loadmat(study['initialStateFile'], simplify_cells=True, variable_names=['report'])['report']
        initial_rms = initial['rmsMerit']; initial_peak = initial['maximumResidual']
    rms = [initial_rms] + [x['main']['rmsMerit'] for x in actual]
    peak = [initial_peak] + [x['main']['maximumResidual'] for x in actual]
    x = np.arange(len(rms))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    a = axes[0, 0]
    a.plot(x, np.array(rms) / rms[0], 'o-', color='#236e9c', label='RMS / initial RMS')
    a.plot(x, np.array(peak) / peak[0], 's-', color='#857046', label='Peak / initial peak')
    a.set_ylabel('Residual relative to current-study start')
    a.set_xlabel('Actual accepted continuation step'); a.set_xticks(x); a.legend()
    a = axes[0, 1]
    number = np.arange(1, len(steps) + 1)
    a.plot(number, [s['directionReport']['linearReport']['relativeCResidual'] for s in steps], 'o-', label='Approximate C')
    a.plot(number, [s['probe']['originalNSDirectionDefect'] if s['probe'] else np.nan for s in steps], 's-', label='Fresh original NS')
    a.axhline(.1, color='#a74b39', ls='--', label='Original NS target 0.1')
    a.axhline(.001, color='#777', ls=':', label='Internal C target 0.001')
    a.set_yscale('log'); a.set_ylabel('Relative direction defect'); a.legend(fontsize=8)
    a.set_xticks(number)
    a = axes[1, 0]
    alphas = []
    bindings = {str(study_file.relative_to(ROOT)): sha(study_file)}
    for step in actual:
        file = __import__('pathlib').Path(step['resultFolder']) / 'wall_velocity_krylov_evidence.mat'
        bindings[str(file.relative_to(ROOT))] = sha(file)
        result = loadmat(file, simplify_cells=True, variable_names=['result'])['result']
        trace = result['trace']
        if not isinstance(trace, dict):
            trace = trace[-1]
        alphas.append(float(trace['alpha']))
    a.bar(np.arange(1, len(alphas) + 1), alphas, color='#236e9c')
    if alphas:
        a.set_yscale('log')
    else:
        a.text(.5, .5, 'No accepted flow update', transform=a.transAxes, ha='center')
    a.set_ylabel('Actually accepted step alpha')
    a.set_xticks(np.arange(1, len(alphas)+1))
    for j, alpha in enumerate(alphas, 1):
        a.annotate(f'{alpha:g}', (j, alpha), xytext=(0, 5), textcoords='offset points', ha='center')
    a = axes[1, 1]
    counts = [len(np.atleast_1d(s['directionReport']['bandFaceIDs'])) if np.size(s['directionReport']['bandFaceIDs']) else 0 for s in steps]
    active = [len(np.atleast_1d(s['directionReport']['activeFaceIDs'])) if np.size(s['directionReport']['activeFaceIDs']) else 0 for s in steps]
    interface = any(s['directionReport'].get('lookaheadSelection', {}).get('selectionPolicy') ==
                    'original-and-predicted-cell-and-interface-faces' for s in steps)
    targeted = interface or any(s['directionReport'].get('lookaheadSelection', {}).get('selectionPolicy') ==
                               'original-and-predicted-required-faces' for s in steps)
    if targeted:
        original = [np.size(s['directionReport']['lookaheadSelection']['originalBandFaceIDs']) for s in steps]
        a.bar(number - .24, original, .22, label='Original interfaces within 10 K')
        a.bar(number, counts, .22, label=('Original + predicted CV/interface required faces' if interface
                                        else 'Original + all predicted required faces'))
        a.bar(number + .24, active, .22, label='Outward rows actually activated')
        a.set_ylim(0, max(original+counts+active)*1.4)
    else:
        a.bar(number - .17, counts, .32, label='Current interfaces within 10 K of bounds')
        a.bar(number + .17, active, .32, label='Outward rows actually activated')
    a.set_ylabel('Interface count'); a.legend(fontsize=8)
    a.set_xticks(number); a.set_xlabel('Current-state direction study step')
    for ax in axes.ravel():
        ax.grid(True, axis='y', alpha=.2)
    fig.suptitle('M5 / 30 km / 20 deg / prescribed Tw=300 K / UNCONVERGED RESEARCH\nCurrent derivatives rebuilt every step; original NS and 80-2000 K gates unchanged; qualified area added = 0', fontsize=11)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'temperature_active_continuation.{ext}', dpi=160)
    plt.close(fig)
    assert all(sha(ROOT / name) == digest for name, digest in bindings.items())
    receipt = dict(status='plotted-actual-current-state-continuation-and-temperature-direction-policy',
                   actualAcceptedSteps=len(actual), acceptedAlpha=alphas,
                   rmsReductionPercent=100 * (1 - rms[-1] / rms[0]),
                   peakReductionPercent=100 * (1 - peak[-1] / peak[0]),
                   inputSHA256=bindings, scriptSHA256=sha(__file__),
                   physicalSolution=False, newProductionSolvedAreaM2=0,
                   browserVisualOrWebGLRuntimeVerified=False)
    target = OUT / 'scientific_plot.json'; assert not target.exists()
    target.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf8')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
