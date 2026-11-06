"""模块用途：用独立积分与纸面物性公式复核本轮证据，并生成科研对照图。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import paper_cp, save


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def main(output_name='verify01'):
    start = time.monotonic()
    assert Path(output_name).name == output_name
    out = ROOT/'outputs/v168_validation'/output_name
    out.mkdir(exist_ok=False)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    review = read(ROOT/'outputs/v168_validation/review02/report.json')
    errors = []
    frames = []
    for row in review['public']:
        data = read(ROOT/'outputs/v159_validation/public01/public_review02'/('frame_%06d.json' % row['localIteration']))
        w = data['wall']
        # 用多项式展开积分 r(s)*(T(s)-Tw)^2，独立于原两点高斯实现。
        s, r, e = (np.array(w[a]) for a in ('S', 'radius', 'fluidT'))
        e -= 311.
        dr, de = np.diff(r), np.diff(e)
        terms = r[:-1]*e[:-1]**2+(2*r[:-1]*e[:-1]*de+dr*e[:-1]**2)/2
        terms += (r[:-1]*de**2+2*dr*e[:-1]*de)/3+dr*de**2/4
        area = 2*np.pi*np.sum(np.diff(s)*(r[:-1]+r[1:])/2)
        rms = np.sqrt(2*np.pi*np.sum(np.diff(s)*terms)/area)
        errors.append(abs(rms-row['areaRMSWallTError_K']))
        frames.append(data)
    assert max(errors) < 1e-11
    final = review['public'][-1]
    assert abs(final['upstreamFluidT_K']-311.) < .001
    assert final['upstreamQ_Wm2']/9300 < .005
    assert abs(final['maxTErrorX_m']-.06) > .2
    accepted = review['inheritedPublicAcceptance']
    assert accepted['passed'] is False
    assert sum(g['passed'] is True for g in accepted['gates']) == 4
    assert len(accepted['failed']) == 21 and len(accepted['pending']) == 1

    # 独立论文cp积分重建总焓差和声速，避免只逆用同一个批处理函数。
    data = ROOT/'tests/data/naca_tn3344'
    contract = read(data/'contract.json')
    c = contract['freeStream']
    T, T0, U = c['staticTemperature_K'], c['totalTemperature_K'], c['speed_m_s']
    R = 8.31451/.02896546
    dh = quad(paper_cp, T, T0, epsabs=1e-8, epsrel=2e-12)[0]
    energy_error = abs(dh-.5*U**2)/(.5*U**2)
    a = np.sqrt(paper_cp(T)/(paper_cp(T)-R)*R*T)
    mach_error = abs(U/a-c['Mach'])
    assert energy_error < 2e-10 and mach_error < 2e-11
    manifest = read(data/'source_manifest.json')
    for entry in manifest['files']+manifest['measurements']:
        assert hashlib.sha256((data/entry['path']).read_bytes()).hexdigest() == entry['sha256']
    heat = np.genfromtxt(str(data/'heat_parameter.csv'), delimiter=',', names=True)
    wall = np.genfromtxt(str(data/'wall_temperature.csv'), delimiter=',', names=True)
    assert len(heat) == 20 and len(wall) == 22
    assert heat['s_over_D'][0] > 0 and wall['s_over_D'][0] == 0
    assert np.all(wall['Tw_K'] > T0)  # 热从壁流向气体；不得把实验正系数当正入壁热流。
    assert np.max(abs(wall['Tw_K']/T0-wall['Tw_over_T0'])) < 1e-14

    # 图只展示既存未合格场与源数据，不把复审画图当作主UI新验收。
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.5), constrained_layout=True)
    w = frames[-1]['wall']; station = frames[-1]['stations']['heat']
    ax[0, 0].semilogy(w['S'], np.maximum(w['q'], .01), color='#286086', label='v159 saved update 1200')
    ax[0, 0].scatter(station['S'], station['measured'], color='#ac5836', s=20, label='Kussoy experiment')
    ax[0, 0].scatter([-.06], [9300], marker='x', color='black', label='Upstream experiment')
    ax[0, 0].set(xlabel='Surface arclength S [m]', ylabel='Heat into wall [W/m²]',
                 title='Public case: heat mismatch remains')
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot(w['S'], np.array(w['fluidT'])-311, color='#286086')
    ax[0, 1].axvline(-.06, ls='--', color='#ac5836', label='Upstream measurement station')
    ax[0, 1].set(xlabel='Surface arclength S [m]', ylabel='Fluid T - given Tw [K]',
                 title='Largest wall-T error is downstream')
    ax[0, 1].legend(fontsize=8)
    for f, label, style in [(frames[0], 'Update 1', '--'), (frames[-1], 'Update 1200', '-')]:
        p = f['profile']
        ax[1, 0].plot(p['predictedU'], 1000*np.array(p['y']), style, label=label)
    p = frames[-1]['profile']
    ax[1, 0].scatter(p['measuredU'], 1000*np.array(p['y']), s=16, color='#ac5836', label='Experiment')
    ax[1, 0].set(xlabel='Axial velocity / measured reference velocity', ylabel='Distance from wall [mm]',
                 ylim=(0, 12), title='Incoming profile is still poorly matched')
    ax[1, 0].legend(fontsize=8)
    rr = review['rocket']
    ax[1, 1].plot([x['seriesUpdates'] for x in rr], [x['maximumWallTError_K'] for x in rr], 'o-', label='Maximum')
    ax[1, 1].plot([x['seriesUpdates'] for x in rr], [x['areaRMSWallTError_K'] for x in rr], 'o-', label='Area RMS')
    ax[1, 1].set(xlabel='Recorded cumulative updates (cold restart)', ylabel='Wall-fluid T error [K]',
                 title='Rocket: improvement is not qualification')
    ax[1, 1].legend(fontsize=8)
    for axes in ax.ravel(): axes.grid(alpha=.18)
    fig.suptitle('v168 evidence review — no new CFD; physical validation remains incomplete', fontsize=13)
    fig.savefig(out/'wall_heat_review.png', dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    ax[0].errorbar(heat['s_over_D'], heat['Nu_s_over_sqrt_Re_s'], yerr=.005, fmt='o', ms=4, capsize=2,
                   color='#286086', label='Digitized experimental circles')
    ax[0].set(xlabel='Meridional arclength s / D', ylabel='Nu_s / sqrt(Re_s)',
              title='Measured heat-transfer parameter (not W/m²)')
    ax[0].legend(fontsize=8)
    ax[1].errorbar(wall['s_over_D'], wall['Tw_K'], yerr=.0015*T0, fmt='^-', ms=4, lw=1, capsize=2, color='#ac5836')
    ax[1].set(xlabel='Meridional arclength s / D', ylabel='Measured wall temperature [K]',
              title='Heated-wall triangles; T0 = 316.667 K')
    for axes in ax:
        axes.axvline(np.pi/4, ls='--', lw=1, color='#777777')
        axes.grid(alpha=.18)
    fig.suptitle('Selected primary benchmark: NACA TN 3344, M=1.97, Re_D=2.28e6\nBars show reading uncertainty only; published Nu uncertainty is ±15%', fontsize=12)
    fig.savefig(out/'naca3344_digitized.png', dpi=170)
    plt.close(fig)
    save(out/'report.json', dict(passed=True, publicIntegralMaxAbsoluteDifference_K=max(errors),
        independentPaperEnergyRelativeError=energy_error, independentPaperMachAbsoluteError=mach_error,
        sourceImageAndMeasurementHashesPassed=True, nativeSolverRuns=0, MATLABRuns=0,
        newPhysicsQualification=False, elapsedSeconds=time.monotonic()-start))
    print(json.dumps(read(out/'report.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-name', default='verify01')
    main(parser.parse_args().output_name)
