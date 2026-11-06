"""模块用途：仅从小历史/限步事件/准确回执绘制v173科研图，不读取场或定位未知最终失败。"""
import argparse
import csv
from pathlib import Path
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from plotNACAAdaptiveHistory import history
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save


def plot(folder, output, termination_review):
    """模块一：锁定已终止任务的小文件，局部/系列编号按原checkpoint语义核对。"""
    started = time.monotonic()
    folder.relative_to(ROOT/'outputs'); output.relative_to(ROOT/'outputs')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed_plot.py').write_bytes(Path(__file__).read_bytes())
    inputs = [Path(__file__).resolve(), ROOT/'scripts/plotNACAAdaptiveHistory.py',
        ROOT/'scripts/su2ThermalRuntime.py', ROOT/'scripts/validateIdealAir.py',
        folder/'history.csv', folder/'conservative_steps.csv', folder/'engine_process.json',
        folder/'case_manifest.json', folder/'step_initialization.json', termination_review,
        folder.parent/'native_real01_process.json']
    assert all(path.is_file() and path.stat().st_size <= 8*1024**2 for path in inputs)
    locks = {str(path): sha(path) for path in inputs}
    save(output/'input_manifest.json', locks)
    for path in inputs:
        target = output/'exact_inputs'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    receipt = read(folder/'engine_process.json'); manifest = read(folder/'case_manifest.json')
    init = read(folder/'step_initialization.json'); termination = read(termination_review)
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 134
    assert all(receipt['inputsStable'].values())
    assert not any(receipt.get(key, False) for key in ('timedOut', 'cancelled', 'outputBudgetExceeded',
        'memoryBudgetExceeded', 'supervisorError'))
    assert init == manifest['initialization']
    assert manifest['requiredRuntimeFlavor'] == 'rah-ideal-air-spatial-step-v1'
    assert termination['passed'] and termination['terminationActualCauseKnown'] is False
    assert termination['recoveredEventsIdentifyFinalFailure'] is False
    assert termination['inputSHA256'][str(folder/'engine_process.json')] == locks[str(folder/'engine_process.json')]
    assert termination['inputSHA256'][str(folder/'conservative_steps.csv')] == locks[str(folder/'conservative_steps.csv')]
    for name in ('history.csv', 'conservative_steps.csv'):
        assert locks[str(folder/name)] == receipt['outputSHA256'][name]
    data = history(folder/'history.csv'); x = data['Inner_Iter'].astype(int)
    assert len(x) == 57 and int(x[-1]) == receipt['checkpoints']['lastFiniteHistoryLocalIteration'] == 56
    offset = int(init['priorSuccessfulUpdates'])+1
    series = offset+x
    assert int(series[-1]) == receipt['checkpoints']['seriesUpdatesAtLastFiniteHistory'] == 1558
    saved = receipt['checkpoints']['complete']
    assert [(row['localIteration'], row['seriesSuccessfulUpdates']) for row in saved] == [(0, 1502)]
    with (folder/'conservative_steps.csv').open(encoding='utf8', newline='') as stream:
        events = list(csv.DictReader(stream))
    assert len(events) == 4 and [int(row['event']) for row in events] == [1, 2, 3, 4]
    assert [int(row['inner_iter']) for row in events] == [46, 56, 56, 56]
    event_rows = []
    for row in events:
        local = int(row['inner_iter'])
        assert row['context_known'] == row['old_valid'] == row['output'] == '1'
        assert float(row['theta']) == .5 and float(row['k']) == 0
        assert 'mode=rhoe phase=node.trial' in row['reason']
        event_rows.append(dict(event=int(row['event']), localIteration=local,
            seriesUpdateNumber=offset+local, localRCMPoint=int(row['local_point']),
            globalPoint=int(row['global_point']), theta=float(row['theta']),
            nativeAlpha=float(row['alpha']), outputStage=True, oldStateValid=True,
            queryMode='rhoe', queryPhase='node.trial', temperatureDirectionClassified=False))

    # 模块二：四面板共享局部编号；系列轴只表示来源更新计数，不表示不中断等价。
    plt.rcParams.update({'font.family': 'Microsoft YaHei', 'axes.unicode_minus': False,
        'font.size': 9.5, 'svg.fonttype': 'path', 'axes.spines.top': False,
        'axes.spines.right': False})
    fig, axes = plt.subplots(4, 1, figsize=(11.8, 11.1), sharex=True,
        gridspec_kw={'height_ratios': [2.1, 2.0, 2.0, 1.7]})
    fig.subplots_adjust(left=.105, right=.975, bottom=.095, top=.81, hspace=.37)
    fig.suptitle('NACA TN3344 · v173 同守恒方向限步 + 空间 Tw', y=.982, fontsize=15, fontweight='bold')
    fig.text(.5, .946, '原生 exit 134，最终原因未知  |  完整场仅 L0 / S1502  |  有限历史至 L56 / S1558',
             ha='center', fontsize=10.5, color='#7c3c2a')
    fig.text(.5, .922, '4 个已恢复节点事件不能定位最终失败；本图不授予收敛、实验或区域物理资格',
             ha='center', fontsize=9.5, color='#444444')
    upper = axes[0].secondary_xaxis('top', functions=(lambda value: value+offset, lambda value: value-offset))
    upper.set_xticks([1502, 1512, 1522, 1532, 1542, 1552, 1558])
    upper.set_xlabel('系列更新编号 S = 1502 + L（冷重启计数）', labelpad=7)

    residuals = [('rms[Rho]', '#24698f', '密度'), ('rms[RhoU]', '#618b70', '动量 x'),
                 ('rms[RhoV]', '#b69937', '动量 y'), ('rms[RhoW]', '#8f73ac', '动量 z'),
                 ('rms[RhoE]', '#bc593e', '总能量')]
    for key, color, label in residuals:
        axes[0].plot(x, data[key], color=color, lw=1.5, label=label)
    axes[0].axhline(-8, color='.35', ls=':', lw=1, label='原 −8 门槛')
    axes[0].set_ylim(-8.6, 3.35)
    axes[0].set_ylabel('log10 RMS 残差')
    axes[0].set_title('a  原生有限历史：残差门槛未过', loc='left', fontsize=10)
    axes[0].legend(ncol=6, loc='lower left', bbox_to_anchor=(0, .035), fontsize=8.2, frameon=False)

    axes[1].fill_between(x, data['Min CFL'], data['Max CFL'], color='#2c7296', alpha=.17,
                         label='原生局部 CFL 范围')
    axes[1].plot(x, data['Avg CFL'], color='#24698f', lw=1.8, label='原生平均 CFL')
    axes[1].axhline(.4, color='.45', ls=':', lw=1, label='初始配置 CFL = 0.4')
    axes[1].set_yscale('log'); axes[1].set_ylabel('下一步 CFL（对数轴）')
    axes[1].set_title('b  实际 CFL 反馈；θ 未反馈给原生 alpha / CFL', loc='left', fontsize=10)
    axes[1].legend(ncol=3, loc='upper left', fontsize=8.3, frameon=False)

    axes[2].plot(x, data['HF'], color='#bc593e', lw=1.8, label='历史 HF，保留正负号')
    axes[2].fill_between(x, data['HF'], 0, where=data['HF'] >= 0, color='#bc593e', alpha=.10)
    axes[2].fill_between(x, data['HF'], 0, where=data['HF'] < 0, color='#24698f', alpha=.08)
    axes[2].axhline(0, color='.35', lw=.85)
    axes[2].set_ylabel('净入壁热功率 / W')
    axes[2].set_title('c  四分之一计算域 signed HF；正号入壁', loc='left', fontsize=10)
    axes[2].legend(loc='upper left', fontsize=8.3, frameon=False)

    for event in event_rows:
        local, number = event['localIteration'], event['event']
        axes[3].scatter([local], [number], color='#a15e3d', marker='D', s=34, zorder=5)
        label = ('L%d / S%d · θ=0.5 · RCM %d / 全局 %d' %
                 (local, event['seriesUpdateNumber'], event['localRCMPoint'], event['globalPoint']))
        axes[3].annotate(label, (local, number), xytext=(-9, 0), textcoords='offset points',
                         ha='right', va='center', fontsize=8.5, color='#5b3e2c',
                         bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=1.5))
    axes[3].set_yticks([1, 2, 3, 4], ['E1', 'E2', 'E3', 'E4'])
    axes[3].set_ylim(4.6, .4)
    axes[3].set_title('d  4 次已恢复 node.trial / rhoe 事件：旧态均合法；不推断最终退出原因',
                      loc='left', fontsize=10)
    axes[3].set_ylabel('事件编号'); axes[3].set_xlabel('本段局部迭代编号 L')
    axes[3].set_xticks([0, 10, 20, 30, 40, 46, 50, 56])
    for ax in axes:
        for local in (46, 56):
            ax.axvline(local, color='#b57c5d', ls='--', lw=.85, alpha=.58, zorder=0)
        ax.grid(alpha=.16); ax.set_xlim(-.6, 57.5)
    fig.text(.105, .041, '仅使用 57 行 history、4 条事件及小型回执；未读取网格、重启或体表场。\n'
             '事件为 Output=true 阶段记录；L46 / S1548、L56 / S1558 均无完整保存场。',
             fontsize=9, color='#444444', va='bottom')
    figures = [output/'spatial_step_history.png', output/'spatial_step_history.svg']
    fig.savefig(str(figures[0]), dpi=170)
    fig.savefig(str(figures[1]))
    plt.close(fig)
    stable = {name: sha(Path(name)) == digest for name, digest in locks.items()}
    assert all(stable.values())
    crossings = [dict(beforeLocal=int(x[i-1]), afterLocal=int(x[i])) for i in range(1, len(x))
                 if data['HF'][i-1]*data['HF'][i] < 0]
    report = dict(passed=True, scope='scientific-small-history-plot-only',
        figures=[dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path)) for path in figures],
        inputSHA256=locks, inputsStable=stable, rows=len(x), nativeExitCode=receipt['exitCode'],
        nativeElapsedSeconds=receipt['elapsedSeconds'], terminationActualCauseKnown=False,
        completeSavedLocalIterations=[0], completeSavedSeriesUpdates=[1502],
        lastFiniteHistoryLocalIteration=int(x[-1]), lastFiniteHistorySeriesUpdates=int(series[-1]),
        localToSeriesUpdateFormula='series = priorSuccessfulUpdates + local + 1 = 1502 + local',
        eventRows=event_rows, recoveredEventsIdentifyFinalFailure=False,
        finalLogResidual={key: float(data[key][-1]) for key, _, _ in residuals},
        actualCFL=dict(meanFirst=float(data['Avg CFL'][0]), meanLast=float(data['Avg CFL'][-1]),
            meanMinimum=float(data['Avg CFL'].min()), meanMaximum=float(data['Avg CFL'].max()),
            localMinimum=float(data['Min CFL'].min()), localMaximum=float(data['Max CFL'].max())),
        signedHeatPowerW=dict(domain='original-quarter-domain', positiveDirection='into-wall',
            first=float(data['HF'][0]), last=float(data['HF'][-1]), minimum=float(data['HF'].min()),
            maximum=float(data['HF'].max()), signChangeBrackets=crossings),
        largeFieldsRead=0, CFDLaunches=0, builds=0, physicalQualified=False,
        eventTemperatureDirectionClassified=False, uninterruptedTrajectoryEquivalent=False,
        elapsedSeconds=time.monotonic()-started)
    save(output/'report.json', report)
    print(dict(figures=[str(path) for path in figures], rows=len(x), events=len(event_rows),
               elapsedSeconds=report['elapsedSeconds']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--termination-review', required=True)
    args = parser.parse_args()
    plot(Path(args.folder).resolve(), Path(args.output).resolve(), Path(args.termination_review).resolve())
