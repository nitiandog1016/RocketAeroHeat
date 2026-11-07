"""独立审核真四面体NK控制；正常发展、数值收敛和物理资格分别记录。"""
from pathlib import Path
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import sha, need, read_field, read_grid, compare_source
from prepareHybridRocketEulerNative import local, save
from eilmerAirThermo import evaluate, bind_model
from hybridEdgeSampling import map_original_wall
from hybridEulerReference import control_reference

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/v197_validation'


def read_diagnostics(path, final_step):
    """模块一：按实际列名读取NK诊断；不借用显式时间推进的列或残差尺度。"""
    names = Path(path).read_text().splitlines()[0].split()
    need(len(names) == len(set(names)), 'NK诊断列重复')
    a = np.loadtxt(path, skiprows=1, ndmin=2)
    need(a.shape == (final_step, len(names)) and np.isfinite(a).all(), 'NK诊断缺行或非有限')
    d = {n: a[:, i] for i, n in enumerate(names)}
    need(np.array_equal(d['step'], np.arange(1, final_step+1)), 'NK步号不连续')
    need(np.all(d['dt'] > 0) and np.all(d['CFL'] > 0), 'NK伪时间非正')
    need(np.all((d['relaxation-factor'] > 0) & (d['relaxation-factor'] <= 1)), 'NK更新松弛错误')
    return d


def snapshot_steps(case, history):
    """模块二：固定源码的restart-metadata行序对应0001起的发展帧；0000是初场。"""
    base = Path(case)/'lmrsim/snapshots'
    data = np.loadtxt(base/'restart-metadata', ndmin=2)
    need(data.shape[1] == 12, 'NK重启元数据列合同改变')
    steps = data[:, 0].astype(int)
    need(np.array_equal(data[:, 0], steps) and np.all(np.diff(steps) > 0), 'NK保存步号错误')
    need(steps[0] > 0 and steps[-1] <= len(history['step']), 'NK保存步超诊断范围')
    need([p.name for p in sorted(base.glob('[0-9][0-9][0-9][0-9]'))] ==
         ['%04d' % i for i in range(len(steps)+1)], 'NK完整帧与元数据不相容')
    # 同一保存步的dt/CFL及五分量残差逐项交叉核对，避免只按目录猜步号。
    cols = ['dt', 'CFL', 'global-residual-abs', 'mass-abs', 'x-mom-abs',
            'y-mom-abs', 'z-mom-abs', 'total-energy-abs']
    for name, index in zip(cols, [1, 2, 5, 7, 8, 9, 10, 11]):
        expected = history[name][steps-1]
        need(np.allclose(data[:, index], expected, rtol=2e-13, atol=0), 'NK保存步诊断不一致: '+name)
    return [0]+steps.tolist()


def normal_terminal(process, text, final_step):
    """模块三：退出、清点、准确复数二进制日志分别验证；exit0不代表收敛。"""
    need(process['passed'] and process['nativeExitCode'] == 0 and process['inputsStable'] and
         process['stopReason'] is None and not process['liveGroupAfterCleanup'], 'NK非正常完整终态')
    need(re.findall(r'^Number-type: (.+)', text, re.M) == ['complex'], '非原生复数NK')
    need(re.findall(r'^FINAL-STEP: (\d+)', text, re.M) == [str(final_step)], 'NK最终步不符')
    need(re.findall(r'^STOP-REASON: (.+)', text, re.M) == ['maximum-steps'], 'NK停止原因不符')


def audit():
    """模块四：完整场同EOS、原控制壁精确映射与事前收益；仅只读已有唯一计算。"""
    start = time.monotonic()
    out = BASE/'audit01'; out.mkdir(exist_ok=False)
    shutil.copyfile(__file__, out/'executed_entry.py')
    case = BASE/'tetra01'; rp = BASE/'run01/report.json'
    record = json.loads(rp.read_text()); process = record['process']
    sources = [Path(__file__), rp, BASE/'preparation01/report.json',
               ROOT/'docs/validation/v197/ENTRY_DECISION.md']
    hashes = {str(p): sha(p) for p in sources}
    r = dict(passed=False, physicalQualification=False, newSolvedHeatAreaM2=0)
    try:
        need(record['processPassed'] and record['profile'] == 'tetra-cea-compression-expansion-nk-120-v1', '控制身份不符')
        need(local(process['cwd']).resolve() == case.resolve(), '运行目录错配')
        for p, digest in record['outputSHA256'].items():
            need(sha(local(p)) == digest, '原生输出改变: '+p)
        pkg = json.loads((ROOT/'tools/runtime/eilmer-nk-geometric-v1/runtime.json').read_text())
        need(process['inputSHA256'][pkg['binary']] == pkg['filesSHA256'][pkg['binary']], '实际执行文件来源错误')
        text = (rp.parent/'native/stdout.txt').read_text()
        normal_terminal(process, text, 120)
        hist = read_diagnostics(case/'lmrsim/diagnostics/nk-diagnostics', 120)
        steps = snapshot_steps(case, hist)
        need(steps == list(range(0, 121, 20)), '缺少事前计划保存帧')
        grid = read_grid(case/'lmrsim/grid/grid-0000'); mapping, _ = compare_source(grid, case/'mesh.su2')
        need(len(grid['cells']) == 14400 and np.all(grid['kinds'] == 10), '非实际14400四面体控制')
        surface = json.loads((case/'surface.json').read_text())
        bind_model(case/'rah-air.gas', sha(case/'rah-air.gas'))
        tri = np.array(surface['vertices'])[np.array(surface['facesOneBased'])-1]
        areas = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)/2
        centers = tri.mean(1); frames = []; walls = []
        for snapshot, step in enumerate(steps):
            f = read_field(case, snapshot, 14400)
            need(all(np.all(f[k] > 0) for k in ['p', 'rho', 'T', 'vol']), '保存场非正')
            need(np.all(f['massf-air'] == 1), '空气组分改变')
            gas = evaluate(f['p'], f['T'])
            errors = {name: float(np.max(abs(f[name]-gas[key])/np.maximum(1, abs(gas[key]))))
                      for name, key in [('rho', 'densityKgm3'), ('e', 'internalEnergyJkg'), ('a', 'soundSpeedMps')]}
            need(max(errors.values()) <= 2e-10, '保存场EOS不相容')
            positions = np.column_stack([f['pos.'+a] for a in 'xyz'])
            location = map_original_wall(surface, grid, mapping['wall'], positions)
            owners = np.array(location['nativeOwnerCellsZeroBased'])
            frames.append(dict(snapshot=snapshot, step=step, temperatureRangeK=[float(f['T'].min()), float(f['T'].max())],
                               minimumPressurePa=float(f['p'].min()), minimumDensityKgm3=float(f['rho'].min()), EOSError=errors))
            walls.append(np.column_stack([f[k][owners] for k in ['p', 'T', 'vel.x', 'vel.y', 'vel.z']]))
        improvement = {key: float(np.log10(hist[col][:20]).mean()-np.log10(hist[col][-20:]).mean())
                       for key, col in [('mass', 'mass-abs'), ('energy', 'total-energy-abs')]}
        relative = {col: float(np.max(hist[col][-20:])) for col in
                    ['mass-rel', 'x-mom-rel', 'y-mom-rel', 'z-mom-rel', 'total-energy-rel']}
        changes = {key: float(np.sum(abs(walls[-1][:, j]-walls[-2][:, j])*areas)/np.sum(abs(walls[-1][:, j])*areas))
                   for j, key in enumerate(['pressure', 'temperature'])}
        ref = control_reference(); save(out/'same_eos_reference.json', ref)
        comparisons = []
        for region, lo, hi, reference in [('compression', .35, .5, ref['states'][1]), ('expansion', .85, .95, ref['states'][2])]:
            use = (-centers[:, 0] >= lo) & (-centers[:, 0] <= hi)
            values = np.average(walls[-1][use], axis=0, weights=areas[use])
            comparisons.append(dict(region=region, count=int(use.sum()), pressurePa=float(values[0]), temperatureK=float(values[1]),
                pressureRelativeError=float(values[0]/reference['pressurePa']-1),
                temperatureRelativeError=float(values[1]/reference['temperatureK']-1),
                status='unconverged development observation; not accuracy qualification'))
        np.savetxt(out/'wall_development.csv', np.column_stack([surface['faceIDs'], centers, areas, walls[-1]]),
            delimiter=',', comments='', header='faceID,x_m,y_m,z_m,area_m2,p_Pa,T_K,ux_mps,uy_mps,uz_mps', fmt='%.17g')
        r.update(passed=True, screeningPassed=all(x >= .05 for x in improvement.values()),
            savedSteps=steps, fields=frames, controlWallFaces=len(areas), controlWallAreaM2=float(areas.sum()),
            meanLog10ImprovementLast20VsFirst20=improvement, minimumRequiredImprovement=.05,
            nativeRelativeResidualLast20Max=relative, nativeRelativeResidualGatePassed=max(relative.values()) <= 1e-8,
            residualScope='NK native component ratios to its initial reference; not explicit solver residual normalization',
            nativeGlobalRelativeResidualFinal=float(hist['global-residual-rel'][-1]),
            numericalConvergencePassed=False, lastPairWallAreaL1=changes, platformComparisons=comparisons,
            nativeMassBalanceFinal=float(hist['mass-balance'][-1]),
            finalCFL=float(hist['CFL'][-1]), relaxationFactorRange=[float(hist['relaxation-factor'].min()), float(hist['relaxation-factor'].max())],
            observedLinearSolveTargetMet=bool(np.all(hist['linear-solve-residual-rel-achieved'] <= hist['linear-solve-residual-rel-target'])),
            processSeconds=process['elapsedSeconds'], peakGroupRSSBytes=process['peakObservedGroupRSSBytes'],
            stableInputCount=len(process['inputStable']), wallHeatFluxIdentity='unknown, not zero', solvedWallTemperatureIdentity='unknown')
    except Exception as error:
        r['error'] = repr(error)
    r.update(sourceSHA256=hashes, inputsStable=all(sha(p) == h for p, h in hashes.items()), elapsedSeconds=time.monotonic()-start)
    save(out/'report.json', r)
    print(json.dumps({k: r.get(k) for k in ['passed', 'screeningPassed', 'error', 'meanLog10ImprovementLast20VsFirst20', 'elapsedSeconds']}))
    return r


if __name__ == '__main__':
    result = audit()
    raise SystemExit(0 if result['passed'] else 1)
