"""复用收敛细网格与粗网格，固定中段做网格敏感性筛选，不授予整箭资格。"""
from pathlib import Path
import json
import hashlib
import numpy as np
from auditEmbeddedSU2FlatPlate import read_wall

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    coarse = ROOT/'outputs/v105_validation/plate03'
    fine = ROOT/'outputs/v108_validation/plate_warm01'
    files = [Path(__file__).resolve()]
    walls = []
    for folder in (coarse, fine):
        receipt = json.loads((folder/'engine_process.json').read_text(encoding='utf8'))
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and not receipt['timedOut']
        assert all(receipt['inputsStable'].values())
        for name in ('history.csv', 'surface_flow.vtk', 'engine_stdout.log'):
            assert sha(folder/name) == receipt['outputSHA256'][name]
            files.append(folder/name)
        history = np.genfromtxt(str(folder/'history.csv'), delimiter=',', names=True)
        assert history['rmsRho'][-1] < -8 and history['rmsRhoE'][-1] < -8
        assert 'All convergence criteria satisfied.' in (folder/'engine_stdout.log').read_text()
        xyz, wall = read_wall(folder/'surface_flow.vtk')
        order = np.argsort(-xyz[:, 0]);s = -xyz[order, 0]
        assert np.all(np.diff(s) > 0)
        walls.append((s, {k: v[order] for k, v in wall.items()}))
        files += [folder/'engine_process.json', folder/'case.cfg']
    s0, w0 = walls[0];s1, w1 = walls[1]
    # 使用原粗网格中段节点，细场线性插值；未裁去前缘物理区域，仅限定比较的有效区间。
    mask = (s0 >= .04) & (s0 <= .16);sample = s0[mask]
    assert min(sample) >= min(s1) and max(sample) <= max(s1)
    changes = {};columns = [sample]
    for name in ('Pressure', 'Heat_Flux'):
        old = w0[name][mask];new = np.interp(sample, s1, w1[name]);delta = new/old-1
        changes[name] = dict(maxAbsoluteRelative=float(max(abs(delta))),
                             minimumRelative=float(min(delta)), maximumRelative=float(max(delta)))
        columns += [old, new, delta]
    pressure_pass = changes['Pressure']['maxAbsoluteRelative'] <= .01
    heat_pass = changes['Heat_Flux']['maxAbsoluteRelative'] <= .05
    physical = json.loads((fine/'physical_audit.json').read_text(encoding='utf8'));files.append(fine/'physical_audit.json')
    assert physical['residualConverged']
    report = dict(status='converged-two-grid-interior-screen; not-grid-independence',
                  interval=[.04, .16], comparisonNodes=len(sample), changes=changes,
                  thresholds=dict(pressure=.01, heatFlux=.05),
                  pressureScreenPassed=bool(pressure_pass), heatFluxScreenPassed=bool(heat_pass),
                  sampleRule='coarse interior nodes; linear fine-grid interpolation; no extrapolation',
                  boundaryLayerRelativeDifferenceRange=[physical['boundaryLayerRelativeDifferenceMin'],
                                                        physical['boundaryLayerRelativeDifferenceMax']],
                  leadingEdgeExcludedFromComparisonOnly=True, gridIndependenceQualified=False,
                  experimentCompared=False, newQualifiedRocketArea=0,
                  inputSHA256={p.relative_to(ROOT).as_posix(): sha(p) for p in files})
    dest = fine/'mesh_comparison.json';assert not dest.exists()
    dest.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf8'))
    np.savetxt(str(fine/'mesh_comparison.csv'), np.column_stack(columns), delimiter=',',
               header='s_m,p_coarse,p_fine_interpolated,p_relative,q_coarse,q_fine_interpolated,q_relative', comments='')
    print(json.dumps(report, indent=2))
    assert pressure_pass and heat_pass, '保留实际失败，不放宽预定筛选门槛'


if __name__ == '__main__':
    main()
