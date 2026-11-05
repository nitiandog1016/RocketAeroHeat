"""只读已保存公开场，记录失败前可观察范围及上游边界层，不推定故障因果。"""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts'))
from su2AxisymmetricAnalysis import mesh_data, fields, wall_values, sample_stations, integrate_scalar
from su2AxisymmetricSurface import surface_frame, compare_volume
from su2BenchmarkAcceptance import stability


def inspect():
    folder = ROOT/'outputs/v152_validation/acc_sbli01'
    manifest = json.loads((folder/'study_manifest.json').read_text(encoding='utf8'))
    contract = json.loads((ROOT/'tests/data/tmr_mach7/contract.json').read_text(encoding='utf8'))
    p, cells, markers = mesh_data(folder/'mesh.su2')
    c = manifest['benchmark']['conditions']
    cp = c['gamma']*c['gasConstant']/(c['gamma']-1)
    h_inf = cp*c['temperatureK']+c['speedMS']**2/2
    volumes = []
    for path in sorted(folder.glob('restart_flow_*.csv')):
        f = fields(path, True)
        wall = wall_values(p, markers, f)
        stations = sample_stations(wall, contract)
        h = cp*f['Temperature']+np.sum(f['velocity']**2, axis=1)/2
        ratio = f['Eddy_Viscosity']/f['Laminar_Viscosity']
        # 在源x=-0.06处逐径向层插值，统计直接来自原网格的边界层状态。
        profile = []
        for row in np.arange(81*101).reshape(101, 81):
            ids = row[np.argsort(p[row, 0])]
            profile.append([np.interp(.06, p[ids, 0], values[ids]) for values in
                            (p[:, 1], -f['velocity'][:, 0], ratio, f['Turb_Kin_Energy'], f['Omega'])])
        profile = np.array(profile)
        peak = int(np.argmax(h))
        volumes.append(dict(localIteration=int(path.stem.split('_')[-1]),
            powerW=integrate_scalar(p, markers['wall'], f['Heat_Flux']),
            pressureRMS=stations['pressure']['rmsRelativeError'], heatRMS=stations['heat']['rmsRelativeError'],
            upstreamPressurePa=float(np.interp(-.06, wall['S'], wall['p'])),
            upstreamHeatWm2=float(np.interp(-.06, wall['S'], wall['q'])),
            yPlus99=float(np.percentile(wall['yplus'], 99)),
            maximumPrimitiveH0Ratio=float(h[peak]/h_inf), maximumPrimitiveH0Node=peak,
            maximumPrimitiveH0Coordinate=p[peak].tolist(),
            maximumEddyViscosityRatio=float(max(ratio)),
            upstreamMaximumEddyViscosityRatio=float(max(profile[:, 2])),
            upstreamProfile=profile.tolist()))
    frames, paired = [], []
    for path in sorted(folder.glob('surface_flow_*.vtk')):
        frame = surface_frame(path, p, markers, contract)
        volume = folder/('restart_flow_%06d.csv' % frame['iteration'])
        if volume.exists():
            paired.append(dict(localIteration=frame['iteration'], errors=compare_volume(frame, volume, p, markers)))
        frames.append(frame)
    last_surface = frames[-1]['iteration']
    steady = stability(frames, last_surface, manifest['benchmark']['geometry']['revolvedWallAreaM2'])
    return dict(scope='saved states only; failure state unavailable; no causal attribution',
                primitiveH0Definition='cp*T_saved + velocity_saved^2/2; excludes SST k',
                upstreamProfileColumns=['radius_m', 'downstream_velocity_m_s', 'mu_t_over_mu', 'k_m2_s2', 'omega_s_inv'],
                volumes=volumes, surfaceVolumePairing=paired,
                lastSurfaceIteration=last_surface,
                lastSurfaceWindowDiagnostic=steady,
                lastSurfaceWindowIsTerminalAcceptance=False,
                surfacePowerHistory=[dict(localIteration=f['iteration'], powerW=f['power']) for f in frames])


if __name__ == '__main__':
    report = inspect()
    target = Path(__file__).with_name('saved_state_diagnostics.json')
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    for row in report['volumes']:
        print(json.dumps({k: v for k, v in row.items() if k != 'upstreamProfile'}))
    print(json.dumps(report['lastSurfaceWindowDiagnostic']))
