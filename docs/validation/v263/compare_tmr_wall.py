"""在 NASA 风洞模型的原测站比较已保存的内部解，同时导出完整截断固壁。"""
import argparse
import csv
import hashlib
import json
import math
from decimal import Decimal
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'tests/data/tmr_mach7'
RUN = ROOT / 'outputs/v67_validation'
HASHES = {
    'index': '657eaf9b521551e2bcd64c631b8022ea76ce6e1a12203854e469b690add84351',
    'contract': '121b83d627f2ac573ea38dbcb155e177e321d2319c4f379474f4698d92935b19',
    'p': 'd39a6c5a82bb3dc1cf9d2151e635a3ee28d9799ea6d04650b00c6af0b0acd116',
    'q': 'be3e93170746b390418583d258be486e2d5bdfc2e39b30d40b6ca045f9e723f1',
}


def check(ok, message):
    if not ok:
        raise ValueError(message)


def locked(path, digest):
    content = path.read_bytes()
    check(hashlib.sha256(content).hexdigest() == digest, str(path))
    return content


def measurements(file, key, reference):
    # 只合并原弧长和原数值均相同的重复行，保留可追溯行号。
    points = {}
    count = 0
    for line in locked(file, HASHES[key]).decode().splitlines():
        if not line.strip() or line.lstrip()[0] not in '-+0123456789.':
            continue
        count += 1
        s_cm, ratio = map(Decimal, line.split())
        record = points.setdefault(s_cm, {'value': float(ratio * reference), 'rows': []})
        check(record['value'] == float(ratio * reference), f'conflicting repeats at {s_cm} cm')
        record['rows'].append(count)
    return points


def wall_from_frame(grid, archive):
    run = archive['runs'][grid]
    path = (ROOT / 'outputs/v66_validation/tmr_wall_thermal/tmr_wall_thermal_evidence.mat'
            if grid == 'baseline' else RUN / 'both/tmr_sensitivity_evidence.mat')
    locked(path, run['sha256'])
    for receipt in run['actualSuccessfulProcesses']:
        rec_path = path.parents[1] / Path(receipt['path']).name
        obs = json.loads(locked(rec_path, receipt['sha256']))
        check(obs['status'] == 'terminal' and obs['exitCode'] == 0 and all(obs['sourceStableAfter'].values()), 'source process')
    saved = loadmat(path, simplify_cells=True, variable_names=['mesh','implicit'])
    result = saved['implicit']
    check(result['steadyConverged'] and max(result['history']['residual'][-3:]) <= 2e-6
          and result['globalBalanceRelative'] <= 5e-10 and not result['physicalSolution'], 'saved steady state')
    mesh = saved['mesh']
    x, r = mesh['X'][:, 0], mesh['R'][:, 0]
    check(np.all(np.diff(x) < 0), 'ambiguous wall ordering')
    w = result['finalInfo']['wall']
    order = np.argsort(-w['centroid'][:, 0])
    check(len(order) == len(x)-1 and len(np.unique(w['faceID'])) == len(order), 'wall face count')
    check(np.all((w['centroid'][order, 0] < x[:-1]) & (w['centroid'][order, 0] > x[1:])), 'wall face mapping')
    a = np.pi * (r[:-1] + r[1:]) * np.hypot(np.diff(x), np.diff(r))
    np.testing.assert_allclose(a, 2*np.pi*w['areaPerRadian'][order], rtol=2e-12, atol=1e-14)
    p, q, tw = (np.asarray(w[key])[order] for key in ('pressure','heatFluxIntoWall','temperatureInput'))
    check(np.all(np.isfinite(p)) and np.all(np.isfinite(q)) and np.all(tw == 311), 'wall fields')
    rows = [dict(segment=i+1, native_face_id=int(w['faceID'][j]), X1_m=x[i], R1_m=r[i],
                 X2_m=x[i+1], R2_m=r[i+1], ring_area_m2=a[i], pressure_Pa=p[i],
                 heat_into_wall_W_m2=q[i], given_Tw_K=tw[i], solved_solid_Tw_K=None)
            for i, j in enumerate(order)]
    return x, rows


def write_csv(path, rows):
    with path.open('x', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare(out):
    out.mkdir(parents=True, exist_ok=False)
    archive = json.loads(locked(RUN/'comparison/sensitivity_comparison.json', HASHES['index']))
    c = json.loads(locked(DATA/'contract.json', HASHES['contract']))
    ref = c['experimental_local_reference']
    raw = {k: measurements(DATA/'raw'/file, k, Decimal(str(ref[name]))) for k, file, name in
           [('p','pw_exp.dat','pressure_Pa'), ('q','qw_exp.dat','reported_wall_heating_magnitude_W_m2')]}
    limits = {k: c['measurement_uncertainties'][name] for k, name in
              [('p','surface_pressure_relative_estimate'), ('q','surface_heat_relative_estimate')]}
    check(len(raw['p']) == 24 and len(raw['q']) == 27, 'source coverage')
    paired = sorted(set(raw['p']) & set(raw['q']))
    union = sorted(set(raw['p']) | set(raw['q']))
    study = {}
    table = []
    for grid in ('baseline','both'):
        x, wall = wall_from_frame(grid, archive)
        old = {(r['field'], r['sourceRow']): r for r in archive['wallPoints'][grid]}
        counters = {'p':0, 'q':0, 'both':0}
        for s in union:
            axial = -float(s/100) * (math.cos(math.radians(c['scope']['flare_angle_deg'])) if s > 0 else 1)
            hits = np.flatnonzero((axial <= x[:-1]) & (axial > x[1:]))
            check(len(hits) == 1, 'station outside wall or ambiguous')
            face = wall[int(hits[0])]
            row = {'grid':grid, 'S_cm':str(s), 'native_face_id':face['native_face_id']}
            for k, variable, field in [('p','pressure','pressure_Pa'),('q','heating','heat_into_wall_W_m2')]:
                entry = raw[k].get(s)
                value = face[field]
                error = (value / entry['value'] - 1) if entry else None
                match = abs(error) <= limits[k] if entry else None
                row.update({k+'_exp_rows':'|'.join(map(str,entry['rows'])) if entry else '',
                            'measured_'+field:entry['value'] if entry else None, 'calculated_'+field:value,
                            k+'_relative_error':error, k+'_within_estimate':match})
                if entry:
                    for source_id in entry['rows']:
                        saved = old[variable,source_id]
                        check(saved['faceID'] == face['native_face_id'] and
                              math.isclose(saved['experimentalValueSI'],entry['value'],rel_tol=2e-12) and
                              math.isclose(saved['computedSignedValueSI'],value,rel_tol=2e-12), 'independent sample mismatch')
                    if s in paired:
                        counters[k] += int(match)
            row['joint_within_estimates'] = bool(row['p_within_estimate'] and row['q_within_estimate']) if s in paired else None
            counters['both'] += int(row['joint_within_estimates'] is True)
            table.append(row)
        study[grid] = dict(gridCells=archive['summaries'][grid]['cellCount'], rings=len(wall),
                           fullTruncatedWallAreaM2=math.fsum(r['ring_area_m2'] for r in wall),
                           insideEstimatesAtPairedStations=counters, sourceMAT_SHA256=archive['runs'][grid]['sha256'])
        if grid == 'both':
            write_csv(out/'wind_tunnel_wall.csv', wall)
    write_csv(out/'wall_stations.csv', table)
    common = {g:{r['S_cm'] for r in table if r['grid']==g and r['joint_within_estimates'] is True} for g in study}
    result = dict(case=c['case_id'], wallTemperatureGivenK=311, indexedPressureRows=sum(len(r['rows']) for r in raw['p'].values()),
                  uniquePressureStations=len(raw['p']), heatStations=len(raw['q']), unionStations=len(union),
                  commonNominalStations=len(paired), commonUpstream=sum(s<0 for s in paired),
                  withinBothEstimatesOnBothGrids=len(common['baseline'] & common['both']),
                  pressureRelativeEstimate=limits['p'], heatRelativeEstimate=limits['q'],
                  method='axisymmetric SA-noft2, Sutherland transport, Dirichlet tangent wall gradient',
                  grids=study, independentPhysicalAreaAddedM2=0, newSolverRuns=0, sources=HASHES)
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    # 同站误差图只画实测重合点，空白坐标不补点。
    fig, axs = plt.subplots(2, 1, figsize=(8.5,5.7), sharex=True, constrained_layout=True)
    for ax, key, title in zip(axs, ['p','q'], ['Wall pressure', 'Heat flux into wall']):
        band=100*limits[key]
        ax.axhspan(-band,band,color='#dce8e4',label='Reported measurement estimate')
        ax.axhline(0,color='#6a7777',linewidth=.8)
        for grid, marker, color in [('baseline','o','#364c6a'),('both','+','#ad532a')]:
            rows=[r for r in table if r['grid']==grid and r['joint_within_estimates'] is not None]
            ax.scatter([float(r['S_cm']) for r in rows],[100*r[key+'_relative_error'] for r in rows],
                       color=color, marker=marker, s=38, label=f"{study[grid]['gridCells']:,} cells")
        ax.set_title(title,loc='left',fontsize=11)
        ax.set_ylabel('Relative error (%)')
        ax.spines[['top','right']].set_visible(False)
    axs[0].legend(fontsize=8,loc='upper right')
    axs[-1].set_xlabel('Surface distance from flare junction, S (cm)')
    fig.suptitle('NASA ACC-SBLI | Experimental stations with both measurements',fontsize=13)
    fig.savefig(out/'paired_errors.png',dpi=155)
    plt.close(fig)
    print(json.dumps(result,ensure_ascii=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('new_output_directory', type=Path)
    compare(parser.parse_args().new_output_directory)
