"""v168模块用途：有限只读复审既存壁温、上游热流与阶段来源，不推进流场。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat

from auditSU2IdealAirState import load_state
from idealAirProperties import IdealAir
from su2SSTInterface import wall_geometry
from su2ThermalRuntime import ROOT, read
from validateIdealAir import save

BASE = ROOT/'outputs/v168_validation'
PUBLIC = ROOT/'outputs/v159_validation/public01'
ROCKET = ROOT/'outputs/v167_validation/real01'


def digest(path):
    """模块一：分块哈希避免将大场副本同时装入内存。"""
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def table(path, rows):
    with path.open('w', encoding='utf8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(output_name='review01'):
    """模块二：绑定历史证据，只产生新目录，不覆盖历史硬链接。"""
    start = time.monotonic()
    assert Path(output_name).name == output_name
    out = BASE/output_name
    out.mkdir(parents=True, exist_ok=False)
    locks = {}

    def lock(path, expected=None):
        path = Path(path).resolve()
        actual = digest(path)
        if expected is not None:
            assert actual == expected, str(path)
        locks[str(path)] = actual
        return path

    def archived(path, version):
        index = read(lock(ROOT/('docs/validation/v%d/evidence_manifest.json' % version)))
        name = Path(path).relative_to(ROOT).as_posix()
        entry = next(item for item in index['files'] if item['source'] == name)
        return lock(path, entry['sha256'])

    deps = [Path(__file__), ROOT/'docs/validation/v168/PLAN.md',
            ROOT/'scripts/auditSU2IdealAirState.py', ROOT/'scripts/idealAirProperties.py',
            ROOT/'scripts/su2SSTInterface.py', ROOT/'scripts/su2ThermalRuntime.py',
            ROOT/'scripts/validateIdealAir.py']
    for path in deps:
        path = lock(path)
        target = out/'executed_sources'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())

    receipt = read(archived(PUBLIC/'engine_process.json', 159))
    old_report = read(archived(PUBLIC/'public_review02/report.json', 159))
    final_review = read(archived(ROOT/'outputs/v159_validation/final_review/report.json', 159))
    assert final_review['passed'] and final_review['originalReportSHA256'] == digest(PUBLIC/'public_review02/report.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
    assert len(receipt['checkpoints']['complete']) == 13
    for name in ('case.cfg', 'case_manifest.json', 'mesh.su2'):
        lock(PUBLIC/name, receipt['inputSHA256'][str(PUBLIC/name)])
    entries = receipt['checkpoints']['complete']
    for entry in entries:
        lock(PUBLIC/entry['volume'], entry['volumeSHA256'])
        lock(PUBLIC/entry['surface'], entry['surfaceSHA256'])

    # 模块三：复用13帧公开剖面，严格区分壁节点极值与插值测站。
    public_rows = []
    public = []
    for entry in entries:
        data = read(archived(PUBLIC/'public_review02'/('frame_%06d.json' % entry['localIteration']), 159))
        w = {key: np.asarray(value) for key, value in data['wall'].items()}
        err = np.abs(w['fluidT']-311.)
        im, iy = int(np.argmax(err)), int(np.argmax(w['yplus']))
        # 轴对称线性壁面上的精确两点高斯积分：r线性，温差平方二次。
        r, s, e = w['radius'], w['S'], w['fluidT']-311.
        weight = np.diff(s)*2*np.pi
        integral = 0.
        for t in (.5-.5/np.sqrt(3), .5+.5/np.sqrt(3)):
            integral += .5*np.sum(weight*((1-t)*r[:-1]+t*r[1:])*((1-t)*e[:-1]+t*e[1:])**2)
        area = float(np.sum(weight*(r[:-1]+r[1:])/2))
        row = dict(localIteration=entry['localIteration'], updates=entry['successfulUpdates'],
            upstreamX_m=.06, upstreamP_Pa=float(np.interp(-.06, s, w['p'])),
            upstreamQ_Wm2=float(np.interp(-.06, s, w['q'])),
            upstreamFluidT_K=float(np.interp(-.06, s, w['fluidT'])),
            upstreamYPlus=float(np.interp(-.06, s, w['yplus'])),
            maximumWallTError_K=float(err[im]), maxTErrorGlobalNode=int(w['ids'][im]),
            maxTErrorX_m=float(-w['x'][im]), maxTErrorS_m=float(s[im]),
            maxTErrorQ_Wm2=float(w['q'][im]), maxTErrorYPlus=float(w['yplus'][im]),
            areaRMSWallTError_K=float(np.sqrt(integral/area)), revolvedArea_m2=area,
            maximumYPlus=float(w['yplus'][iy]), maxYPlusGlobalNode=int(w['ids'][iy]),
            maxYPlusX_m=float(-w['x'][iy]), yPlusP99=float(data['yPlus99']),
            upstreamVelocityRatioRMS=float(data['profile']['velocityRatioRMSError']),
            heatPower_W=float(data['revolvedHeatPowerW']))
        assert abs(row['upstreamQ_Wm2']-data['localReference']['heatFluxWm2']) < 1e-9
        assert abs(row['maximumWallTError_K']-data['maximumWallFluidTemperatureDifferenceK']) < 1e-10
        public_rows.append(row)
        public.append(data)
    table(out/'public_wall_evolution.csv', public_rows)

    # 模块四：只读取公开末体场，独立定位法向邻点与实际热梯度。
    geom = wall_geometry(PUBLIC/'mesh.su2')
    _, fields, nd = load_state(PUBLIC/entries[-1]['volume'], True, captured=True, spatial=True)
    assert nd == 2
    w = public[-1]['wall']
    ids = np.asarray(w['ids'], dtype=int)
    for a, b in [('p', 'Pressure'), ('q', 'Heat_Flux'), ('fluidT', 'Temperature'), ('yplus', 'Y_Plus')]:
        assert np.array_equal(np.array(w[a]), fields[b][ids])
    model = IdealAir()
    locks.update(model.locks)
    selected = np.argsort(abs(geom['points'][geom['wallIDs'], 0]-.06))[:2]
    selected = np.unique(np.r_[selected, np.where(geom['wallIDs'] == public_rows[-1]['maxTErrorGlobalNode'])[0]])
    details = []
    for index in selected:
        i, j = int(geom['wallIDs'][index]), int(geom['normalNeighbor'][index])
        p, T = fields['Pressure'][i], fields['Temperature'][i]
        cp = float(model.state('PT', p, T)['cp'][0])
        lam = cp*(fields['Laminar_Viscosity'][i]/.72+fields['Eddy_Viscosity'][i]/.9)
        d = float(geom['thermalDistance'][index])
        q = float(lam*(fields['Temperature'][j]-311.)/d)
        assert abs(q-fields['Heat_Flux'][i])/max(1, abs(q)) < 2e-10
        details.append(dict(globalNode=i, neighborGlobalNode=j, projectX_m=float(geom['points'][i, 0]),
            radius_m=float(geom['points'][i, 1]), givenTw_K=311., fluidT_K=float(T),
            neighborFluidT_K=float(fields['Temperature'][j]), thermalDistance_m=d,
            conductivity_W_mK=float(lam), recomputedQ_Wm2=q, savedQ_Wm2=float(fields['Heat_Flux'][i]),
            yPlus=float(fields['Y_Plus'][i]),
            savedTotalResidualEnergy=float(fields['Residual_Energy'][i]),
            savedTotalResidualDensity=float(fields['Residual_Density'][i]),
            nativeSplitBoundaryRHSAvailable=False, nativeHistoricalIncrementAvailable=False))
    table(out/'public_final_wall_details.csv', details)
    prof = public[-1]['profile']
    table(out/'public_final_upstream_profile.csv', [dict(y_m=y, predictedU_over_Uref=u, measuredU_over_Uref=v)
        for y, u, v in zip(prof['y'], prof['predictedU'], prof['measuredU'])])

    # 模块五：复用三维已审核壁面CSV，按原三角面积核对误差并定位极值。
    rocket_report = read(archived(ROOT/'outputs/v167_validation/review01/report.json', 167))
    rr = read(archived(ROCKET/'engine_process.json', 167))
    geometry_path = ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
    vol = loadmat(str(lock(geometry_path, rocket_report['inputSHA256'][str(geometry_path)])),
                  variable_names=['vertices', 'sourceFaces'], squeeze_me=True)
    F = vol['sourceFaces'].astype(int)-1
    V = vol['vertices'][:9878]
    assert F.shape == (19752, 3)
    area = np.linalg.norm(np.cross(V[F[:, 1]]-V[F[:, 0]], V[F[:, 2]]-V[F[:, 0]]), axis=1)/2
    rocket_rows = []
    for entry, old in zip(rr['checkpoints']['complete'], rocket_report['frames']):
        step = entry['localIteration']
        assert step == old['localIteration']
        audit_dir = ROCKET/('audit%03d' % step)
        audit = read(archived(audit_dir/'report.json', 167))
        path = lock(audit_dir/'wall.csv', audit['wallCSV_SHA256'])
        wall = np.genfromtxt(str(path), delimiter=',', names=True)
        assert np.array_equal(wall['nodeID'], np.arange(len(V)))
        # 旧审核CSV坐标源自7位VTK；按已审核全局ID取精确网格坐标，保留原1微米映射门槛。
        csv_xyz = np.column_stack([wall[a+'_m'] for a in 'xyz'])
        map_error = float(np.max(np.linalg.norm(csv_xyz-V, axis=1)))
        assert map_error < 1e-6
        assert np.all(wall['given_Tw_K'] == 300.)
        T, q, yp = wall['fluid_T_K'], wall['q_W_m2'], wall['yPlus']
        e = T[F]-300.
        rms = float(np.sqrt(area@((np.sum(e, axis=1)**2+np.sum(e**2, axis=1))/12)/sum(area)))
        assert abs(rms-old['areaRMSWallTemperatureErrorK']) < 1e-11
        im, iy = int(np.argmax(abs(T-300))), int(np.argmax(yp))
        row = dict(localIteration=step, seriesUpdates=old['seriesSuccessfulUpdates'],
            maximumWallTError_K=float(abs(T[im]-300)), areaRMSWallTError_K=rms,
            maxTErrorGlobalNode=im, maxTErrorX_m=float(V[im, 0]), maxTErrorY_m=float(V[im, 1]),
            maxTErrorZ_m=float(V[im, 2]), maxTErrorQ_Wm2=float(q[im]), maxTErrorYPlus=float(yp[im]),
            maxYPlusGlobalNode=iy, maxYPlusX_m=float(V[iy, 0]), maximumYPlus=float(yp[iy]),
            yPlusP99=float(np.quantile(yp, .99)), roundedVTKCoordinateMaxDistance_m=map_error,
            faceMeanTErrorOver1KAreaFraction=float(sum(area[abs(np.mean(e, axis=1)) > 1.])/sum(area)),
            faceMeanTErrorOver5KAreaFraction=float(sum(area[abs(np.mean(e, axis=1)) > 5.])/sum(area)))
        assert abs(row['maximumWallTError_K']-old['maximumWallTemperatureErrorK']) < 1e-11
        rocket_rows.append(row)
    table(out/'rocket_wall_evolution.csv', rocket_rows)

    # 模块六：准确编译树中的函数片段，当前文件须哈希相同才能作为旧执行来源。
    source_root = ROOT/'tools/runtime/su2-thermal-build/source'
    build = ROOT/'outputs/v159_validation/build01'
    tree = read(archived(build/'source_tree.json', 159))
    code_records = []
    choices = [('SU2_CFD/src/solvers/CNSSolver.cpp', 630, 792),
               ('SU2_CFD/src/solvers/CNSSolver.cpp', 140, 178),
               ('SU2_CFD/src/variables/CNSVariable.cpp', 140, 194),
               ('SU2_CFD/src/iteration/CFluidIteration.cpp', 70, 155)]
    for name, a, b in choices:
        path = build/'executed_sources'/name
        if not path.exists():
            path = source_root/name
        lock(path, tree[str(source_root/name)])
        lines = path.read_text(encoding='utf8').splitlines()
        target = out/('%s_%d_%d.txt' % (Path(name).stem, a, b))
        target.write_text('\n'.join('%d: %s' % (i+1, lines[i]) for i in range(a-1, b))+'\n', encoding='utf8')
        code_records.append(dict(source=str(path), compiledSHA256=tree[str(source_root/name)],
                                 firstLine=a, lastLine=b, excerpt=target.name))

    # 模块七：新统计不改变原26门槛；输出只作研究判断，不授予物理资格。
    stable = {path: digest(path) == value for path, value in locks.items()}
    assert all(stable.values())
    report = dict(schema='RocketAeroHeat.wall-heat-evidence-review.v1', passed=True,
        publicFrames=13, rocketFrames=3, public=public_rows, rocket=rocket_rows,
        finalPublicDetails=details, exactSourceExcerpts=code_records,
        inheritedPublicAcceptance=final_review['acceptanceWithSpecificPendingReason'],
        experimentalTotalTemperatureDenominatorK=None, newNativeUpdates=0,
        originalThresholdsChanged=False, sourceSHA256=locks, inputsStable=stable,
        elapsedSeconds=time.monotonic()-start, physicsQualified=False, newQualifiedArea_m2=0)
    save(out/'report.json', report)
    size = sum(path.stat().st_size for path in out.rglob('*') if path.is_file())
    assert size <= 16*1024**2
    print(json.dumps(dict(passed=True, elapsedSeconds=report['elapsedSeconds'], outputBytes=size,
                         publicFinal=public_rows[-1], rocketFinal=rocket_rows[-1]), ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-name', default='review01')
    main(parser.parse_args().output_name)
