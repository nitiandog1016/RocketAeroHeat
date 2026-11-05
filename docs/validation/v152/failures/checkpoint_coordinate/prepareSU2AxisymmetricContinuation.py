"""公开域同模型检查点续算，固定系列上限并绑定准确的历史执行源。"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json
from prepareSU2CheckpointContinuation import config
from prepareSU2AxisymmetricInitialState import run_prepared
from su2AxisymmetricAnalysis import mesh_data
from su2AxisymmetricSurface import surface_frame, compare_volume
from su2SSTInterface import state


def series_budget(completed, requested):
    """同一相容初场系列最多12000更新，拒绝段号重置扩展总预算。"""
    assert isinstance(completed, int) and isinstance(requested, int)
    assert completed > 0 and requested > 0 and completed+requested <= 12000
    return dict(completedUpdatesBeforeRun=completed, requestedUpdates=requested,
                maximumSeriesUpdates=12000)


def prepare(source, folder, plan):
    """完成源终态、逐字状态、历史执行源和采样检查后才创建新运行目录。"""
    source, folder, plan = (Path(p).resolve() for p in (source, folder, plan))
    source.relative_to(ROOT/'outputs')
    folder.relative_to(ROOT/'outputs')
    assert not folder.exists() and plan.is_file()
    receipt = json.loads((source/'engine_process.json').read_text(encoding='utf8'))
    old_manifest = json.loads((source/'study_manifest.json').read_text(encoding='utf8'))
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and receipt['stopReason'] is None
    assert all(receipt['inputsStable'].values())
    assert old_manifest['sourceLocalIteration'] is None and old_manifest['independentInitialization']['accepted']
    assert 'continuation' not in old_manifest, '本阶段只接受既定独立初场系列的唯一续算'
    original_sources = dict(receipt['sourceSHA256'], **old_manifest['sourceSHA256'])
    verified_sources = {}
    for name, expected in original_sources.items():
        path = Path(name)
        if sha(path) != expected:
            relative = path.relative_to(ROOT)
            archived = ROOT/'docs/validation/v151/executed_sources'/relative
            if not archived.is_file() or sha(archived) != expected:
                archived = ROOT/'docs/validation/v152/prior_sources'/relative
            assert archived.is_file() and sha(archived) == expected, name
            path = archived
        verified_sources[str(path)] = expected
    for name, expected in receipt['outputSHA256'].items():
        assert sha(source/name) == expected, name
    history = np.genfromtxt(str(source/'history.csv'), delimiter=',', names=True)
    assert np.array_equal(history['Inner_Iter'], np.arange(len(history)))
    step = int(history['Inner_Iter'][-1])
    assert step == 1199 and len(history) == 1200
    checkpoint = source/('restart_flow_%06d.csv' % step)
    assert sha(checkpoint) == receipt['outputSHA256'][checkpoint.name] == sha(source/'restart_flow.csv')
    assert sha(source/('surface_flow_%06d.vtk' % step)) == sha(source/'surface_flow.vtk')
    a, _, nd, _, _, _ = state(checkpoint, True, True, True)
    points, cells, markers = mesh_data(source/'mesh.su2')
    assert nd == 2 and len(a) == 16261 and len(cells) == 16000
    assert np.array_equal(a[:, 1:3], points)
    contract = ROOT/'tests/data/tmr_mach7/contract.json'
    measurements = json.loads(contract.read_text(encoding='utf8'))
    pairing = []
    for saved in sorted(source.glob('restart_flow_*.csv')):
        index = int(saved.stem.split('_')[-1])
        frame = surface_frame(source/('surface_flow_%06d.vtk' % index), points, markers, measurements)
        errors = compare_volume(frame, saved, points, markers)
        pairing.append(dict(iteration=index, coordinateErrorM=frame['coordinateErrorM'], errors=errors))
    assert len(pairing) == 12
    old_config = config(source/'case.cfg')
    assert old_config['OUTPUT_FILES'] == '(SURFACE_PARAVIEW_ASCII, RESTART_ASCII)'
    changes = dict(ITER='10800', OUTPUT_WRT_FREQ='100, 1000', SCREEN_WRT_FREQ_INNER='500')
    current = dict(old_config)
    current.update(changes)
    assert current['RESTART_SOL'] == 'YES' and current['SOLUTION_FILENAME'] == 'initial_flow'
    assert {k: v for k, v in current.items() if k not in changes} == {
        k: v for k, v in old_config.items() if k not in changes}
    lineage = series_budget(len(history), int(current['ITER']))
    lineage.update(sourceDirectory=str(source), sourceSavedLocalIteration=step,
                   sourceCheckpoint=str(checkpoint), sourceStateSHA256=sha(checkpoint),
                   seriesRootDirectory=str(source), independentSeriesInitialization=old_manifest['independentInitialization'],
                   pseudoTimeHistoryRestarts=True, physicalTimeContinuation=False,
                   method='exact native conservative checkpoint; original global IDs')
    dependencies = [Path(__file__).resolve(), plan, contract,
                    ROOT/'scripts/prepareSU2AxisymmetricInitialState.py',
                    ROOT/'scripts/su2AxisymmetricAnalysis.py', ROOT/'scripts/su2AxisymmetricSurface.py',
                    ROOT/'scripts/auditSU2SSTState.py', ROOT/'scripts/auditSU2AxisymmetricStudy.py',
                    ROOT/'scripts/su2BenchmarkAcceptance.py',
                    source/'engine_process.json', source/'study_manifest.json',
                    source/'case.cfg', source/'mesh.su2', source/'history.csv', checkpoint]
    folder.mkdir(parents=True)
    os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
    os.link(str(checkpoint), str(folder/'initial_flow.csv'))
    (folder/'case.cfg').write_bytes(('\n'.join(k+'= '+v for k, v in current.items())+'\n').encode('ascii'))
    save_json(folder/'surface_pairing.json', dict(passed=True, frames=pairing, bound=6e-7))
    save_json(folder/'continuation.json', lineage)
    dependencies += [folder/'surface_pairing.json', folder/'continuation.json']
    locks = dict(verified_sources)
    locks.update({str(p): sha(p) for p in dependencies})
    manifest = dict(kind=old_manifest['kind'], benchmark=old_manifest['benchmark'],
                    sourceDirectory=str(source), sourceLocalIteration=step, continuation=lineage,
                    sourceFieldSpatialOutput=True, sourceOriginalStartupPassed=True,
                    controlsChanged=changes, samePhysicsAndNumerics=True,
                    maximumSeconds=720, maximumRSSBytes=4*1024**3, maximumOutputBytes=192*1024**2,
                    validateSavedState=True, stabilitySampling='native-surface-100',
                    sourceSHA256=locks, historicalSourceEvidence=verified_sources,
                    plan=str(plan), planSHA256=sha(plan), physicsQualified=False, newQualifiedArea=0)
    save_json(folder/'study_manifest.json', manifest)
    print(json.dumps(dict(prepared=True, sourceLocalIteration=step,
                         completedUpdatesBeforeRun=1200, pairedFrames=len(pairing))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source')
    parser.add_argument('--folder', required=True)
    parser.add_argument('--plan')
    parser.add_argument('--run-prepared', action='store_true')
    args = parser.parse_args()
    if args.run_prepared:
        raise SystemExit(run_prepared(args.folder))
    prepare(args.source, args.folder, args.plan)
