"""稀疏检查点的容量计划和终态清单，分开记录完整场、失败场及实际更新次数。"""
from pathlib import Path
import csv
import re

import numpy as np

from runEmbeddedSU2 import sha
from su2ThermalRuntime import read

STAGE = 'captured-primitive-k-and-post-turbulence-conservative-state-v1'


def output_plan(nodes, wall_nodes, wall_faces, iterations, every, maximum_bytes, extra_columns=0):
    """模块一。按最长ASCII数值、末次强制场和别名估计文件长度上界。"""
    for value in (nodes, wall_nodes, wall_faces, iterations, every, maximum_bytes):
        assert isinstance(value, int) and not isinstance(value, bool) and value > 0
    assert every <= iterations and wall_nodes <= nodes
    assert type(extra_columns) is int and 0 <= extra_columns <= 16
    # 3D SST含阶段列共27列，每标量预留32字节。表面向量包含三个分量。
    volume = nodes*32*(27+extra_columns)+4096
    surface = wall_nodes*32*(27+extra_columns)+wall_faces*128+8192
    samples = min(iterations, 2+(iterations-1)//every)
    upper = (samples+1)*(volume+surface)+iterations*4096+1024**2
    assert upper <= maximum_bytes, '稀疏输出估计超过本次容量预算'
    result = dict(schemaVersion=1, stage=STAGE, every=every, maximumUpdates=iterations,
        firstLocalIteration=0, firstSuccessfulUpdate=1, maximumNamedPairs=samples,
        estimatedVolumeBytes=volume, estimatedSurfaceBytes=surface,
        estimatedMaximumBytes=upper, maximumOutputBytes=maximum_bytes,
        units='file lengths including final aliases; physical allocation unknown',
        monitor='poll at 0.25 seconds; a write in progress can exceed the stop threshold')
    if extra_columns:
        result.update(schemaVersion=2, extraColumns=extra_columns)
    return result


def output_length(folder, input_paths):
    """仅统计本任务生成文件，既有初场、网格和输入清单单独保留。"""
    return sum(path.stat().st_size for path in folder.iterdir()
               if path.is_file() and str(path) not in input_paths and path.name != 'engine_process.json')


def inventory(folder, receipt):
    """模块二。进程终止后审核每对具名文件，缺帧或截断文件保留拒绝原因。"""
    from auditSU2IdealAirState import load_state
    from auditSU2SSTState import vtk_tokens
    from su2SSTInterface import mesh

    manifest = read(folder/'case_manifest.json')
    if manifest.get('caseKind') == 'naca-tn3344-spatial-step-cold-restart-v1':
        from su2NACASpatialStepInterface import inventory as spatial_step_inventory
        return spatial_step_inventory(folder, receipt)
    if manifest.get('caseKind') == 'naca-tn3344-adaptive-cold-restart-v1':
        from su2NACAAdaptiveInterface import inventory as adaptive_naca_inventory
        return adaptive_naca_inventory(folder, receipt)
    if manifest.get('caseKind') == 'naca-tn3344-spatial-wall-laminar-v1':
        from su2NACASpatialInterface import inventory as spatial_naca_inventory
        return spatial_naca_inventory(folder, receipt)
    if manifest.get('caseKind') == 'naca-tn3344-quarter-3d-laminar-v1':
        # 三维层流没有SST阶段；多壁温marker由专用清单审核，不冒用旧wall字段。
        from su2NACALaminarInterface import inventory as naca_inventory
        return naca_inventory(folder, receipt)
    contract = manifest['idealAir']['checkpointContract']
    points, _, wall = mesh(folder/'mesh.su2')
    wall_ids = np.unique([index for face in wall for index in face[1:]])
    names = receipt['outputSHA256']
    matches = [re.fullmatch(r'(?:restart_flow|surface_flow)_(\d+)\.(?:csv|vtk)', name) for name in names]
    steps = sorted({int(match[1]) for match in matches if match is not None})
    history = {}
    history_path = folder/'history.csv'
    if history_path.exists():
        with history_path.open(encoding='utf8') as stream:
            rows = csv.reader(stream, skipinitialspace=True)
            keys = [key.strip() for key in next(rows, [])]
            for row in rows:
                try:
                    values = np.array(row, float)
                    if len(values) == len(keys) and np.isfinite(values).all():
                        index = values[keys.index('Inner_Iter')]
                        if index == int(index):
                            history[int(index)] = dict(zip(keys, values.tolist()))
                except (ValueError, OverflowError):
                    continue
    complete, rejected = [], []
    offset = manifest['initialization']['priorSuccessfulUpdates']
    for step in steps:
        paths = [folder/('restart_flow_%06d.csv' % step), folder/('surface_flow_%06d.vtk' % step)]
        try:
            assert step in history and 0 <= step < contract['maximumUpdates'], '缺少同索引有限历史行'
            assert all(index in history for index in range(step+1)), '检查点前历史更新不连续'
            for path in paths:
                assert path.is_file() and sha(path) == names[path.name], '文件缺失或身份变化'
            spatial = manifest['requiredRuntimeFlavor'] == 'rah-ideal-air-axisymmetric-v1'
            data, fields, nd = load_state(paths[0], True, captured=True, spatial=spatial)
            assert data.shape[0] == len(points) and nd == points.shape[1]
            assert np.allclose(data[:, 1:nd+1], points, rtol=2e-14, atol=1e-14)
            xyz, surface = vtk_tokens(paths[1])
            required = {'Density', 'Momentum', 'Energy', 'Turb_Kin_Energy', 'Omega', 'Pressure',
                'Temperature', 'Mach', 'Pressure_Coefficient', 'Velocity', 'Laminar_Viscosity',
                'Skin_Friction_Coefficient', 'Heat_Flux', 'Y_Plus', 'Eddy_Viscosity', 'Flow_Primitive_K'}
            if spatial:
                required.update(('Residual_Density', 'Residual_Momentum', 'Residual_Energy', 'Residual_TKE', 'Residual_Omega'))
            assert set(surface) == required, '表面字段缺失或未经声明'
            assert len(xyz) == len(wall_ids) and np.isfinite(xyz).all()
            assert all(np.isfinite(value).all() for value in surface.values())
            complete.append(dict(localIteration=step, successfulUpdates=step+1,
                seriesSuccessfulUpdates=offset+step+1, volume=paths[0].name, surface=paths[1].name,
                volumeSHA256=sha(paths[0]), surfaceSHA256=sha(paths[1]),
                stage=STAGE, physicalQualified=False))
        except (AssertionError, ValueError, KeyError, StopIteration, RuntimeError, OSError) as error:
            rejected.append(dict(localIteration=step, reason=type(error).__name__+': '+str(error)))
    last = complete[-1] if complete else None
    last_history = max(history) if history else None
    normal = receipt['exitCode'] == 0 and not any(receipt.get(key, False)
        for key in ('timedOut', 'cancelled', 'outputBudgetExceeded', 'supervisorError'))
    expected = sorted(set(range(0, last_history+1, contract['every'])) | {last_history}) if last_history is not None else []
    missing = sorted(set(expected)-{row['localIteration'] for row in complete})
    return dict(stage=STAGE, complete=complete, rejected=rejected, lastComplete=last,
        lastFiniteHistoryLocalIteration=last_history,
        historyRows=len(history), historyContiguous=sorted(history) == list(range(len(history))),
        seriesUpdatesAtLastFiniteHistory=offset+last_history+1 if last_history is not None else offset,
        normalTerminal=normal, normalFinalPairPresent=bool(normal and last and last['localIteration'] == last_history),
        expectedLocalIterations=expected, missingScheduledLocalIterations=missing,
        normalScheduleComplete=bool(normal and last and not missing and not rejected),
        abnormalInstantStateAvailable=False, physicsQualified=False)
