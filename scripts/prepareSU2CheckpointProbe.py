"""从已有真实检查点建立一次原生空间残差诊断；不覆写旧结果，不延长旧任务。"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    source = (ROOT / args.source).resolve()
    dest = (ROOT / args.out).resolve()
    dest.relative_to(ROOT / 'outputs')
    assert not dest.exists()
    receipt = json.loads((source / 'engine_process.json').read_text())
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    for name in ('restart_flow.csv', 'surface_flow.vtk', 'history.csv'):
        assert sha(source / name) == receipt['outputSHA256'][name]
    audit = json.loads((ROOT / 'docs/validation/v109/independent_audit.json').read_text())
    for name in ('restart_flow.csv', 'surface_flow.vtk', 'history.csv'):
        assert sha(source / name) == audit['inputSHA256'][(source / name).relative_to(ROOT).as_posix()]
    config = (source / 'case.cfg').read_text()
    changes = {'ITER': '1', 'OUTPUT_WRT_FREQ': '1', 'SCREEN_WRT_FREQ_INNER': '1',
               'OUTPUT_FILES': '(SURFACE_PARAVIEW_ASCII, PARAVIEW_ASCII, RESTART_ASCII)',
               'VOLUME_OUTPUT': '(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL, MESH_QUALITY)',
               'HISTORY_OUTPUT': '(ITER, RMS_RES, MAX_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)'}
    for key, value in changes.items():
        config, count = re.subn(r'^' + key + r'\s*=.*$', key + '= ' + value, config, flags=re.M)
        assert count == 1, key
    # 8.6.0由VOLUME_OUTPUT中的MESH_QUALITY自动启用，未提供同名WRT配置开关。
    dest.mkdir(parents=True)
    for old, new in [('mesh.su2', 'mesh.su2'), ('restart_flow.csv', 'initial_flow.csv'),
                     ('case_manifest.json', 'source_case_manifest.json')]:
        shutil.copyfile(str(source / old), str(dest / new))
    (dest / 'case.cfg').write_bytes(config.encode('ascii'))
    provenance = dict(sourceDirectory=source.relative_to(ROOT).as_posix(), checkpointIteration=audit['checkpointIteration'],
                      diagnosticIterations=1, changes=changes, samePhysicsAndGeometry=True,
                      residualIdentity='one new native assembly from saved state; not original terminal residual',
                      inputSHA256={p.relative_to(ROOT).as_posix(): sha(p) for p in
                                   [source / n for n in ('case.cfg', 'mesh.su2', 'restart_flow.csv', 'engine_process.json')]
                                   + [Path(__file__).resolve(), ROOT / 'docs/validation/v110/PLAN.md']})
    (dest / 'probe_manifest.json').write_bytes((json.dumps(provenance, indent=2) + '\n').encode('utf8'))
    request = dict(workDirectory=str(dest), configFile='case.cfg', meshFile='mesh.su2', threads=2,
                   maximumSeconds=120, sourceFiles=['scripts/prepareSU2CheckpointProbe.py',
                   (dest / 'probe_manifest.json').relative_to(ROOT).as_posix()])
    (dest / 'engine_request.json').write_bytes((json.dumps(request, indent=2) + '\n').encode('utf8'))
    print(json.dumps({'checkpoint': provenance['checkpointIteration'], 'request': str(dest / 'engine_request.json')}))


if __name__ == '__main__':
    main()
