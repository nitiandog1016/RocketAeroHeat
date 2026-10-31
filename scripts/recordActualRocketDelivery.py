"""封存v118真实完整域准备、预算取消/重试和原面科研结果；大体场只登记不复制。"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/v118_validation'
DEST = ROOT / 'docs/validation/v118'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        assert sha(source) == sha(target), str(target)
    else:
        shutil.copyfile(str(source), str(target))


def archive(path):
    digest = sha(path)
    target = DEST / 'source_versions' / (path.name + '.' + digest[:16] + '.raw')
    copy(path, target)
    return dict(sha256=digest, archived=target.relative_to(ROOT).as_posix())


def main():
    target = DEST / 'delivery_manifest.json'
    assert not target.exists()
    common = ['scripts/runMatlabResearch.py', 'buildSU2HybridVolume.m', 'readSU2HybridVolume.m',
              'writeSU2RocketCase.m', 'writeSU2HybridInitialState.m', 'buildSU2WallProjection.m',
              'resolveSU2ThermalBoundary.m', 'scripts/meshSU2HybridVolume.py']
    view = ['createSU2ResearchView.m', 'readSU2SurfaceResult.m', 'mapSU2WallResult.m',
            'attachSU2ThermalBoundary.m', 'auditSU2ResultConvergence.m', 'identifySU2NativeCheckpoint.m',
            'describeSU2Runtime.m', 'describeSU2Convergence.m', 'createAerothermalGeometryView.m',
            'updateAerothermalGeometryView.m', 'web/geometry-su2.html', 'web/geometry.html',
            'web/geometry.css', 'web/geometry.js', 'web/ale-section-core.js', 'web/su2-section-view.js']
    phases = {
        'prepare01': common + ['tests/prepareActualRocketFullDomain.m'],
        'controls01': common + ['tests/testSU2InitializationWorkbench.m', 'rocketAeroHeatWorkbench.m',
                               'estimateSU2CheckpointStorage.m', 'aeroHeatWorkflowState.m'],
        'review01': common + view + ['tests/reviewActualRocketFullDomain.m']}
    sources = {}
    jobs = {}
    for name, used in phases.items():
        receipt = read(OUT / (name + '_process.json'))
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
        assert all(receipt['sourceStableAfter'].values())
        for rel in used:
            assert sha(ROOT / rel) == receipt['sourceSHA256'][rel], (name, rel)
            sources[rel] = archive(ROOT / rel)
        jobs[name] = dict(exitCode=0, sourceStable=True, sourceBindings={rel:sources[rel]['sha256'] for rel in used})
        for suffix in ['_process.json', '.log', '_stdout.log']:
            copy(OUT / (name+suffix), DEST / 'processes' / (name+suffix))
    for rel in ['scripts/prepareSU2BudgetedRetry.py', 'scripts/runEmbeddedSU2.py', 'scripts/su2ThermalRuntime.py',
                'tests/auditActualRocketPreparation.py', 'tests/auditActualRocketResult.py', 'tests/auditSU2V110.py',
                'tests/summarizeActualRocketProgress.py', 'scripts/recordActualRocketDelivery.py']:
        sources[rel] = archive(ROOT / rel)
    mesh = OUT / 'actual01/mesh'
    rec = read(mesh / 'mesh_process.json')
    assert rec['status'] == 'terminal' and rec['exitCode'] == 0 and rec['inputStable'] and rec['sourceStable']
    assert sha(Path(rec['input'])) == rec['inputSHA256']
    assert sha(ROOT / 'scripts/meshSU2HybridVolume.py') == rec['scriptSHA256']
    for rel, digest in rec['outputSHA256'].items():
        assert sha(mesh / rel) == digest
    for name in ['mesh_process.json', 'mesh_report.json', 'gmsh.log']:
        copy(mesh / name, DEST / 'mesh' / name)
    for name in ['mesh_input.mat', 'surface.mat', 'preparation.json']:
        copy(OUT / 'actual01' / name, DEST / 'inputs' / name)
    engines = {}
    for name in ['engine', 'budget_run']:
        folder = OUT / 'actual01' / name
        rec = read(folder / 'engine_process.json')
        assert rec['status'] == 'terminal' and all(rec['inputsStable'].values())
        assert not rec['timedOut']
        assert rec['cancelled'] if name == 'engine' else rec['exitCode'] == 0 and not rec['cancelled']
        for rel, digest in rec['inputSHA256'].items():
            assert sha(Path(rel)) == digest, rel
        for rel, digest in rec['outputSHA256'].items():
            assert sha(folder / rel) == digest, rel
        for file in ['case.cfg', 'case_manifest.json', 'initialization_manifest.json', 'engine_request.json',
                     'engine_process.json', 'engine_stdout.log', 'history.csv']:
            copy(folder / file, DEST / 'cases' / name / file)
        if name == 'engine':
            assert not list(folder.glob('surface_flow*')) and not list(folder.glob('restart_flow*'))
            copy(folder / 'engine_cancel.json', DEST / 'cases' / name / 'engine_cancel.json')
        else:
            for step in [150, 299]:
                file = 'surface_flow_%06d.vtk' % step
                copy(folder / file, DEST / 'cases' / name / file)
        engines[name] = dict(exitCode=rec['exitCode'], cancelled=rec['cancelled'], timedOut=False,
                             seconds=rec['elapsedSeconds'], progress=rec['progress'], physicsQualified=False)
    for name in ['preflight_audit.json', 'budget_revision.json', 'independent_audit.json', 'checkpoint_progress.json']:
        copy(OUT / name, DEST / 'checks' / name)
    for name in ['controls01', 'review01']:
        for file in (OUT / name).iterdir():
            if file.is_file() and file.suffix in ['.mat', '.csv', '.json']:
                copy(file, DEST / 'checks' / name / file.name)
    # 非本轮工作及正式运行时保持。仅读取旧清单，不覆盖其原始身份。
    preserved = read(ROOT / 'outputs/v117_validation/preflight_preservation.json')
    for rel, digest in preserved.items():
        assert sha(ROOT / rel) == digest, rel
    preservation = dict(reference='outputs/v117_validation/preflight_preservation.json',
                        checked=len(preserved), allUnchanged=True, sha256=sha(ROOT / 'outputs/v117_validation/preflight_preservation.json'))
    (DEST / 'checks/preservation.json').write_text(json.dumps(preservation, indent=2)+'\n', encoding='utf-8')
    # st_dev/st_ino按身份去重，统计的是文件长度，绝不称NTFS实际占用或可回收空间。
    files = [p for p in OUT.rglob('*') if p.is_file()]
    identities = {}
    inventory = []
    for p in files:
        st = p.stat()
        key = (st.st_dev, st.st_ino)
        assert st.st_ino != 0
        identities.setdefault(key, st.st_size)
        inventory.append(dict(path=p.relative_to(ROOT).as_posix(), bytes=st.st_size,
                              fileID=str(st.st_ino), sha256=sha(p)))
    storage = dict(files=len(files), uniqueFileObjects=len(identities),
                   pathLengthBytes=sum(p.stat().st_size for p in files), uniqueFileLengthBytes=sum(identities.values()),
                   goalBytes=100*1024**2, exceeded=sum(identities.values()) > 100*1024**2,
                   actualDiskAllocationMeasured=False, uniqueScientificResultsDeleted=False,
                   v117BlockedCleanupRetried=False, inventory=inventory)
    (DEST / 'checks/local_storage.json').write_text(json.dumps(storage, indent=2)+'\n', encoding='utf-8')
    result = dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=str(ROOT)).decode().strip(),
                  status='tested-original-domain-numerical-startup; not-physically-qualified',
                  physicsQualified=False, newQualifiedArea=0, sources=sources, matlabJobs=jobs, engineJobs=engines,
                  evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
                  localOnly=['full hybrid volume and SU2 mesh', 'initial/restart conservative fields',
                             'official/controlled runtimes and build caches', 'original historical baselines'])
    target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(result['evidence']), sources=len(sources), matlabJobs=len(jobs),
                         engineJobs=len(engines), uniqueLocalFileLengthBytes=storage['uniqueFileLengthBytes'])))


if __name__ == '__main__':
    main()
