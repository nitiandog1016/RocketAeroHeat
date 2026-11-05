"""归档本轮小型证据并核对准确来源，大型数值场保留本地。"""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/v149_validation'
DEST = ROOT/'docs/validation/v149'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_bytes((json.dumps(value, ensure_ascii=False, indent=2)+'\n').encode('utf-8'))


def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        assert sha(target) == sha(source)
    else:
        shutil.copyfile(str(source), str(target))


def archive():
    """原生运行只读，归档包含原失败和全部独立数值矩阵。"""
    assert not (DEST/'delivery_manifest.json').exists()
    selected = []
    jacobian = OUT/'jacobian01'
    for p in jacobian.iterdir():
        if p.name not in ('native_jacobian.exe', 'execute.stdout'):
            selected.append(p)
    case = OUT/'acc_sbli02'
    selected += [case/n for n in ('case.cfg', 'benchmark_manifest.json', 'engine_process.json',
                                  'engine.log', 'history.csv', 'audit01.json', 'audit02.json',
                                  'wall_quantities.csv', 'surface_flow.vtk')]
    selected += list(OUT.glob('*request*.json'))
    selected += [OUT/'run_bounded.py']
    processes = {}
    for folder in OUT.iterdir():
        if folder.is_dir() and 'supervisor' in folder.name:
            r = json.loads((folder/'process.json').read_text(encoding='utf-8'))
            assert r['status'] == 'terminal' and r['stopReason'] is None
            assert all(r['sourceStableAfter'].values())
            assert sha(folder/'stdout.log') == r['stdoutSHA256']
            processes[folder.name] = dict(exitCode=r['exitCode'], seconds=r['seconds'],
                                          peakFamilyRSSBytes=r['peakFamilyRSSBytes'], sourceStable=True)
            selected += [folder/'process.json', folder/'stdout.log']
    for p in selected:
        copy(p, DEST/'evidence'/p.relative_to(OUT))

    # 逐字固定官方源码片段，同时保存整文件 SHA 与原行号。
    source_zip = ROOT/'tools/runtime/su2-thermal-build/su2.zip'
    expected_zip = 'c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    assert sha(source_zip) == expected_zip
    excerpts = []
    specs = [
        ('SU2_CFD/src/solvers/CNSSolver.cpp', 'void CNSSolver::BC_Isothermal_Wall_Generic(',
         'void CNSSolver::BC_Isothermal_Wall('),
        ('SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
         'void CFVMFlowSolverBase<V, R>::Friction_Forces(', None),
    ]
    with zipfile.ZipFile(str(source_zip)) as z:
        for suffix, begin, end in specs:
            name = next(n for n in z.namelist() if n.endswith(suffix))
            data = z.read(name)
            text = data.decode('utf-8')
            if 'CFVM' in suffix:
                pos = text.index('HF_Visc_Local += HeatFlux[iMarker][iVertex] * Area;')
                start = text.rfind('      if ((geometry->nodes->GetDomain(iPoint))', 0, pos)
                stop = text.index('    END_SU2_OMP_FOR', pos)
            else:
                start, stop = text.index(begin), text.index(end)
            assert 0 <= start < stop
            fragment = text[start:stop]
            target = DEST/'official_source'/(Path(suffix).name+'.excerpt.txt')
            target.parent.mkdir(exist_ok=True)
            target.write_bytes(fragment.encode('utf-8'))
            excerpts.append(dict(member=name, fileSHA256=hashlib.sha256(data).hexdigest(),
                                 startLine=text[:start].count('\n')+1,
                                 endLine=text[:stop].count('\n')+1,
                                 excerpt=target.name, excerptSHA256=sha(target)))
        licenses = [n for n in z.namelist() if n.rsplit('/', 1)[-1] in ('LICENSE.md', 'COPYING')]
        for name in licenses:
            if name.count('/') == 1:
                (DEST/'official_source'/Path(name).name).write_bytes(z.read(name))
    save(DEST/'official_source/manifest.json', dict(archiveSHA256=expected_zip,
         commit='7048fd35f56d4bf52e586429d6746a36a10255f8', excerpts=excerpts))

    # 审计使用的脚本版本和计划版本必须能按原哈希定位。
    local_candidates = list((DEST/'failures').glob('*')) + list((DEST/'evidence').rglob('*'))
    hash_to_path = {sha(p): str(p.relative_to(ROOT)) for p in local_candidates if p.is_file()}
    versions = []
    for receipt in DEST.glob('evidence/*supervisor*/process.json'):
        r = json.loads(receipt.read_text(encoding='utf-8'))
        for name, expected in r['sourceSHA256'].items():
            path = Path(name)
            if sha(path) == expected:
                retained = str(path.relative_to(ROOT))
            else:
                assert expected in hash_to_path, '准确执行源未保留 '+name
                retained = hash_to_path[expected]
            versions.append(dict(receipt=str(receipt.relative_to(ROOT)), source=name,
                                 sha256=expected, retainedAt=retained))
    owned = [ROOT/'scripts'/n for n in (
        'verifySU2SSTFlowJacobian.py', 'reviewSU2SSTFlowJacobian.py',
        'prepareSU2AxisymmetricBenchmark.py', 'runSU2AxisymmetricBenchmark.py',
        'auditSU2AxisymmetricBenchmark.py')]
    owned += [ROOT/n for n in ('.gitattributes', 'AGENTS.md', 'docs/PROJECT_GOALS.md',
              'docs/SST_JACOBIAN_PUBLIC_BENCHMARK_V149.md',
              'docs/reviews/TECHNICAL_ROUTE_REVIEW_2026-11-01.md')]
    own_data = list((ROOT/'tests/data/tmr_mach7_su2').rglob('*'))
    own_data = [p for p in own_data if p.is_file()]
    output_files = [p for p in OUT.rglob('*') if p.is_file()]
    evidence_files = [p for p in DEST.rglob('*') if p.is_file()]
    size = lambda files: sum(p.stat().st_size for p in files)
    total = size(output_files + evidence_files + owned + own_data)
    assert total < 60*1024**2
    report = dict(evidenceArchivePassed=True, originalStartupPassed=False,
                  thermalReadbackPassed=True, physicsQualified=False, newQualifiedArea=0,
                  previousGoalTurn='progress', processes=processes, sourceVersionRetention=versions,
                  outputsBytes=size(output_files), evidenceBytesBeforeManifest=size(evidence_files),
                  countedBytesBeforeManifest=total, budgetBytes=60*1024**2,
                  byteAccounting='file lengths including archives, excludes previous shared runtime; allocation not measured',
                  uniqueFilesDeleted=0, wholeRocketExecutions=0, publicBenchmarkExecutions=1,
                  sources={str(p.relative_to(ROOT)): sha(p) for p in owned+own_data+evidence_files},
                  localOutputs={str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha(p))
                                for p in output_files})
    save(DEST/'delivery_manifest.json', report)
    print(json.dumps(dict(evidenceArchivePassed=True, countedBytes=total,
                          archivedFiles=len(evidence_files), processes=processes)))


if __name__ == '__main__':
    archive()
