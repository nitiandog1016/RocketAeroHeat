"""仅归档小证据和准确执行源，原大场、网格及运行包留在本地。"""
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'outputs/v159_validation'
DEST = ROOT/'docs/validation/v159/evidence'
assert not DEST.exists()
selected = [SOURCE/'native_path_review.json', SOURCE/'run_case.py', Path(__file__).resolve()]
selected += [p for p in (SOURCE/'build01').iterdir() if p.is_file()]
selected.append(SOURCE/'build01/executed_sources/SU2_CFD/src/solvers/CEulerSolver.cpp')
for case in ('small01', 'cancel01', 'public01'):
    folder = SOURCE/case
    selected += [p for p in folder.iterdir() if p.is_file() and p.suffix in ('.json', '.cfg', '.log')]
    selected.append(folder/'history.csv')
    selected += [p for p in (folder/'executed_sources').rglob('*') if p.is_file()]
    if case != 'public01':
        selected += [folder/'mesh.su2', folder/'initial_flow.csv']
        selected += list(folder.glob('restart_flow_*.csv'))+list(folder.glob('surface_flow_*.vtk'))
    else:
        selected += [folder/'surface_flow_000000.vtk', folder/'surface_flow_001199.vtk']
for sub in ('small01/review02', 'public01/public_review02', 'final_review'):
    selected += [p for p in (SOURCE/sub).rglob('*') if p.is_file()]
selected += [p for p in (SOURCE/'public01/ideal_air_audit').iterdir()
             if p.is_file() and (p.suffix in ('.json', '.log', '.py') or p.name == 'wall.csv')]
selected.append(SOURCE/'contracts01/report.json')
selected.append(SOURCE/'contracts01/truncated/inventory.json')
for case in ('planar', 'wrongModel', 'wrongInlet', 'wrongVelocity', 'missingPhase'):
    for name in ('case.cfg', 'case_manifest.json', 'ideal_air_contract.json'):
        selected.append(SOURCE/'contracts01'/case/name)
entries = []
for source in sorted(set(selected)):
    relative = source.relative_to(ROOT)
    target = DEST/relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(source), str(target))
    raw = source.read_bytes()
    assert raw == target.read_bytes()
    entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=relative.as_posix(),
                        bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
manifest = dict(files=entries, bytes=sum(row['bytes'] for row in entries),
    count=len(entries), physicalAllocationMeasured=False, deletedUniqueFiles=0,
    localOnly=['public volume fields and initial state', 'full public mesh', 'native runtime and build objects'],
    inheritedNativeSources='v156 phase-v2 exact sources; v159 changes only CEulerSolver.cpp',
    physicalSolutionQualified=False)
(DEST.parent/'evidence_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(dict(count=manifest['count'], bytes=manifest['bytes']))
