"""模块用途。归档v166小源码和检查证据，保留大场本地身份及失败记录。"""
import json
import shutil
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

BASE = ROOT/'outputs/v166_validation'
DEST = ROOT/'docs/validation/v166/evidence'


def main():
    assert not DEST.exists()
    assert read(BASE/'review01/report.json')['passedChecks'] and read(BASE/'export_check.json')['passed']
    ui = read(BASE/'matlab/main01_process.json')
    assert ui['exitCode'] == 0 and all(ui['sourceStableAfter'].values())
    assert read(BASE/'real02/engine_process.json')['status'] == 'terminal'
    selected, local = [], []
    suffixes = ('.json','.py','.m','.md','.log','.cfg','.csv','.vtk','.txt','.su2')
    for path in BASE.rglob('*'):
        if not path.is_file(): continue
        if path.relative_to(BASE).parts[0] == 'collaborator_sources' or path.suffix.lower() not in suffixes or path.stat().st_size > 2*1024**2:
            local.append(path)
        else: selected.append(path)
    entries = []
    for source in sorted(selected):
        target = DEST/source.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(source), str(target)); assert sha(source) == sha(target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source.relative_to(ROOT).as_posix(),
                            bytes=source.stat().st_size, sha256=sha(source)))
    report = dict(count=len(entries), bytes=sum(e['bytes'] for e in entries), files=entries,
        localOnly=[dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(local)],
        localOnlyAreNotGitHubBackups=True, uniqueResultsDeleted=0, cleanupBytes=0,
        developmentBenefitPassed=read(BASE/'review01/report.json')['developmentBenefitPassed'],
        physicsQualified=False, newQualifiedArea=0)
    save(DEST.parent/'evidence_manifest.json', report)
    print(json.dumps(dict(files=len(entries), bytes=report['bytes'], localOnlyFiles=len(local))))


if __name__ == '__main__': main()
