"""模块用途。交付后仅合并本轮逐字初场、网格和终态别名，不删除独有证据。"""
from pathlib import Path
import json
import os
import shutil
import subprocess
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

BASE = ROOT/'outputs/v166_validation'


def main():
    target = BASE/'storage_receipt.json'; assert not target.exists()
    local = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    remote = subprocess.check_output(['git','ls-remote','origin','refs/heads/main'], cwd=ROOT, text=True).split()[0]
    assert local == remote
    assert subprocess.check_output(['git','show','HEAD:docs/FULL_DOMAIN_DEVELOPMENT_V166.md'], cwd=ROOT)
    receipt, ui = read(BASE/'real02/engine_process.json'), read(BASE/'matlab/main01_process.json')
    assert receipt['exitCode'] == 3221225501 and ui['exitCode'] == 0 and receipt['status'] == ui['status'] == 'terminal'
    assert all(receipt['inputsStable'].values()) and all(ui['sourceStableAfter'].values())
    assert read(BASE/'export_check.json')['passed']
    folder = (BASE/'real02').resolve()
    assert receipt['checkpoints']['lastComplete']['localIteration'] == 100
    pairs = [(folder/'restart_flow_000100.csv', folder/'restart_flow.csv'),
             (folder/'surface_flow_000100.vtk', folder/'surface_flow.vtk'),
             (BASE/'real01/initial_flow.csv', folder/'initial_flow.csv')]
    files = sorted(p for p in folder.rglob('*') if p.is_file())
    before = {str(p):sha(p) for p in files}
    for keep, alias in pairs:
        keep=keep.resolve(); alias=alias.resolve(); keep.relative_to((ROOT/'outputs').resolve())
        assert alias.parent == folder and not keep.is_symlink() and not alias.is_symlink()
        assert sha(keep) == sha(alias); before[str(keep)] = sha(keep)
    free_before = shutil.disk_usage(ROOT).free; changed = []
    for keep, alias in pairs:
        row = dict(keptPath=str(keep), replacedAlias=str(alias), bytes=alias.stat().st_size,
                   sha256=sha(alias), alreadyLinked=os.path.samefile(keep,alias)); changed.append(row)
        if row['alreadyLinked']: continue
        temporary = alias.with_name(alias.name+'.rah-v166-link')
        assert temporary.parent == folder and not temporary.exists()
        os.link(str(keep), str(temporary)); os.replace(str(temporary), str(alias))
        assert os.path.samefile(keep,alias)
    assert all(sha(Path(name)) == digest for name,digest in before.items())
    free_after = shutil.disk_usage(ROOT).free
    save(target, dict(passed=True, deliveredGitCommit=local, checkedFiles=len(before), originalSHA256=before,
        pairs=changed, allOriginalPathsAndSHA256Preserved=True, uniqueResultsDeleted=0,
        replacedCopyBytes=sum(r['bytes'] for r in changed if not r['alreadyLinked']),
        observedDiskFreeIncreaseBytes=free_after-free_before, diskFreeBefore=free_before, diskFreeAfter=free_after,
        physicalDiskUsageClaimed=False, largeResultsBackedUpToGitHub=False,
        scriptSHA256=sha(Path(__file__).resolve())))
    result = read(target)
    print(json.dumps({k:result[k] for k in ('passed','checkedFiles','replacedCopyBytes','observedDiskFreeIncreaseBytes')}))


if __name__ == '__main__': main()
