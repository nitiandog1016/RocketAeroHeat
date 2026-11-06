"""仅合并本轮两份已终止任务的逐字初场副本，保留原路径和全部独有文件。"""
import os
from pathlib import Path
import shutil
import sys

import psutil

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import atomic_write_json, sha
from su2ThermalRuntime import read

BASE = ROOT/'outputs/v162_validation'
receipt_path = BASE/'storage_receipt.json'
assert not receipt_path.exists()
source = (ROOT/'outputs/v161_validation/control/initial_flow.csv').resolve()
assert source.is_file() and not source.is_symlink()
aliases = [(BASE/label/'initial_flow.csv').resolve() for label in ('debug01', 'debug02')]
for alias in aliases:
    alias.relative_to(BASE.resolve())
    run = read(alias.parent/'debugger_process.json')
    assert run['status'] == 'terminal' and all(run['inputsStable'].values())
    assert not any(psutil.pid_exists(pid) for pid in run['observedProcessIDs'])
    assert sha(alias) == sha(source) == run['inputSHA256'][str(alias)]
files = [p for p in BASE.rglob('*') if p.is_file()]
before = {str(p):sha(p) for p in files}
free_before = shutil.disk_usage(ROOT).free
pairs = []
for alias in aliases:
    assert not alias.is_symlink() and not os.path.samefile(source, alias)
    pending = alias.with_name(alias.name+'.link-pending')
    assert not pending.exists() and pending.parent == alias.parent
    os.link(str(source), str(pending))
    os.replace(str(pending), str(alias))
    assert os.path.samefile(source, alias)
    pairs.append(dict(source=str(source), alias=str(alias), bytes=alias.stat().st_size, sha256=sha(alias)))
assert all(sha(Path(name)) == digest for name, digest in before.items())
atomic_write_json(receipt_path, dict(passed=True, pairs=pairs, originalSHA256=before,
    checkedFiles=len(before), allOriginalPathsAndSHA256Preserved=True,
    replacedCopyBytes=sum(p['bytes'] for p in pairs),
    observedDiskFreeIncreaseBytes=shutil.disk_usage(ROOT).free-free_before,
    uniqueResultsDeleted=0, physicalDiskUsageClaimed=False))
print('Deduplicated', sum(p['bytes'] for p in pairs), 'bytes across', len(before), 'unchanged files')
