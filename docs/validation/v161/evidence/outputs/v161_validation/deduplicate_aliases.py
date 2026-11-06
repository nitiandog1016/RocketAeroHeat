import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"scripts"))
"""合并本轮已经终止的相同输出别名；保留路径、字节和独有场。"""
from pathlib import Path
import json
import os
import shutil
from su2ThermalRuntime import ROOT, read, sha
BASE = ROOT / "outputs/v161_validation"
def write(path, value):
    path.write_text(json.dumps(value, indent=2)+"\n", encoding="utf8")


def main():
    # 模块一：仅允许本轮终态目录，不扫描或更改历史共享存储。
    base = BASE.resolve()
    assert base == (ROOT / 'outputs/v161_validation').resolve()
    target = base / 'storage_receipt.json'
    assert not target.exists()
    ui = read(base / 'matlab/linear_view01_process.json')
    assert ui['status'] == 'terminal' and ui['exitCode'] == 0
    assert all(ui['sourceStableAfter'].values())
    assert read(base / 'export_check.json')['passed']
    files = sorted(p for p in base.rglob('*') if p.is_file())
    before = {str(p): sha(p) for p in files}
    def distinct_length():
        return sum({(p.stat().st_dev, p.stat().st_ino): p.stat().st_size for p in files}.values())
    original_length = distinct_length()
    free_before = shutil.disk_usage(ROOT).free
    pairs = []
    # 模块二：显式限定两组的最终检查点和相同别名，先建立硬链接再原子替换。
    for label in ['control']:
        folder = base / label
        receipt = read(folder / 'engine_process.json')
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 3221225501
        assert all(receipt['inputsStable'].values())
        import psutil
        assert not psutil.pid_exists(receipt['processID']) and not psutil.pid_exists(receipt['supervisorProcessID'])
        for stem, ext in [('restart_flow', '.csv'), ('surface_flow', '.vtk')]:
            source = (folder / (stem + '_000000' + ext)).resolve()
            alias = (folder / (stem + ext)).resolve()
            assert source.parent == folder and alias.parent == folder and source != alias
            assert not source.is_symlink() and not alias.is_symlink()
            assert sha(source) == sha(alias) == receipt['outputSHA256'][alias.name]
            row = dict(source=str(source), alias=str(alias), sha256=sha(alias),
                       bytes=alias.stat().st_size, alreadyLinked=os.path.samefile(source, alias))
            pairs.append(row)
            if row['alreadyLinked']:
                continue
            pending = alias.with_name(alias.name + '.rah-link-temp')
            assert not pending.exists() and pending.parent == folder
            os.link(str(source), str(pending))
            os.replace(str(pending), str(alias))
            assert os.path.samefile(source, alias)
    # 模块三：核对本轮所有原文件的字节；文件长度与磁盘空闲观测分开报告。
    for p, expected in before.items():
        assert sha(Path(p)) == expected, p
    free_after = shutil.disk_usage(ROOT).free
    result = dict(passed=True, checkedFiles=len(before), originalSHA256=before, pairs=pairs,
                  allOriginalPathsAndSHA256Preserved=True, uniqueResultsDeleted=0,
                  replacedCopyBytes=sum(r['bytes'] for r in pairs if not r['alreadyLinked']),
                  distinctFileLengthBefore=original_length, distinctFileLengthAfter=distinct_length(),
                  observedDiskFreeIncreaseBytes=free_after-free_before,
                  diskFreeBefore=free_before, diskFreeAfter=free_after,
                  physicalDiskUsageClaimed=False, largeResultsBackedUpToGitHub=False,
                  scriptSHA256=sha(Path(__file__).resolve()))
    write(target, result)
    print(json.dumps({k: result[k] for k in ['passed', 'checkedFiles', 'replacedCopyBytes', 'observedDiskFreeIncreaseBytes']}))


if __name__ == '__main__':
    main()
