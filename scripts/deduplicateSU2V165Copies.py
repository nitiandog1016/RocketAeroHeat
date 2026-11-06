"""模块用途。GitHub交付后合并本轮逐字副本，保留所有路径与独有结果。"""
from pathlib import Path
import json
import os
import shutil
import subprocess

from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

BASE=ROOT/'outputs/v165_validation'


def main():
    target=BASE/'storage_receipt.json'
    assert not target.exists()
    local=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    remote=subprocess.check_output(['git','ls-remote',
                                   'origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]
    assert local==remote, '源码交付尚未同步'
    assert subprocess.check_output(['git','show','HEAD:docs/CONSERVATIVE_STEP_V165.md'],cwd=ROOT)
    receipt=read(BASE/'real01/engine_process.json');ui=read(BASE/'matlab/main01_process.json')
    assert receipt['exitCode']==ui['exitCode']==0 and receipt['status']==ui['status']=='terminal'
    assert all(receipt['inputsStable'].values()) and all(ui['sourceStableAfter'].values())
    assert read(BASE/'export_check.json')['passed']
    folder=(BASE/'real01').resolve()
    pairs=[(folder/'restart_flow_000049.csv',folder/'restart_flow.csv'),
           (folder/'surface_flow_000049.vtk',folder/'surface_flow.vtk'),
           (ROOT/'outputs/v164_validation/real01/initial_flow.csv',folder/'initial_flow.csv')]
    files=sorted(p for p in folder.rglob('*') if p.is_file())
    before={str(p):sha(p) for p in files}
    for keep,alias in pairs:
        keep=keep.resolve();alias=alias.resolve()
        keep.relative_to((ROOT/'outputs').resolve())
        assert alias.parent==folder and not keep.is_symlink() and not alias.is_symlink()
        assert sha(keep)==sha(alias)
        before[str(keep)]=sha(keep)
    free_before=shutil.disk_usage(ROOT).free
    changed=[]
    for keep,alias in pairs:
        row=dict(keptPath=str(keep),replacedAlias=str(alias),bytes=alias.stat().st_size,
                 sha256=sha(alias),alreadyLinked=os.path.samefile(keep,alias))
        changed.append(row)
        if row['alreadyLinked']:
            continue
        temporary=alias.with_name(alias.name+'.rah-v165-link')
        assert temporary.parent==folder and not temporary.exists()
        os.link(str(keep),str(temporary))
        os.replace(str(temporary),str(alias))
        assert os.path.samefile(keep,alias)
    for name,digest in before.items():
        assert sha(Path(name))==digest,name
    free_after=shutil.disk_usage(ROOT).free
    result=dict(passed=True,deliveredGitCommit=local,checkedFiles=len(before),originalSHA256=before,
        pairs=changed,allOriginalPathsAndSHA256Preserved=True,uniqueResultsDeleted=0,
        replacedCopyBytes=sum(row['bytes'] for row in changed if not row['alreadyLinked']),
        observedDiskFreeIncreaseBytes=free_after-free_before,diskFreeBefore=free_before,diskFreeAfter=free_after,
        physicalDiskUsageClaimed=False,largeResultsBackedUpToGitHub=False,
        scriptSHA256=sha(Path(__file__).resolve()))
    save(target,result)
    print(json.dumps({k:result[k] for k in ('passed','checkedFiles','replacedCopyBytes','observedDiskFreeIncreaseBytes')}))


if __name__=='__main__':
    main()
