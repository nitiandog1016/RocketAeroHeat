"""仅合并v125两份逐字重复终态别名；不删原场、失败、初态或其他历史结果。"""
from pathlib import Path
import os,json,shutil
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    base=(ROOT/'outputs/v125_validation').resolve();case=base/'phase01';target=base/'storage_receipt.json'
    assert not target.exists() and read(case/'engine_process.json')['status']=='terminal'
    for name in ['sst_phase_view01','sst_phase_view02']:
        r=read(base/(name+'_process.json'));assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
    files=sorted(p for p in base.rglob('*') if p.is_file());before={str(p):sha(p) for p in files}
    def length():return sum({(p.stat().st_dev,p.stat().st_ino):p.stat().st_size for p in files}.values())
    prior=length();free=shutil.disk_usage(ROOT).free;plan=[]
    for stem,suffix in [('restart_flow','.csv'),('surface_flow','.vtk')]:
        source=(case/(stem+'_000000'+suffix)).resolve();alias=(case/(stem+suffix)).resolve()
        source.relative_to(base);alias.relative_to(base)
        assert source!=alias and sha(source)==sha(alias)
        row=dict(source=str(source),destination=str(alias),sha256=sha(alias),bytes=alias.stat().st_size,
                 alreadyLinked=os.path.samefile(source,alias));plan.append(row)
        if row['alreadyLinked']:continue
        pending=alias.with_name(alias.name+'.rah-link-temp');assert not pending.exists()
        os.link(str(source),str(pending));os.replace(str(pending),str(alias));assert os.path.samefile(source,alias)
    for p,h in before.items():assert sha(Path(p))==h
    finalfree=shutil.disk_usage(ROOT).free
    r=dict(checkedFiles=len(files),allOriginalPathsAndSHA256Preserved=True,uniqueResultsDeleted=0,plan=plan,
        distinctFileLengthBefore=prior,distinctFileLengthAfter=length(),
        priorSharedInputFileLength=sum((case/n).stat().st_size for n in ['mesh.su2','initial_flow.csv']),
        replacedCopyBytes=sum(p['bytes'] for p in plan if not p['alreadyLinked']),
        observedDiskFreeIncreaseBytes=finalfree-free,diskFreeBefore=free,diskFreeAfter=finalfree,
        physicalDiskUsageClaimed=False,largeResultsBackedUpToGitHub=False,originalSHA256=before,
        scriptSHA256=sha(Path(__file__).resolve()))
    target.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in r.items() if k not in ['originalSHA256','plan']},indent=2))

if __name__=='__main__':main()
