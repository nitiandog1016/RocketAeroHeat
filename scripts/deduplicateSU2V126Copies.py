"""仅合并v126终态逐字副本，保持所有路径和科学字节；不碰历史共享对象。"""
from pathlib import Path
import os,json,shutil
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    base=(ROOT/'outputs/v126_validation').resolve()
    case=Path(read(base/'main01/engine_folder.json')['folder']).resolve()
    case.relative_to(ROOT/'outputs/embedded_workbench')
    target=base/'storage_receipt.json';assert not target.exists()
    for name in ['sst_output_main01','sst_output_main02']:
        r=read(base/(name+'_process.json'));assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
    folders=[base/'format01',base/'restart02',case]
    for folder in folders:assert read(folder/'engine_process.json')['status']=='terminal'
    assert read(case/'workflow/pipeline_process.json')['status']=='terminal'
    files=sorted(p for root in [base,case] for p in root.rglob('*') if p.is_file())
    before={str(p):sha(p) for p in files}
    def length():return sum({(p.stat().st_dev,p.stat().st_ino):p.stat().st_size for p in files}.values())
    prior=length();free=shutil.disk_usage(ROOT).free;pairs=[];plan=[]
    for folder,step in zip(folders,[2,1,49]):
        for stem,ext in [('restart_flow','.csv'),('surface_flow','.vtk')]:
            pairs.append((folder/(stem+'_%06d'%step+ext),folder/(stem+ext)))
    pairs.append((ROOT/'outputs/v125_validation/phase01/mesh.su2',case/'mesh.su2'))
    for source,alias in pairs:
        source=source.resolve();alias=alias.resolve()
        source.relative_to(ROOT/'outputs');alias.relative_to(ROOT/'outputs')
        assert source!=alias and (base in alias.parents or case in alias.parents)
        assert sha(source)==sha(alias)
        row=dict(source=str(source),destination=str(alias),sha256=sha(alias),bytes=alias.stat().st_size,
                 alreadyLinked=os.path.samefile(source,alias));plan.append(row)
        if row['alreadyLinked']:continue
        pending=alias.with_name(alias.name+'.rah-link-temp');assert not pending.exists()
        os.link(str(source),str(pending));os.replace(str(pending),str(alias));assert os.path.samefile(source,alias)
    for p,h in before.items():assert sha(Path(p))==h
    finalfree=shutil.disk_usage(ROOT).free
    r=dict(checkedFiles=len(files),allOriginalPathsAndSHA256Preserved=True,uniqueResultsDeleted=0,plan=plan,
        distinctFileLengthBefore=prior,distinctFileLengthAfter=length(),
        directoryLogicalBytesBefore=sum(Path(p).stat().st_size for p in before),
        priorSharedInputFileLength=(case/'initial_flow.csv').stat().st_size+(case/'mesh.su2').stat().st_size,
        replacedCopyBytes=sum(p['bytes'] for p in plan if not p['alreadyLinked']),
        observedDiskFreeIncreaseBytes=finalfree-free,diskFreeBefore=free,diskFreeAfter=finalfree,
        physicalDiskUsageClaimed=False,largeResultsBackedUpToGitHub=False,originalSHA256=before,
        scriptSHA256=sha(Path(__file__).resolve()))
    target.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in r.items() if k not in ['originalSHA256','plan']},indent=2))


if __name__=='__main__':main()
