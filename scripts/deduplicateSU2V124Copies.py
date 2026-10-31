"""仅合并本轮终态中的逐字重复副本，保留被拒场/检查点和全部路径。"""
from pathlib import Path
import json,os,shutil
from prepareSU2CheckpointContinuation import sha,read,ROOT
BASE=ROOT/'outputs/v124_validation'
def within(p,r):
    try:p.relative_to(r);return True
    except ValueError:return False

def main():
    report=BASE/'deduplication.json';assert not report.exists()
    case=Path(read(BASE/'workbench03/engine_folder.json')['folder']).resolve();case.relative_to(ROOT/'outputs/embedded_workbench')
    assert read(case/'workflow/pipeline_process.json')['status']=='terminal'
    assert read(case/'engine_process.json')['status']=='terminal'
    for p in BASE.glob('sst_*_process.json'):
        r=read(p);assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
    files=sorted(p.resolve() for root in [BASE,case] for p in root.rglob('*') if p.is_file())
    before={str(p):(sha(p),p.stat().st_size,p.stat().st_ino) for p in files}
    def distinct():
        ids=set();total=0
        for p in files:
            s=p.stat();key=(s.st_dev,s.st_ino)
            if key not in ids:ids.add(key);total+=s.st_size
        return total
    lengthBefore=distinct();old=ROOT/'outputs/v120_validation/study01/candidate/mesh.su2'
    existing={sha(old):old};plan=[]
    for p in files:
        h,length,_=before[str(p)]
        if length<1024**2:continue
        if h not in existing:existing[h]=p;continue
        source=existing[h].resolve();source.relative_to(ROOT/'outputs');p.relative_to(ROOT/'outputs')
        assert any(within(p,root) for root in [BASE,case]) and source!=p
        plan.append(dict(source=str(source),destination=str(p),sha256=h,bytes=length,alreadyLinked=os.path.samefile(source,p)))
    free=shutil.disk_usage(ROOT).free
    for row in plan:
        if row['alreadyLinked']:continue
        source=Path(row['source']);dest=Path(row['destination']);temp=dest.with_name(dest.name+'.rah-link-temp')
        assert not temp.exists() and sha(source)==sha(dest)==row['sha256']
        os.link(source,temp);os.replace(temp,dest);assert os.path.samefile(source,dest)
    for p,(h,_,_) in before.items():assert sha(Path(p))==h
    after=shutil.disk_usage(ROOT).free
    r=dict(plan=plan,allOriginalPathsAndSHA256Preserved=True,checkedFiles=len(files),uniqueResultsDeleted=0,
        replacedCopyBytes=sum(x['bytes'] for x in plan if not x['alreadyLinked']),diskFreeBefore=free,diskFreeAfter=after,
        observedDiskFreeIncreaseBytes=after-free,distinctFileLengthBefore=lengthBefore,distinctFileLengthAfter=distinct(),
        physicalDiskUsageClaimed=False,largeResultsBackedUpToGitHub=False,originalSHA256={p:v[0] for p,v in before.items()},
        currentDependenciesPreserved=True,scriptSHA256=sha(Path(__file__).resolve()))
    report.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    dest=ROOT/'docs/validation/v124/storage_receipt.json';assert not dest.exists();shutil.copyfile(report,dest)
    print(json.dumps({k:v for k,v in r.items() if k not in ['plan','originalSHA256']},indent=2))

if __name__=='__main__':main()
