"""仅合并本轮已终态、逐字相同的最终检查点副本；保留全部原路径与独有状态。"""
from pathlib import Path
import hashlib,json,os,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v120_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'storage_receipt.json';assert not target.exists();pairs=[]
    for name in ['control','candidate']:
        folder=OUT/'study01'/name;r=json.loads((folder/'engine_process.json').read_text());assert r['status']=='terminal' and r['exitCode']==0
        for stem,suffix in [('restart_flow','csv'),('surface_flow','vtk')]:
            a=(folder/(stem+'.'+suffix)).resolve();b=(folder/(stem+'_000199.'+suffix)).resolve()
            a.relative_to(OUT);b.relative_to(OUT);h=sha(a);assert h==sha(b)==r['outputSHA256'][a.name]==r['outputSHA256'][b.name]
            pairs.append(dict(target=str(a),retained=str(b),sha256=h,bytes=a.stat().st_size,oldFileID=str(a.stat().st_ino)))
    freeBefore=shutil.disk_usage(str(ROOT)).free
    for pair in pairs:
        a=Path(pair['target']);b=Path(pair['retained']);assert sha(a)==sha(b)==pair['sha256']
        a.unlink();os.link(str(b),str(a))
        assert sha(a)==pair['sha256'] and a.stat().st_ino==b.stat().st_ino
        pair['newFileID']=str(a.stat().st_ino)
    freeAfter=shutil.disk_usage(str(ROOT)).free
    result=dict(pairs=pairs,replacedDuplicateFileBytes=sum(p['bytes'] for p in pairs),diskFreeBefore=freeBefore,diskFreeAfter=freeAfter,
                observedDiskFreeIncrease=freeAfter-freeBefore,uniqueStatesDeleted=False,originalPathsPreserved=True,
                note='window free-space change includes unrelated activity; file length is not exact NTFS allocation')
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='pairs'}))
if __name__=='__main__':main()
