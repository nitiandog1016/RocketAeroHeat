"""仅将本轮逐字节相同的完成态副本合并硬链接；独有场、失败和检查点均保留。"""
from pathlib import Path
import hashlib,json,os,shutil
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v122_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=BASE/'deduplication.json';assert not target.exists()
    for p in BASE.rglob('engine_process.json'):assert json.loads(p.read_text())['status']=='terminal'
    pairs=[]
    for folder in ['format_property_equivalence','linked_property_equivalence']:
        for mode in ['PT','DE']:pairs.append((BASE/'actual01/initial_thermodynamics'/(mode+'_output.csv'),BASE/folder/(mode+'_output.csv')))
    for folder in ['guard_control','format_control','linked_control']:
        for name in ['restart_flow.csv','restart_flow_000019.csv','surface_flow.vtk','surface_flow_000019.vtk']:
            pairs.append((ROOT/'outputs/v121_validation/thermal20'/name,BASE/folder/name))
    for name in ['SU2_CFD.exe','RAH_CoolProp_State.exe']:
        pairs.append((ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-recon-v1/bin'/name,ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-recon-v2/bin'/name))
    plan=[]
    for source,dest in pairs:
        source=source.resolve();dest=dest.resolve();source.relative_to(ROOT);dest.relative_to(ROOT)
        assert source!=dest and source.is_file() and dest.is_file() and sha(source)==sha(dest)
        plan.append(dict(source=str(source),destination=str(dest),sha256=sha(source),bytes=dest.stat().st_size,alreadyLinked=os.path.samefile(source,dest)))
    before=shutil.disk_usage(ROOT).free
    for row in plan:
        if row['alreadyLinked']:continue
        source=Path(row['source']);dest=Path(row['destination']);temporary=dest.with_name(dest.name+'.rah-link-temp')
        assert not temporary.exists() and sha(source)==row['sha256'] and sha(dest)==row['sha256']
        os.link(str(source),str(temporary));os.replace(str(temporary),str(dest))
        assert os.path.samefile(source,dest) and sha(dest)==row['sha256']
    result=dict(pairs=plan,allPathsAndSHA256Preserved=True,uniqueResultsDeleted=0,bytesOfReplacedCopies=sum(r['bytes'] for r in plan if not r['alreadyLinked']),
        diskFreeBefore=before,diskFreeAfter=shutil.disk_usage(ROOT).free,scope='only explicitly listed immutable v122 copies; no historical shared objects',largeDataBackedUpToGitHub=False)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='pairs'},indent=2))
if __name__=='__main__':main()
