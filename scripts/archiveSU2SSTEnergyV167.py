"""模块用途。保存v167准确构建/原生/界面小证据，大场和运行库保持本地。"""
import json
import shutil
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save

BASE=ROOT/'outputs/v167_validation';DEST=ROOT/'docs/validation/v167/evidence'


def main():
    assert not DEST.exists()
    for n in ('review01/report.json','export_check.json','common_frame_comparison.json'):
        r=read(BASE/n);assert r.get('passed',r.get('passedChecks'))
    ui=read(BASE/'matlab/main01_process.json');assert ui['exitCode']==0 and all(ui['sourceStableAfter'].values())
    assert read(BASE/'real01/engine_process.json')['status']=='terminal'
    entries=[];local=[]
    suffixes={'.json','.py','.m','.md','.log','.cfg','.csv','.vtk','.txt','.su2','.cpp','.hpp','.diff'}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        if p.relative_to(BASE).parts[0]=='collaborator_sources' or p.suffix.lower() not in suffixes or p.stat().st_size>2*1024**2:
            local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)));continue
        target=DEST/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(p),str(target));assert sha(p)==sha(target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),
                            bytes=p.stat().st_size,sha256=sha(p)))
    save(DEST.parent/'evidence_manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries,
        localOnly=local,localOnlyAreNotGitHubBackups=True,uniqueResultsDeleted=0,
        developmentBenefitPassed=False,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),localOnlyFiles=len(local))))


if __name__=='__main__':main()
