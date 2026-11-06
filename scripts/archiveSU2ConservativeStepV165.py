"""模块用途。归档本轮小型源码/证据，真实大场及未授权并发文件保持本地。"""
import json
import shutil
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save

BASE=ROOT/'outputs/v165_validation'
DEST=ROOT/'docs/validation/v165/evidence'


def main():
    assert not DEST.exists()
    assert read(BASE/'review01/report.json')['passed'] and read(BASE/'export_check.json')['passed']
    ui=read(BASE/'matlab/main01_process.json')
    assert ui['exitCode']==0 and all(ui['sourceStableAfter'].values())
    parent=read(ROOT/'outputs/v164_validation/build01/source_tree.json')
    deltas={}
    for label in ('fixture01','fixture02','build01'):
        current=read(BASE/label/'source_tree.json')
        differences={n:dict(parentSHA256=parent.get(n),currentSHA256=s) for n,s in current.items() if parent.get(n)!=s}
        assert len(differences)==3 and not set(parent)-set(current)
        deltas[label]=dict(fullLocalTreeSHA256=sha(BASE/label/'source_tree.json'),differences=differences)
    save(BASE/'source_deltas.json',dict(parentTree='docs/validation/v164/evidence/outputs/v164_validation/build01/source_tree.json',
        parentTreeSHA256=sha(ROOT/'outputs/v164_validation/build01/source_tree.json'),deltas=deltas))
    selected=[];local=[]
    for path in BASE.rglob('*'):
        if not path.is_file():continue
        if path.relative_to(BASE).parts[0]=='collaborator_sources' or path.suffix.lower() not in ('.json','.py','.cpp','.hpp','.inc','.m','.md','.log','.cfg','.csv','.vtk','.raw','.txt','.diff','.su2') or path.stat().st_size>2*1024**2:
            local.append(path)
        else:selected.append(path)
    entries=[]
    for source in sorted(selected):
        target=DEST/source.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(source),str(target));assert sha(source)==sha(target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=source.relative_to(ROOT).as_posix(),
                            bytes=source.stat().st_size,sha256=sha(source)))
    report=dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries,
        localOnly=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(local)],
        localOnlyAreNotGitHubBackups=True,uniqueResultsDeleted=0,cleanupBytes=0,
        numericalCandidateAccepted=False,physicsQualified=False,newQualifiedArea=0)
    save(DEST.parent/'evidence_manifest.json',report)
    print(json.dumps(dict(files=len(entries),bytes=report['bytes'],localOnlyFiles=len(local))))


if __name__=='__main__':main()
