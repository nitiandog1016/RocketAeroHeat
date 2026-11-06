"""模块用途。保存v164小型证据与准确源码，排除大场、运行包及协作者修改。"""
from pathlib import Path
import json
import shutil

from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

BASE = ROOT/'outputs/v164_validation'
DEST = ROOT/'docs/validation/v164/evidence'


def main():
    assert not DEST.exists()
    assert read(BASE/'review01/report.json')['passed'] and read(BASE/'export_check.json')['passed']
    ui=read(BASE/'matlab/main02_process.json')
    assert ui['exitCode']==0 and all(ui['sourceStableAfter'].values())
    # 模块一。限定本轮目录；小型完整输入/输出保留，真实大数组仅保存其散列。
    selected=[];local_only=[]
    for path in BASE.rglob('*'):
        if not path.is_file():
            continue
        relative=path.relative_to(BASE)
        if relative.parts[0]=='collaborator_sources':
            local_only.append(path);continue
        if path.suffix.lower() not in ('.json','.py','.cpp','.hpp','.m','.md','.log','.cfg','.csv','.vtk','.raw','.txt','.diff','.su2'):
            local_only.append(path);continue
        if path.stat().st_size>2*1024**2:
            local_only.append(path);continue
        selected.append(path)
    # 三次完整源树是相同构建输入，父树已在v163归档；提交差分证明避免重复数MB。
    parent=read(ROOT/'outputs/v163_validation/build03/source_tree.json')
    source_deltas={}
    for label in ('fixture01','fixture02','build01'):
        current=read(BASE/label/'source_tree.json')
        differences={name:dict(parentSHA256=parent.get(name),currentSHA256=value)
                     for name,value in current.items() if parent.get(name)!=value}
        assert len(differences)==2 and not set(parent)-set(current)
        source_deltas[label]=dict(fullLocalTreeSHA256=sha(BASE/label/'source_tree.json'),differences=differences)
    save(BASE/'source_deltas.json',dict(parentTree='docs/validation/v163/evidence/outputs/v163_validation/build03/source_tree.json',
        parentTreeSHA256=sha(ROOT/'outputs/v163_validation/build03/source_tree.json'),deltas=source_deltas))
    selected.append(BASE/'source_deltas.json')
    # 模块二。逐字归档并回读；没有任何删除或对历史对象的原位写入。
    entries=[]
    for source in sorted(set(selected)):
        target=DEST/source.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(source),str(target))
        assert sha(source)==sha(target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=source.relative_to(ROOT).as_posix(),
                            bytes=source.stat().st_size,sha256=sha(source)))
    report=dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries,
        localOnly=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(local_only)],
        localOnlyAreNotGitHubBackups=True,uniqueResultsDeleted=0,cleanupBytes=0,
        originalSmallAllColumnGatePassed=False,physicsQualified=False,newQualifiedArea=0)
    save(DEST.parent/'evidence_manifest.json',report)
    print(json.dumps(dict(files=report['count'],bytes=report['bytes'],localOnlyFiles=len(local_only))))


if __name__=='__main__':
    main()
