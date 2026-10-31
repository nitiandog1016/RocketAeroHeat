"""绑定v111已结束的数值筛选、原生检查点测试与存储复核；不重跑求解。"""
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v111_validation';DEST=ROOT/'docs/validation/v111'
FILES=['identifySU2NativeCheckpoint.m','auditSU2ResultConvergence.m','writeSU2RocketCase.m',
       'readSU2SurfaceResult.m','mapSU2WallResult.m','buildSU2WallProjection.m',
       'scripts/analyzeSU2DualVolumes.py','scripts/prepareSU2ReconstructionStudy.py',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/diagnoseSU2Checkpoint.py',
       'scripts/mergeSU2V111Duplicates.ps1','scripts/bindSU2V111Evidence.py',
       'tests/auditSU2V111.py','tests/auditSU2V110.py','tests/testSU2NativeCheckpoint.m',
       'tests/testSU2NativeCheckpointRejection.m','docs/validation/v111/PLAN.md']


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))


def save(data,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert path.read_bytes()==data,str(path)
    else:path.write_bytes(data)
    return path.relative_to(ROOT).as_posix()


def copy(path,relative):return save(path.read_bytes(),DEST/relative)


def main():
    assert not (DEST/'delivery_manifest.json').exists(),'Already bound; do not overwrite.'
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    bindings=[]

    def check(path,expected,owner):
        p=Path(path)
        if not p.is_absolute():p=ROOT/p
        if sha(p)==expected:
            bindings.append(dict(owner=owner,path=str(path),sha256=expected,verifiedAgainst='current-file'))
        else:
            assert expected in versions,(str(path),expected)
            bindings.append(dict(owner=owner,path=str(path),sha256=expected,verifiedAgainst=versions[expected]))

    jobs=[]
    for p in sorted(OUT.glob('*_process.json')):
        s=read(p);assert s['status']=='terminal' and s['exitCode']==0 and all(s['sourceStableAfter'].values())
        selected={name:h for name,h in s['sourceSHA256'].items() if name in FILES}
        for name,h in selected.items():
            if h not in versions:
                data=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT))
                assert hashlib.sha256(data).hexdigest()==h,(name,h)
                versions[h]=save(data,DEST/'source_versions'/(Path(name).name+'.'+h[:16]+'.raw'))
        for name in s['explicitInputs']:check(name,s['sourceSHA256'][name],p.name)
        name=p.name[:-len('_process.json')]
        copy(p,Path('processes')/p.name)
        copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))
        jobs.append(dict(name=name,exitCode=s['exitCode'],sourceStable=True,
                         sources={name:versions[h] for name,h in selected.items()}))
    assert len(jobs)==2

    engines=[]
    for name in ('first_probe01','control01','cfl01'):
        folder=OUT/name;s=read(folder/'engine_process.json')
        assert s['status']=='terminal' and s['exitCode']==0 and all(s['inputsStable'].values())
        assert s['request']['threads']==2 and s['elapsedSeconds']<s['request']['maximumSeconds']<=720
        for path,h in s['inputSHA256'].items():check(path,h,name)
        for path,h in s['outputSHA256'].items():check(folder/path,h,name)
        study=read(folder/'study_manifest.json')
        for path,h in study['inputSHA256'].items():check(path,h,name+'/study_manifest')
        for p in sorted(folder.iterdir()):
            if p.suffix in ('.json','.cfg','.log') or p.name in ('history.csv','surface_flow.vtk'):
                copy(p,Path('engine')/name/p.name)
        engines.append(dict(name=name,exitCode=0,elapsedSeconds=s['elapsedSeconds'],
                            threads=2,timedOut=s['timedOut'],physicsQualified=False))
    assert sum(x['elapsedSeconds'] for x in engines)<600
    for name in ('normalization.json','comparison.json'):
        report=read(OUT/name)
        for path,h in report['inputSHA256'].items():check(path,h,name)
        copy(OUT/name,name)
    report=read(OUT/'comparison.json')
    assert report['firstOrderWarmupRejected'] and not report['cflCap50ScreenPassed']
    assert not report['defaultCFLChanged'] and not report['physicsQualified'] and report['newQualifiedArea']==0
    native=read(OUT/'native01/report.json')
    assert [a['checkpointIteration'] for a in native['audits']]==[0,599,599]
    assert all(a['checkpointIdentityConfirmed'] and not a['checkpointResidualPassed'] and not a['productionQualified'] for a in native['audits'])
    negative=read(OUT/'native_reject01/report.json')
    assert all(negative[k] for k in ('missingIndexedReceiptRemainsUnknown','changedIndexedHashRejected','legacyConfigurationKeepsOldPath','originalFilesUnmodified'))
    for name in ('native01','native_reject01'):
        for p in sorted((OUT/name).rglob('*')):
            if p.is_file() and p.suffix in ('.json','.mat','.cfg'):
                copy(p,Path('checks')/name/p.relative_to(OUT/name))

    merge=read(OUT/'duplicate_merge.json');assert len(merge['pairs'])==8
    for pair in merge['pairs']:
        assert os.path.samefile(pair['target'],pair['retained'])
        assert sha(Path(pair['target']))==sha(Path(pair['retained']))==pair['sha256']
    for name in ('duplicate_merge.json','storage_audit.json'):copy(OUT/name,name)
    output_files=[p for p in OUT.rglob('*') if p.is_file()];unique={}
    for p in output_files:
        st=p.stat();unique[(st.st_dev,st.st_ino)]=st.st_size
    assert sum(unique.values())<100*1024**2
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
              for p in sorted(DEST.rglob('*')) if p.is_file()]
    documentation={p:sha(ROOT/p) for p in ('AGENTS.md','docs/PROJECT_GOALS.md','docs/SU2_NORMALIZATION_CHECKPOINT_V111.md')}
    result=dict(status='tested-negative-numerical-screen-and-native-checkpoint-identity; three-dimensional-physics-incomplete',
                matlabJobs=jobs,engineJobs=engines,totalSolverSeconds=sum(x['elapsedSeconds'] for x in engines),
                sourceVersions=versions,inputAndOutputBindings=bindings,evidence=evidence,documentationSHA256=documentation,
                outputPathBytes=sum(p.stat().st_size for p in output_files),outputInodeUniqueBytes=sum(unique.values()),
                evidenceBytes=sum(x['bytes'] for x in evidence),preMergeUniqueBytes=131128192,
                transientUniqueBytesExceeded100MB=True,storageIncludesHistoricalInputLinks=True,
                actualFreeBytesIncrease=merge['freeBytesChange'],allOriginalPathsPreserved=True,
                sourceCheckpointIteration=3750,newQualifiedArea=0,physicsQualified=False,
                mainUINamedCheckpointToggleImplemented=False,mainUITimeoutReopenImplemented=False,
                newBrowserValidationPerformed=False,largeVolumeAndRestartBackedUpOnGitHub=False)
    save((json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode('utf-8'),DEST/'delivery_manifest.json')
    print(json.dumps({k:result[k] for k in ('totalSolverSeconds','outputInodeUniqueBytes','evidenceBytes','newQualifiedArea')}))


if __name__=='__main__':main()
