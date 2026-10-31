"""绑定v112真实数值/UI进程、缺失终态、修复证据与存储复核，不执行求解。"""
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v112_validation';DEST=ROOT/'docs/validation/v112'
FILES=['rocketAeroHeatWorkbench.m','aeroHeatWorkflowState.m','embeddedSU2Job.m','estimateSU2CheckpointStorage.m',
       'readSU2TimeoutCheckpoint.m','identifySU2NativeCheckpoint.m','auditSU2ResultConvergence.m',
       'writeSU2RocketCase.m','readSU2SurfaceResult.m','mapSU2WallResult.m','buildSU2WallProjection.m',
       'buildSU2RocketVolume.m','buildSU2NormalLayerCoordinates.m','su2WorkbenchMeshOptions.m','buildRocketTestBody.m',
       'createCaseConfig.m','createSU2ResearchView.m','describeSU2Convergence.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m',
       'web/geometry.js','web/geometry.css','web/geometry-su2.html','web/ale-section-core.js','web/su2-section-view.js',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/prepareSU2MultigridStudy.py',
       'scripts/mergeSU2V112Duplicates.ps1','scripts/bindSU2V112Evidence.py',
       'tests/testSU2CheckpointWorkbench.m','tests/testSU2ReceiptContention.py','tests/reviewSU2V112.m',
       'tests/testSU2ResearchReplay.m','tests/auditSU2V112.py','tests/auditSU2V110.py','docs/validation/v112/PLAN.md']


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def save(data,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==data,str(p)
    else:p.write_bytes(data)
    return p.relative_to(ROOT).as_posix()


def copy(p,dest):return save(p.read_bytes(),DEST/dest)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    bindings=[]
    def check(path,h,owner):
        p=Path(path)
        if not p.is_absolute():p=ROOT/p
        matches=sha(p)==h
        assert matches or h in versions,(str(p),h)
        bindings.append(dict(owner=owner,path=str(path),sha256=h,verifiedAgainst='current-file' if matches else versions[h]))
    jobs=[]
    for p in sorted(OUT.glob('*_process.json')):
        s=read(p);assert s['status']=='terminal' and all(s['sourceStableAfter'].values())
        selected={name:h for name,h in s['sourceSHA256'].items() if name in FILES}
        for name,h in selected.items():
            if h not in versions:
                b=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT));assert hashlib.sha256(b).hexdigest()==h
                versions[h]=save(b,DEST/'source_versions'/(Path(name).name+'.'+h[:16]+'.raw'))
        for name in s['explicitInputs']:check(name,s['sourceSHA256'][name],p.name)
        name=p.name[:-len('_process.json')];copy(p,Path('processes')/p.name)
        copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))
        jobs.append(dict(name=name,exitCode=s['exitCode'],sourceStable=True,sources={n:versions[h] for n,h in selected.items()}))
    assert len(jobs)==3 and sum(j['exitCode']==0 for j in jobs)==2
    ui=read(OUT/'ui02/report.json');failed=ROOT/'outputs/embedded_workbench/d7a55815-b446-4345-b076-cbac806f8b3e'
    folders=[('mg01',OUT/'mg01'),('ui_failed',failed),('ui_timeout',Path(ui['actualTimeoutFolder']))]
    engines=[]
    for name,folder in folders:
        s=read(folder/'engine_process.json')
        for path,h in s['inputSHA256'].items():check(path,h,name)
        if name!='ui_failed':
            assert s['status']=='terminal' and all(s['inputsStable'].values())
            for path,h in s['outputSHA256'].items():check(folder/path,h,name)
            engines.append(dict(name=name,status='terminal',exitCode=s['exitCode'],timedOut=s['timedOut'],
                                elapsedSeconds=s['elapsedSeconds'],threads=s['request']['threads']))
        else:
            observation=read(OUT/'failed_supervisor_observation.json')
            assert s['status']=='running' and 'exitCode' not in s and not observation['childPresent'] and not observation['supervisorPresent']
            assert 'PermissionError' in (folder/'supervisor.log').read_text(encoding='utf-8',errors='replace')
            engines.append(dict(name=name,status='original-terminal-receipt-missing; independently-observed-stopped',
                                originalExitCodeKnown=False,postHocInputSHAValidated=True,
                                observedFiles={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))
        if (folder/'study_manifest.json').exists():
            for path,h in read(folder/'study_manifest.json')['inputSHA256'].items():check(path,h,name+'/study')
        for p in sorted(folder.iterdir()):
            if p.is_file() and (p.suffix in ('.json','.cfg','.log') or p.name in ('history.csv','surface_flow.vtk','input_snapshot.mat')):
                copy(p,Path('engine')/name/p.name)
    for name in ('independent_audit.json','failed_supervisor_observation.json','storage_plan.json','storage_receipt.json','storage_audit.json'):
        copy(OUT/name,name)
    audit=read(OUT/'independent_audit.json')
    for path,h in audit['inputSHA256'].items():check(path,h,'independent-audit')
    assert not audit['multigridScreenPassed'] and audit['uiCheckpointIteration']==100 and audit['uiTerminalIteration']==103
    for name in ('ui01','ui02','review01','receipt_lock01','receipt_lock02'):
        for p in sorted((OUT/name).rglob('*')):
            if p.is_file():copy(p,Path('checks')/name/p.relative_to(OUT/name))
    for path,h in read(OUT/'receipt_lock02/report.json')['sourceSHA256'].items():check(path,h,'actual-Windows-sharing-lock')
    merge=read(OUT/'storage_receipt.json');assert len(merge['pairs'])==11
    for pair in merge['pairs']:
        assert os.path.samefile(pair['target'],pair['retained'])
        assert sha(Path(pair['target']))==sha(Path(pair['retained']))==pair['sha256']
    allfiles=[p for p in OUT.rglob('*') if p.is_file()]
    for _,folder in folders[1:]:allfiles += [p for p in folder.rglob('*') if p.is_file()]
    unique={(p.stat().st_dev,p.stat().st_ino):p.stat().st_size for p in allfiles}
    assert sum(unique.values())<100*1024**2
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.rglob('*')) if p.is_file()]
    docs={p:sha(ROOT/p) for p in ('AGENTS.md','docs/PROJECT_GOALS.md','docs/SU2_MULTIGRID_CHECKPOINT_UI_V112.md')}
    manifest=dict(status='tested-main-timeout-checkpoint-and-receipt-contention-fix; multigrid-screen-failed',
                  baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
                  matlabJobs=jobs,engineJobs=engines,sourceVersions=versions,inputAndOutputBindings=bindings,
                  evidence=evidence,documentationSHA256=docs,evidenceBytes=sum(x['bytes'] for x in evidence),
                  outputPathBytes=sum(p.stat().st_size for p in allfiles),outputInodeUniqueBytes=sum(unique.values()),
                  includesHistoricalInputLinks=True,storageFreeIncrease=merge['observedFreeIncrease'],
                  newQualifiedArea=0,physicsQualified=False,mainUIActualTimeoutValidated=True,
                  UIScreenshotCondition='M4 alpha20 beta5; actual timeout saved100 monitored103; not M5 multigrid',
                  largeRestartAndVolumeBackedUpOnGitHub=False)
    save((json.dumps(manifest,ensure_ascii=False,indent=2)+'\n').encode('utf-8'),DEST/'delivery_manifest.json')
    print(json.dumps({k:manifest[k] for k in ('outputInodeUniqueBytes','evidenceBytes','storageFreeIncrease','newQualifiedArea')}))


if __name__=='__main__':main()
