"""核验v101实际进程、精确MATLAB/UI版本与小型复现证据，不把全快照当调用闭包。"""
from pathlib import Path
import hashlib
import json
import subprocess
from scipy.io import loadmat
from auditThermalPressureSurfaceFlow import sha

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v101_validation';DEST=ROOT/'docs/validation/v101'
JOBS=['shockFitKernel101a','shockFitKernel101b']+['shockFitExport101'+c for c in 'abcde']+['shockFitLegacy101']


def preserve(path,content):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert path.read_bytes()==content,str(path)
    else:path.write_bytes(content)


def main():
    assert not (DEST/'provenance.json').exists()
    candidates={};versions={};processes=[]
    for path in DEST.rglob('*.raw'):candidates.setdefault(sha(path),path)
    for path in (VAL).glob('view*/workbench/**/*.js'):candidates.setdefault(sha(path),path)
    for path in (VAL).glob('view*/workbench/**/*.css'):candidates.setdefault(sha(path),path)
    for folder in VAL.glob('view*/workbench'):
        candidates.setdefault(sha(folder/'index.html'),folder/'index.html')
        # 几何HTML只在</body>前追加两个脚本；可准确恢复原模板，不推测缺失字节。
        html=(folder/'geometry/index.html').read_bytes()
        append=b'<script src="../ale-data.js"></script><script src="ale-study-geometry.js"></script>'
        assert html.count(append)==1
        raw=html.replace(append,b'');digest=hashlib.sha256(raw).hexdigest()
        if digest not in candidates:
            path=DEST/'raw/UI_versions'/('geometry.html.'+digest[:16]+'.raw');preserve(path,raw);candidates[digest]=path

    def resolve(relative,expected):
        relative=relative.replace('\\','/');key=(relative,expected)
        if key in versions:return
        current=ROOT/relative
        candidate=current if current.exists() and sha(current)==expected else candidates.get(expected)
        if candidate is None:
            old=subprocess.run(['git','show','HEAD:'+relative],cwd=ROOT,capture_output=True)
            assert old.returncode==0 and hashlib.sha256(old.stdout).hexdigest()==expected,relative
            candidate=DEST/'raw/source_versions'/(Path(relative).name+'.'+expected[:16]+'.raw');preserve(candidate,old.stdout)
        elif str(candidate).startswith(str(VAL)):
            archived=DEST/'raw/UI_versions'/(Path(relative).name+'.'+expected[:16]+'.raw');preserve(archived,candidate.read_bytes());candidate=archived
        versions[key]=dict(restorePath=relative,sha256=expected,exactBytesAt=candidate.relative_to(ROOT).as_posix())

    for job in JOBS:
        file=VAL/(job+'_process.json');process=json.loads(file.read_text(encoding='utf-8-sig'))
        code=1 if job=='shockFitKernel101a' else 0
        assert process['status']=='terminal' and process['exitCode']==code
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        preserve(DEST/'raw/processes'/(file.name+'.raw'),file.read_bytes())
        log=VAL/(job+'_stdout.log');preserve(DEST/'raw/logs'/(log.name+'.raw'),log.read_bytes())
        required={p:h for p,h in process['sourceSHA256'].items() if p.endswith('.m')}
        if job.startswith('shockFitExport'):
            required.update({p:h for p,h in process['sourceSHA256'].items() if p in {
                'web/geometry.js','web/geometry.css','web/geometry.html','web/ale-study.js',
                'web/ale-study.css','web/ale-study.html','web/ale-study-geometry.js'}})
        for path,digest in required.items():resolve(path,digest)
        processes.append(dict(job=job,status='terminal',exitCode=code,engineProcessIDs=process['engineProcessIDs'],
            wholeSnapshotStable=True,processSHA256=sha(file),exactMATLABAndUsedUIVersionsVerified=len(required)))
    for path in ['kernel02/shock_fitted_kernel_summary.json','kernel02/independent_audit.json',
                 'ale_ui_node_view05.json','ui_export_independent.json','legacy_node.json']:
        source=VAL/path;preserve(DEST/Path(path).name,source.read_bytes())
    evidence=VAL/'kernel02/shock_fitted_kernel_evidence.mat';s=loadmat(evidence,simplify_cells=True)
    reference=dict(kind='numeric-regression-reference-not-experimental-data',evidenceSHA256=sha(evidence),
        upstream={k:s['up'][k].tolist() if hasattr(s['up'][k],'tolist') else s['up'][k] for k in ['pressure','temperature','density','velocityVector']},
        movingShock=[],piston=[])
    for gas,group,piston in zip(s['gases'],s['jumpCases'],s['pistons']):
        b=group['base'];reference['movingShock'].append(dict(mode=gas['thermalPropertyMode'],normal=b['normal'].tolist(),
            shockNormalSpeedMS=b['shockNormalSpeedMS'].tolist(),pressurePa=b['state']['pressure'].tolist(),
            temperatureK=b['state']['temperature'].tolist(),velocityMS=b['state']['velocity'].tolist()))
        reference['piston'].append(dict(mode=gas['thermalPropertyMode'],relativeShockSpeedMS=float(piston['upstreamRelativeShockSpeedMS']),
            pressurePa=float(piston['jump']['state']['pressure']),temperatureK=float(piston['jump']['state']['temperature']),
            ALEAdvanced=bool(piston['ALEAdvanced'])))
    target=ROOT/'tests/data/shock_fitted_kernel/reference.json';preserve(target,(json.dumps(reference,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    tests=['tests/auditShockFittedALEKernel.py','tests/auditShockFittedExport.py','tests/testShockFittedStudy.js',
           'tests/auditShockFittedProvenance.py','tests/data/shock_layer_transport/parabolic_shock_layer_evidence.mat']
    result=dict(status='passed-terminal-processes-exact-MATLAB-and-used-UI-source-versions',processes=processes,
        sourceVersions=list(versions.values()),currentAuditAndInputSHA256={p:sha(ROOT/p) for p in tests},
        scientificEvidenceSHA256=sha(evidence),allUnrelatedSnapshotBytesClaimedArchived=False,
        productionQualified=False,newProductionSolvedAreaM2=0,goalRemainsActive=True)
    preserve(DEST/'provenance.json',(json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(status='passed',terminalProcesses=len(processes),exactSourceVersions=len(versions),physicalQualification=False)))


if __name__=='__main__':main()
