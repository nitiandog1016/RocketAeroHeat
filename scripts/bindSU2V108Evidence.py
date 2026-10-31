"""绑定v108终态、源版本、小型初场来源和数值核对；不复制工程或运行流场。"""
from pathlib import Path
import hashlib
import gzip
import io
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v108_validation';DEST=ROOT/'docs/validation/v108'
FILES=['scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/prepareSU2FlatPlateWarmStart.py',
       'buildSU2FlatPlateCase.m','buildSU2RocketVolume.m','buildSU2NormalLayerCoordinates.m','writeSU2RocketCase.m',
       'buildRocketTestBody.m','buildSU2WallProjection.m','mapSU2WallResult.m','readSU2SurfaceResult.m',
       'regionalFlow/axisymmetric/buildRocketAxisymmetricDomain.m','regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m',
       'rocketAeroHeatWorkbench.m','createCaseConfig.m','su2WorkbenchMeshOptions.m','aeroHeatWorkflowState.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','web/geometry.js','web/geometry.html','web/geometry.css',
       'tests/auditSU2WarmStart.py','tests/auditSU2WarmConvergence.py','tests/auditSU2NormalLayers.py','tests/auditSU2V108.py',
       'tests/auditEmbeddedSU2FlatPlate.py','tests/studySU2NormalLayers.m','tests/testSU2NormalImport.m','tests/testSU2MeshModeUI.m',
       'scripts/bindSU2V108Evidence.py']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def save_bytes(data,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():assert dst.read_bytes()==data,str(dst)
    else:dst.write_bytes(data)
    return dst.relative_to(ROOT).as_posix()


def copy(p,target):return save_bytes(p.read_bytes(),DEST/target)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    # 前序作业曾加载旧版界面/配置等文件；保存准确字节，不用最终版本冒名。
    jobs=[]
    for receipt in sorted(OUT.glob('*_process.json')):
        state=json.loads(receipt.read_text());assert state['status']=='terminal' and all(state['sourceStableAfter'].values())
        chosen={p:h for p,h in state['sourceSHA256'].items() if p in FILES}
        for path,digest in chosen.items():
            if digest not in versions:
                candidate=subprocess.check_output(['git','show','HEAD:'+path],cwd=str(ROOT))
                if hashlib.sha256(candidate).hexdigest()!=digest:
                    candidates=[p for version in ('v106','v107') for p in (ROOT/'docs/validation'/version/'source_versions').glob('*.raw') if sha(p)==digest]
                    assert candidates,(path,digest);candidate=candidates[0].read_bytes()
                assert hashlib.sha256(candidate).hexdigest()==digest
                versions[digest]=save_bytes(candidate,DEST/'source_versions'/(Path(path).name+'.'+digest[:16]+'.raw'))
        name=receipt.name[:-len('_process.json')]
        copy(receipt,Path('processes')/receipt.name);copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))
        jobs.append(dict(name=name,exitCode=state['exitCode'],sources={p:versions[h] for p,h in chosen.items()}))
    assert len(jobs)==4 and sum(j['exitCode']==0 for j in jobs)==3
    engines=[]
    for sub,target in [('plate_warm01','warm'),('normal01/native_input','normal')]:
        folder=OUT/sub;s=json.loads((folder/'engine_process.json').read_text());assert s['status']=='terminal' and s['exitCode']==0
        assert not s['timedOut'] and not s['cancelled'] and all(s['inputsStable'].values())
        for path,digest in s['inputSHA256'].items():assert sha(Path(path))==digest or digest in versions,(path,digest)
        for name,digest in s['outputSHA256'].items():assert sha(folder/name)==digest,(sub,name)
        for name in ('case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log','surface_flow.vtk'):
            copy(folder/name,Path('engines')/target/name)
        # 大迭代历史留本地；公开原始末两行和完整SHA，不能把摘录当全历史。
        history=folder/'history.csv';lines=history.read_bytes().splitlines(keepends=True)
        save_bytes(b''.join([lines[0]]+lines[-2:]),DEST/'engines'/target/'history_final_excerpt.txt')
        engines.append(dict(folder=folder.relative_to(ROOT).as_posix(),exitCode=s['exitCode'],elapsedSeconds=s['elapsedSeconds'],
                            historySHA256=sha(history),historyRows=len(lines)-1,historyPublication='header-and-final-two-rows-only'))
    for name in ('plate_warm01/initial_audit.json','plate_warm01/warm_start_manifest.json','plate_warm01/physical_audit.json',
                 'plate_warm01/physical_audit.csv','plate_warm01/mesh_comparison.json','plate_warm01/mesh_comparison.csv',
                 'plate_warm01/case_before_warm_start.cfg','normal01/geometry_report.json','normal01/independent_geometry_audit.json',
                 'import01/report.json','import01/mapped_startup.mat','ui02/report.json'):
        copy(OUT/name,Path('checks')/name)
    coarse=ROOT/'outputs/v105_validation/plate03';origin=json.loads((coarse/'engine_process.json').read_text())
    for path,h in origin['inputSHA256'].items():
        p=Path(path)
        if p.name in ('case.cfg','mesh.su2'):assert sha(p)==h
    assert sha(coarse/'restart_flow.csv')==origin['outputSHA256']['restart_flow.csv']
    copy(coarse/'restart_flow.csv',Path('initialization_seed')/'coarse_restart_flow.csv')
    history_bytes=(coarse/'history.csv').read_bytes();buffer=io.BytesIO()
    with gzip.GzipFile(fileobj=buffer,mode='wb',mtime=0) as stream:stream.write(history_bytes)
    compressed=buffer.getvalue()
    assert gzip.decompress(compressed)==history_bytes
    save_bytes(compressed,DEST/'initialization_seed'/'coarse_history.csv.gz')
    shared=[]
    for name,prior in [('mesh.su2','inputs/plate_mesh.su2'),('case.cfg','cases/plate03/case.cfg'),
                       ('case_manifest.json','cases/plate03/case_manifest.json'),('engine_process.json','cases/plate03/engine_process.json'),
                       ('surface_flow.vtk','cases/plate03/surface_flow.vtk'),
                       ('engine_stdout.log','cases/plate03/engine_stdout.log')]:
        p=ROOT/'docs/validation/v105'/prior;assert p.is_file() and sha(p)==sha(coarse/name),(name,prior)
        shared.append(dict(name=name,path=p.relative_to(ROOT).as_posix(),sha256=sha(p)))
    audits=[OUT/'plate_warm01/mesh_comparison.json',DEST/'independent_mapping_audit.json']
    for p in audits:
        audit=json.loads(p.read_text())
        for path,h in audit['inputSHA256'].items():assert sha(ROOT/path)==h,(path,h)
    output_bytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file());assert output_bytes<100*1024*1024
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.rglob('*')) if p.is_file()]
    report=dict(status='tested-internal-warm-start-and-normal-mesh-main-UI; rocket-startup-unqualified',
                matlabJobs=jobs,engineJobs=engines,sourceVersions=versions,sharedCoarseInputs=shared,evidence=evidence,
                outputBytes=output_bytes,evidenceBytes=sum(p['bytes'] for p in evidence),newQualifiedRocketArea=0,
                completeRocketConverged=False,gridIndependenceQualified=False,experimentCompared=False,
                compressedCoarseHistory=dict(path='docs/validation/v108/initialization_seed/coarse_history.csv.gz',
                    uncompressedSHA256=sha(coarse/'history.csv'),uncompressedBytes=len(history_bytes)),
                fullLargeOutputsBackedUpOnGitHub=False)
    (DEST/'delivery_manifest.json').write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:report[k] for k in ('outputBytes','evidenceBytes','newQualifiedRocketArea')}))


if __name__=='__main__':main()
