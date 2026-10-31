"""绑定v105实际MATLAB/SU2终态、准确源码和小型证据；不再启动数值计算。"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'outputs/v105_validation'
TARGET = ROOT/'docs/validation/v105'
FILES = ['buildSU2FlatPlateCase.m','runEmbeddedSU2Case.m','readSU2SurfaceResult.m',
         'tests/studyEmbeddedSU2FlatPlate.m','tests/testEmbeddedSU2Contracts.m',
         'tests/auditEmbeddedSU2FlatPlate.py','scripts/installSU2Runtime.py',
         'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def copy(file, relative):
    dest = TARGET/relative
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        assert sha(dest)==sha(file), str(dest)
    else:
        dest.write_bytes(file.read_bytes())
    return dict(path=dest.relative_to(ROOT).as_posix(),sha256=sha(dest),bytes=dest.stat().st_size)


def main():
    assert not (TARGET/'delivery_manifest.json').exists(), '不覆盖历史交付绑定'
    sources = {}
    for path in FILES + ['scripts/bindEmbeddedSU2Evidence.py']:
        file=ROOT/path
        item=copy(file,Path('source_versions')/(file.name+'.'+sha(file)[:16]+'.raw'))
        sources[sha(file)]=item['path']
    for file in (TARGET/'source_versions').glob('*.raw'):
        sources[sha(file)]=file.relative_to(ROOT).as_posix()

    jobs=[]
    for process in sorted(SOURCE.glob('*_process.json')):
        record=json.loads(process.read_text(encoding='utf-8'))
        assert record['status']=='terminal' and all(record['sourceStableAfter'].values())
        relevant={p:h for p,h in record['sourceSHA256'].items() if p in FILES}
        assert all(h in sources for h in relevant.values()), [p for p,h in relevant.items() if h not in sources]
        name=process.name[:-len('_process.json')]
        jobs.append(dict(name=name,exitCode=record['exitCode'],
                         lockedProjectSourceVersions={p:sources[h] for p,h in relevant.items()},
                         receipt=copy(process,Path('processes')/process.name),
                         log=copy(SOURCE/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))))
    assert len(jobs)==5 and sum(j['exitCode']==0 for j in jobs)==4

    engines=[]
    for relative in ['plate01','plate02','plate03','contracts01/invalid_config','contracts01/time_limit']:
        folder=SOURCE/relative
        record=json.loads((folder/'engine_process.json').read_text(encoding='utf-8'))
        assert record['status']=='terminal' and all(record['inputsStable'].values())
        for path,expected in record['inputSHA256'].items():
            assert sha(Path(path))==expected, path
        for name,expected in record['outputSHA256'].items():
            assert sha(folder/name)==expected, name
        items=[]
        for name in ['case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log']:
            items.append(copy(folder/name,Path('cases')/relative/name))
        if (folder/'surface_flow.vtk').exists():
            items.append(copy(folder/'surface_flow.vtk',Path('cases')/relative/'surface_flow.vtk'))
        engines.append(dict(name=relative,exitCode=record['exitCode'],timedOut=record['timedOut'],
                            elapsedSeconds=record['elapsedSeconds'],files=items))
    copy(SOURCE/'plate03/mesh.su2',Path('inputs')/'plate_mesh.su2')
    copy(SOURCE/'contracts01/report.json',Path('contracts')/'report.json')
    copy(ROOT/'tools/runtime/su2/8.6.0/runtime.json',Path('runtime')/'runtime.json')
    copy(SOURCE/'runtime_help.log',Path('runtime')/'help.log')

    # 从MATLAB保存结果与独立Python输出逐值比较，确保两个读取端指向同一真实壁面。
    data=loadmat(str(SOURCE/'plate03/matlab_receipt.mat'),squeeze_me=True,struct_as_record=False)['data']
    rows=np.genfromtxt(str(TARGET/'plate03_audit_final.csv'),delimiter=',',skip_header=1)
    assert np.array_equal(data.pressureRaw,rows[:,1])
    assert np.array_equal(data.temperatureRaw,rows[:,2])
    assert np.array_equal(data.heatFluxRaw,rows[:,3])
    history=np.genfromtxt(str(SOURCE/'plate03/history.csv'),delimiter=',',names=True)
    hf=history['HF'][-501:]
    integrated=float(np.trapz(rows[:,3],rows[:,0]))
    relative=abs(integrated/float(hf[-1])-1)
    assert relative < 2e-6
    for name in ['02','03']:
        audit=json.loads((TARGET/('plate'+name+'_audit_final.json')).read_text(encoding='utf-8'))
        assert audit['auditSHA256'] in sources
        for file,expected in audit['inputSHA256'].items():
            assert sha(SOURCE/('plate'+name)/file)==expected

    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
              for p in sorted(TARGET.rglob('*')) if p.is_file()]
    result=dict(status='tested-integration-and-limited-numerical-comparison',
                matlabJobs=jobs,engineJobs=engines,sourceVersions=sources,
                finalMatlabAndPythonValuesIdentical=True,plate03SurfaceNodes=65,
                plate03IntegratedHeatPerUnitSpan=integrated,plate03EngineHeatPerUnitSpan=float(hf[-1]),
                integratedHeatRelativeDifference=relative,
                last501HeatHistoryRelativeRange=float(np.ptp(hf)/abs(hf[-1])),
                historyPrecisionCaveat='Zero range at saved decimal precision is not exact physical stationarity.',
                evidence=evidence, researchOutputBytes=sum(p.stat().st_size for p in SOURCE.rglob('*') if p.is_file()),
                runtimeBytes=sum(p.stat().st_size for p in (ROOT/'tools/runtime/su2/8.6.0').rglob('*') if p.is_file()),
                newQualifiedRocketArea=0,productionQualified=False,meshAndDomainIndependenceTested=False,
                publicExperimentCompared=False,mainWorkbenchIntegrated=False)
    (TARGET/'delivery_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['matlabJobs','engineJobs','sourceVersions','evidence']},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
