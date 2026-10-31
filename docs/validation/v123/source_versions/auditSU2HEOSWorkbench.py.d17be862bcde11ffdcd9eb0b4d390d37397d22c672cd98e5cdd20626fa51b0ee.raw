"""独立按原VTK核对主工作台CSV、真实分阶段终态和原生拾取值。"""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    base=ROOT/'outputs/v123_validation';target=base/'independent_audit.json';assert not target.exists()
    engine=Path(read(base/'ui01/engine_folder.json')['folder']);p=read(engine/'workflow/pipeline_process.json');r=read(engine/'engine_process.json')
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['inputsStable'].values())
    assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values()) and not r['cancelled'] and not r['timedOut']
    assert r['request']['maximumSeconds']==180 and r['request']['threads']==2 and r['progress']['Inner_Iter']==1
    assert p['nativeReceiptSHA256']==sha(engine/'engine_process.json') and p['auditSHA256']==sha(engine/'thermodynamics/audit.json')
    assert [v['stage'] for v in p['phases']]==['preparation','solver','eosAudit']
    assert all(v['elapsedSeconds']<=v['maximumSeconds'] for v in p['phases'])
    results=[]
    for folder,mat,var,table_name,native,iteration in [(base/'ui01','short_result.mat','short','short_fields.csv',engine,1),
            (base/'ui02','result.mat','r','fields.csv',ROOT/'outputs/v122_validation/actual04',19)]:
        result=loadmat(str(folder/mat),squeeze_me=True,struct_as_record=False)[var];V=result.V;F=result.F.astype(int)-1
        xyz,f=vtk_tokens(native/'surface_flow.vtk');distance,ids=cKDTree(V).query(xyz)
        assert max(distance)<1e-6 and len(set(ids))==len(V)==9878
        node=np.zeros((len(V),3));node[ids]=np.column_stack([f[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']]);face=node[F].mean(axis=1)
        table=np.genfromtxt(str(folder/table_name),names=True,delimiter=',');assert np.array_equal(table['originalFaceID'],np.arange(1,19753))
        errors={}
        for j,key in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):
            error=float(max(abs(table[key]-face[:,j])/np.maximum(abs(face[:,j]),1)));errors[key]=error;assert error<1e-12
        assert (table['givenWallTemperature_K']==300).all() and np.isnan(table['solvedWallTemperature_K']).all() and not table['productionQualified'].any()
        meta=read(folder/(table_name+'.json'));audit=meta['convergenceAudit']
        assert audit['checkpointIdentityConfirmed'] and audit['checkpointIteration']==iteration and audit['sourceCumulativeIteration']==799
        assert 'cumulativeCheckpointIteration' not in audit and not audit['checkpointResidualPassed']
        assert meta['flight']['referenceMach']==5 and abs(meta['flight']['Mach']-4.998034974966865)<1e-12
        assert meta['algorithmRuntime']['flavor']=='rah-coolprop-recon-v3'
        if iteration==19:
            with (folder/'section.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.reader(stream))[1:]
            assert len(rows)==92
            for row in rows:
                assert abs(float(row[0])+1.5)<1e-14 and abs(float(row[10])-face[int(row[2])-1,1])<1e-9*max(abs(float(row[10])),1)
            probe=read(ROOT/'docs/validation/v123/native_ui/observations.json');k=probe['nativeProbeFace']-1
            assert abs(face[k,0]-probe['probePressurePa'])<.0001 and abs(face[k,1]-probe['probeHeatFluxWm2'])<.0001
        results.append(dict(nativeFolder=str(native),localIteration=iteration,originalFaces=len(F),maximumRelativeErrors=errors))
    matlab=[]
    for name,expected in [('heos_ui01',1),('heos_ui02',0)]:
        a=read(base/(name+'_process.json'));assert a['status']=='terminal' and a['exitCode']==expected and all(a['sourceStableAfter'].values())
        matlab.append(dict(name=name,exitCode=expected,sourceStable=True))
    processes=[]
    for path in sorted((ROOT/'outputs/embedded_workbench').glob('*/workflow/pipeline_process.json')):
        a=read(path)
        if a['startedUTC']<p['startedUTC'] and a['startedUTC'][:10]!=p['startedUTC'][:10]:continue
        # 本轮只有三个后台任务；旧输出没有此新流程文件。
        assert a['status']=='terminal'
        if not a['engineStarted']:assert not (path.parent.parent/'engine_process.json').exists()
        processes.append(dict(folder=str(path.parent.parent),exitCode=a['exitCode'],engineStarted=a['engineStarted'],stage=a['stage']))
    assert sorted(q['exitCode'] for q in processes)==[0,1,125]
    report=dict(passed=True,matlab=matlab,pipelines=processes,nativeEngineRuns=1,engineSeconds=r['elapsedSeconds'],
        pipelinePhases=p['phases'],fieldChecks=results,nativeProbeVerified=True,sectionRows=92,
        physicsQualified=False,newQualifiedArea=0,sources={str(q):sha(q) for q in [Path(__file__).resolve(),engine/'engine_process.json',
        engine/'workflow/pipeline_process.json',base/'ui01/short_fields.csv',base/'ui02/fields.csv',base/'ui02/section.csv']})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
