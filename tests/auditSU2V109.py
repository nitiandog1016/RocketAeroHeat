"""独立核对v109初场机制、超时检查点身份及实际截面导出，不追加数值运行。"""
from pathlib import Path
import json
import hashlib
import csv
import io
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from auditEmbeddedSU2FlatPlate import read_wall

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v109_validation'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    folder=OUT/'rocket_seed01';receipt=json.loads((folder/'engine_process.json').read_text())
    assert receipt['status']=='terminal' and receipt['timedOut'] and not receipt['cancelled'] and all(receipt['inputsStable'].values())
    for name,h in receipt['outputSHA256'].items():assert sha(folder/name)==h
    history=np.genfromtxt(str(folder/'history.csv'),delimiter=',',names=True)
    baseline=np.genfromtxt(str(ROOT/'outputs/v108_validation/normal01/native_input/history.csv'),delimiter=',',names=True)
    improvement=float(baseline['rmsRhoE'][1]-history['rmsRhoE'][1]);assert improvement>=.5
    rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1)
    assert rest.shape==(59122,9) and np.array_equal(rest[:,0],np.arange(len(rest)))
    rho=rest[:,4];velocity=rest[:,5:8]/rho[:,None];internal=rest[:,8]/rho-(velocity**2).sum(axis=1)/2
    T=internal*.4/287.05287;pressure=rho*287.05287*T
    assert np.all(np.isfinite(rest)) and np.all(rho>0) and np.all(internal>0)
    xyz,wall=read_wall(folder/'surface_flow.vtk');distance,ids=cKDTree(rest[:,1:4]).query(xyz)
    assert len(np.unique(ids))==len(xyz) and max(distance)<1e-6
    errors=dict(pressure=float(max(abs(wall['Pressure']/pressure[ids]-1))),temperature=float(max(abs(wall['Temperature']/T[ids]-1))))
    assert max(errors.values())<2e-6
    saved=loadmat(str(ROOT/'outputs/v108_validation/normal01/normal_volume.mat'),squeeze_me=True,struct_as_record=False)
    projection=saved['projection'];native_power=float(wall['Heat_Flux']@projection.engineNodeArea[ids])
    # stdout在强制超时前有缓冲未刷出，不能靠日志最后一个写入表猜检查点时刻。
    # 用实际VTK节点热负荷与所有250步写出候选的高精度HF唯一匹配。
    candidates=np.flatnonzero(np.mod(history['Inner_Iter'],250)==0)
    differences=abs(history['HF'][candidates]-native_power)/max(abs(native_power),1)
    matches=candidates[differences<2e-6];assert len(matches)==1,'检查点历史身份不唯一或输出不完整'
    checkpoint=int(matches[0]);assert history['Inner_Iter'][checkpoint]<history['Inner_Iter'][-1]
    mapped=loadmat(str(OUT/'checkpoint01/checkpoint_original_faces.mat'),squeeze_me=True,struct_as_record=False)['result']
    fields=np.zeros((len(rest),3));fields[ids]=np.column_stack([wall[k] for k in ('Pressure','Heat_Flux','Temperature')])
    expected=projection.weights@fields
    actual=np.column_stack((mapped.pressureRaw,mapped.heatFluxRaw,mapped.fluidNodeTemperatureMeanRaw))
    mapping=float(np.max(abs(expected-actual)/np.maximum(1,abs(expected))));assert mapping<2e-12
    assert abs(mapped.conservativeHeatLoadPowerRaw/native_power-1)<2e-12
    columns=['rmsRho','rmsRhoU','rmsRhoV','rmsRhoW','rmsRhoE','HF'];tail=history[-251:]
    heat_change=float((tail['HF'][-1]-tail['HF'][0])/max(abs(tail['HF'][0]),1))
    # 15组完整导出必须逐字匹配网页回执，且每个原面数值/状态来自同一快照。
    records=json.loads((OUT/'sections02/section_records.json').read_text(encoding='utf8'))
    mesh=json.loads((OUT/'sections02/case_for_node.json').read_text(encoding='utf8'));count=0
    for i,record in enumerate(records):
        file=OUT/'sections02'/('%s_x%02d.csv'%(record['fieldKey'],i%3+1))
        assert file.read_bytes().decode('utf8')==record['csv']
        rows=list(csv.reader(io.StringIO(record['csv'].lstrip('\ufeff'))))
        for row in rows[1:]:
            face=int(row[2])-1;key=record['fieldKey'];assert float(row[0])==record['x']
            if key=='wallSolved':assert row[10]=='' and row[11]=='未求解'
            elif key=='wallInput':assert float(row[10])==300 and row[11]=='给定条件'
            else:assert float(row[10])==mesh[key][face] and row[11]=='研究计算值；未授予物理资格'
            count+=1
    final=json.loads((OUT/'final_ui01/report.json').read_text(encoding='utf8'))
    assert (OUT/'final_ui01/final_current_section.csv').read_bytes().decode('utf8')==final['section']['csv']
    paths=[Path(__file__).resolve(),folder/'engine_process.json',folder/'history.csv',folder/'restart_flow.csv',folder/'surface_flow.vtk',
           OUT/'checkpoint01/checkpoint_original_faces.mat',OUT/'sections02/section_records.json',OUT/'sections02/case_for_node.json',
           OUT/'sections02/node_audit.json',OUT/'final_ui01/report.json',OUT/'final_ui01/final_current_section.csv']
    report=dict(startupEnergyLogImprovement=improvement,startupScreenPassed=True,residualConverged=False,
                terminalIteration=int(history['Inner_Iter'][-1]),terminalFields={k:float(history[k][-1]) for k in columns},
                checkpointIteration=int(history['Inner_Iter'][checkpoint]),checkpointFields={k:float(history[k][checkpoint]) for k in columns},
                checkpointMatchRelativeError=float(min(differences)),checkpointIdentityMethod='unique native-area HF match among output-cadence history; EOS verifies restart/VTK consistency',
                checkpointDensityMinimum=float(min(rho)),checkpointTemperatureRange=[float(min(T)),float(max(T))],
                checkpointWallEOSRelativeErrors=errors,originalFaceMappingMaxNormalizedError=mapping,
                last250IterationHeatPowerRelativeChange=heat_change,originalFaces=2448,sectionCSVRows=count,
                actualFinalSliderX=final['section']['x'],physicsQualified=False,newQualifiedRocketArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    dest=ROOT/'docs/validation/v109/independent_audit.json';assert not dest.exists()
    dest.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps({k:v for k,v in report.items() if k!='inputSHA256'},indent=2))


if __name__=='__main__':main()
