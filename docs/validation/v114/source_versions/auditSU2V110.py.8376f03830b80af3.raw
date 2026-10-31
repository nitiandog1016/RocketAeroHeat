"""独立流式读取SU2 VTK，核对空间残差统计、实际审核与主工作台导出。"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v110_validation'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def vtk_tokens(path):
    # 与生产诊断的正则块解析独立：按VTK词法依次消耗点、单元和字段。
    with path.open() as stream:
        assert stream.readline().startswith('# vtk')
        stream.readline(); assert stream.readline().strip()=='ASCII'
        assert stream.readline().strip()=='DATASET UNSTRUCTURED_GRID'
        words=iter(stream.read().split())
    assert next(words)=='POINTS'; n=int(next(words)); next(words)
    xyz=np.fromiter((float(next(words)) for _ in range(3*n)),float).reshape(n,3)
    assert next(words)=='CELLS'; nc=int(next(words)); size=int(next(words))
    for _ in range(size):next(words)
    assert next(words)=='CELL_TYPES' and int(next(words))==nc
    for _ in range(nc):next(words)
    assert next(words)=='POINT_DATA' and int(next(words))==n
    fields={}
    for tag in words:
        name=next(words);next(words)
        if tag=='SCALARS':
            width=int(next(words));assert next(words)=='LOOKUP_TABLE' and next(words)=='default'
        else:assert tag=='VECTORS';width=3
        fields[name]=np.fromiter((float(next(words)) for _ in range(n*width)),float).reshape(n,width)
    return xyz,fields


def main():
    paths=[OUT/'spatial_report_final.json',OUT/'probe02/volume_flow.vtk',OUT/'probe02/history.csv',
           OUT/'audit05/timeout_audit.json',OUT/'audit05/plate_audit.json',OUT/'audit05/negative_cases.json',
           OUT/'audit05/report.json',OUT/'audit05/current_fields.csv.json',OUT/'audit05/current_section.csv.json']
    report=json.loads(paths[0].read_text());xyz,fields=vtk_tokens(paths[1])
    with paths[2].open() as stream:
        reader=csv.reader(stream,skipinitialspace=True);keys=[x.strip() for x in next(reader)];row=dict(zip(keys,map(float,next(reader))));assert next(reader,None) is None
    residual=np.column_stack((fields['Residual_Density'],fields['Residual_Momentum'],fields['Residual_Energy']))
    columns=['rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]']
    errors=[abs(np.log10(np.linalg.norm(residual[:,i])/np.sqrt(len(xyz)))-row[k]) for i,k in enumerate(columns)]
    assert max(errors)<2e-6
    squared=residual[:,-1]**2;fraction=float(squared[xyz[:,0]<0].sum()/squared.sum())
    assert abs(fraction-report['downstreamOfBody']['rmsRhoE'])<1e-12
    assert np.argmax(squared)==report['topEnergyNodes'][0]['node']
    assert abs(fields['Orthogonality'].min()-report['minimumOrthogonality']['valueDeg'])<1e-12
    assert abs(sum(report[k]['rmsRhoE'] for k in ('wall','normalShell','outer'))-1)<2e-12
    timeout=json.loads(paths[3].read_text());plate=json.loads(paths[4].read_text());negative=json.loads(paths[5].read_text())
    assert timeout['checkpointIteration']==3750 and timeout['terminalIteration']==3921 and timeout['checkpointIdentityConfirmed']
    assert not timeout['terminalResidualPassed'] and not timeout['checkpointResidualPassed']
    assert plate['terminalResidualPassed'] and plate['terminalIteration']==23477 and plate['checkpointIteration'] is None
    assert plate['heatLoadCandidateCount']==17 and not plate['checkpointIdentityConfirmed']
    assert negative['changedOutputRejected']['status']=='evidence-rejected'
    assert negative['sparseHistoryRejected']['status']=='evidence-rejected'
    assert negative['duplicateHeatLoadUnconfirmed']['heatLoadCandidateCount']>=2
    ui=json.loads(paths[6].read_text(encoding='utf8'))
    for p in paths[7:]:assert json.loads(p.read_text(encoding='utf8'))['convergenceAudit']==ui['audit']
    assert ui['audit']['checkpointIteration']==1 and not ui['audit']['checkpointResidualPassed']
    assert not ui['solverRun'] and ui['rendered']['webglError']==0
    for p,h in report['inputSHA256'].items():assert sha(ROOT/p)==h
    out=ROOT/'docs/validation/v110/independent_audit.json';assert not out.exists()
    result=dict(nativeResidualRMSMaximumLogDifference=max(errors),downstreamIntegratedEnergySquareFraction=fraction,
                nativeNodes=len(xyz),checkpointIteration=3750,terminalIteration=3921,
                convergedPlateAmbiguousCandidates=17,mainUISnapshot='v107 M4 alpha20 beta5, iteration1',
                fieldAndSectionMetadataMatch=True,physicsQualified=False,newQualifiedArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths+[Path(__file__).resolve()]})
    out.write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
