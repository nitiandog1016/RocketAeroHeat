"""独立复核v106外域对照、真实3D启动与原面映射证据，不重新求解。"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.io import loadmat
from auditEmbeddedSU2FlatPlate import read_wall

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def points(path):
    lines=path.read_text(encoding='ascii').splitlines()
    index=next(i for i,line in enumerate(lines) if line.startswith('NPOIN='))
    count=int(lines[index].split('=')[1])
    return np.array([[float(v) for v in line.split()[:2]] for line in lines[index+1:index+1+count]])


def main():
    base=ROOT/'outputs/v105_validation/plate03'
    outer=ROOT/'outputs/v106_validation/plate_outer01'
    p0,p1=points(base/'mesh.su2'),points(outer/'mesh.su2')
    assert len(p1)>len(p0) and np.array_equal(p0,p1[:len(p0)])
    x0,w0=read_wall(base/'surface_flow.vtk');x1,w1=read_wall(outer/'surface_flow.vtk')
    assert np.array_equal(x0,x1)
    mask=(-x0[:,0]>=.04)&(-x0[:,0]<=.16)
    changes={name:float(np.max(np.abs(w1[name][mask]/w0[name][mask]-1))) for name in ['Pressure','Heat_Flux']}
    audit=json.loads((ROOT/'docs/validation/v106/plate_outer_audit.json').read_text(encoding='utf-8'))
    assert audit['residualConverged'] and changes['Pressure']<.01 and changes['Heat_Flux']<.05
    fine=ROOT/'outputs/v106_validation/plate_fine01'
    process=json.loads((fine/'engine_process.json').read_text())
    assert process['timedOut'] and not (fine/'surface_flow.vtk').exists()
    history=np.genfromtxt(str(fine/'history.csv'),delimiter=',',names=True)
    # 被终止任务的末条监测量仅是进程证据，不用于面热流误差或网格无关验收。
    lastFine={key:float(history[key][-1]) for key in ['Inner_Iter','rmsRho','rmsRhoE','HF','maxHF']}

    geometry=loadmat(str(ROOT/'outputs/v106_validation/rocket_volume02/volume_and_source.mat'),squeeze_me=True,struct_as_record=False)
    v,surface=geometry['volume'],geometry['surface'];V=v.vertices
    area=0.;bodyVol=0.;edge=[]
    for face in v.bodyWall:
        ids=np.asarray(face,dtype=int)-1
        edge.extend([tuple(sorted((ids[i],ids[(i+1)%len(ids)]))) for i in range(len(ids))])
        for j in range(1,len(ids)-1):
            a,b,c=V[ids[[0,j,j+1]]]
            area+=np.linalg.norm(np.cross(b-a,c-a))/2
            bodyVol+=a@np.cross(b,c)/6
    _,counts=np.unique(np.asarray(edge),axis=0,return_counts=True)
    assert np.all(counts==2)
    assert abs(area/surface.area.sum()-1)<2e-12
    assert np.all(v.cellVolume>0) and len(v.cells)==11712
    resultDir=ROOT/'outputs/v106_validation/rocket_projection02'
    report=json.loads((resultDir/'report.json').read_text(encoding='utf-8'))
    state=loadmat(str(resultDir/'projection_and_result.mat'),squeeze_me=True,struct_as_record=False)
    projection,result=state['projection'],state['result']
    assert report['originalFaceCount']==len(surface.F)==2448
    assert np.max(np.abs(np.asarray(projection.weights.sum(axis=1)).ravel()-1))<2e-8
    assert projection.weights.data.min()>=0
    physicalPower=float(result.faceArea@result.heatFluxRaw)
    nativePower=float(result.faceArea@result.conservativeHeatLoadDensityRaw)
    assert abs(nativePower/float(report['engineHeatPowerRaw'])-1)<2e-6
    paths=[base/'surface_flow.vtk',base/'mesh.su2',outer/'surface_flow.vtk',outer/'mesh.su2',fine/'engine_process.json',
           fine/'history.csv',resultDir/'report.json',resultDir/'projection_and_result.mat']
    result=dict(status='limited-sensitivity-and-interface-checks-only',
                oldNearWallNodesExactlyPreserved=len(p0),outerNodeCount=len(p1),outerHeight=float(p1[:,1].max()),
                outerRelativeChanges=changes,comparisonInterval=[.04,.16],
                thresholds=dict(pressure=.01,heatFlux=.05),outerScreenPassed=True,
                fineMeshTimedOut=True,fineMeshLastRecorded=lastFine,fineMeshComparisonAvailable=False,
                threeD=dict(cells=len(v.cells),originalFaces=len(surface.F),bodyArea=float(area),bodyVolume=float(bodyVol),
                            wallClosed=True,physicalAreaPower=physicalPower,nativeMappedHeatLoadPower=nativePower,
                            SU2HistoryHF=report['engineHeatPowerRaw'],
                            areaConventionPowerRelativeDifference=physicalPower/float(report['engineHeatPowerRaw'])-1,
                            mappingReport=report),
                newQualifiedRocketArea=0,meshIndependenceQualified=False,publicExperimentCompared=False,
                sourceSHA256=sha(Path(__file__)),inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in paths})
    dest=ROOT/'docs/validation/v106/independent_audit.json';assert not dest.exists()
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
