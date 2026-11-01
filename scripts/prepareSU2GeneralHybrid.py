"""把已验收的完整PLC混合网格接回工作台；只转换元数据，不重生成几何。"""
from pathlib import Path
import json,os,time
import numpy as np
from scipy.io import loadmat,savemat
from prepareSU2CheckpointContinuation import ROOT,sha,read

METHOD='original-triangle-general-prisms-and-tetgen-plc-exterior'

def prepare(source,target):
    source=Path(source).resolve();target=Path(target).resolve();target.relative_to(ROOT/'outputs')
    assert not target.exists();start=time.monotonic()
    r=read(source/'mesh_process.json');c=read(source/'candidate_report.json');g=read(source/'geometry_audit.json')
    proof=read(source/'independent_review.json')
    assert r['status']=='terminal' and r['exitCode']==0 and r['sourceStable'] and all(r['sourceStableAfter'].values())
    assert c['candidateAccepted'] and c['geometryPassed'] and c['qualityGatePassed'] and c['sizeGatePassed']
    assert proof['passed'] and proof['geometryCandidateAccepted'] and not proof['physicsQualified']
    for name in ['hybrid_volume.mat','candidate_report.json','geometry_audit.json','preparation_report.json']:
        assert sha(source/name)==r['outputSHA256'][name]
    fan=ROOT/'outputs/v133_validation/fan_workbench';old=read(fan/'mesh_report.json')
    assert sha(fan/'mesh_report.json')==r['sourceSHA256'][str(fan/'mesh_report.json')]
    locks=[source/n for n in ['hybrid_volume.mat','mesh_process.json','candidate_report.json','geometry_audit.json','independent_review.json','preparation_report.json']]
    locks += [fan/'mesh_report.json',fan/'mesh_input.mat',Path(__file__).resolve()]
    locked={str(p):sha(p) for p in locks}
    d=loadmat(str(source/'hybrid_volume.mat'));ids=d['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1
    valid=ids>=0;assert len(np.unique(ids[valid]))==sum(valid)
    inverse=np.full(int(max(ids))+1,-1,np.int64);inverse[ids[valid]]=np.flatnonzero(valid)
    report=dict(old)
    for name in ['sideNodeIDs','baseNodeIDs']:
        mapped=inverse[np.asarray(old[name],np.int64)];assert mapped.min()>=0;report[name]=mapped.tolist()
    report.update(method=METHOD,coordinateMode=str(d['columnCoordinateMode'][0]),normalLayers=32,
        normalHeightMode='per-column reference coordinates, not path length',geometryAccepted=True,
        hasTangentialOffsets=True,meshSourceDirectory=str(source),meshSourceSHA256=sha(source/'hybrid_volume.mat'),
        nodes=c['nodes'],prisms=c['prisms'],tetrahedra=c['tetrahedra'],pyramids=0,
        minimumTetSICN=c['minimumTetSICN'],tetsBelowSICN01=c['tetsBelowSICN01'],
        bodyArea=c['originalWallAreaM2'],bodyVolume=g['bodyVolume'],farfieldVolume=g['outerVolume'],
        minimumTetVolume=g['minimumTetVolume'],minimumPrismJacobian=g['minimumPrismJacobian'],
        exteriorAlgorithm=c['exteriorAlgorithm'],maximumTangentialSlope=proof['maximumOuterColumnTangentialSlope'],
        geometryOnly=True,physicsQualified=False,newQualifiedArea=0)
    report.pop('gmshVersion',None)
    target.mkdir(parents=True);os.link(str(source/'hybrid_volume.mat'),str(target/'hybrid_volume.mat'))
    inp=loadmat(str(fan/'mesh_input.mat'));savemat(str(target/'mesh_input.mat'),dict(V=d['sourceVertices'],F=d['sourceFaces'],
        optionsJSON=inp['optionsJSON'],sourceMeshSHA256=sha(source/'hybrid_volume.mat'),coordinateMode=d['columnCoordinateMode']),do_compression=True)
    (target/'mesh_report.json').write_text(json.dumps(report,indent=2)+'\n')
    stable={p:sha(Path(p))==h for p,h in locked.items()};assert all(stable.values())
    receipt=dict(status='terminal',exitCode=0,timedOut=False,inputStable=True,sourceStable=True,sourceSHA256=locked,
        sourceStableAfter=stable,input=str(target/'mesh_input.mat'),inputSHA256=sha(target/'mesh_input.mat'),
        outputSHA256={n:sha(target/n) for n in ['mesh_report.json','hybrid_volume.mat']},
        elapsedSeconds=time.monotonic()-start,geometryRegenerated=False,physicsQualified=False)
    (target/'mesh_process.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--target',required=True);a=p.parse_args()
    print(json.dumps(prepare(a.source,a.target),indent=2))
