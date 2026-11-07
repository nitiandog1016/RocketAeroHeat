"""独立全域核对内点移动：边界/连接/有向体积/SICN及所有受影响单元。"""
from pathlib import Path
import json,sys,time,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareOpenFoamControl import save,sha
from openFoamConservation import geometry
from auditOpenFoamRocketGeometry import internal_area_determinant


def main():
    start=time.monotonic();base=ROOT/'outputs/v200_validation';out=base/'rocket_relocation_audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py');case=base/'rocket03'
    report=json.loads((base/'rocket_relocation01/report.json').read_text());assert report['passed']
    for p,h in report['inputSHA256'].items():assert sha(p)==h
    old=geometry(base/'rocket02');g=geometry(case);free=np.array(report['freeVerticesZeroBased']);fixed=np.ones(len(g['points']),bool);fixed[free]=False
    assert np.array_equal(g['points'][fixed],old['points'][fixed])
    for k in ['cells','faces','owner','neighbor']:assert np.array_equal(g[k],old[k]),k
    assert g['patches']==old['patches'];assert np.array_equal(g['sf'][len(g['neighbor']):],old['sf'][len(g['neighbor']):])
    affected=np.flatnonzero(np.isin(g['cells'],free).any(axis=1));assert affected.tolist()==report['affectedCellsZeroBased']
    xyz=g['points'][g['cells']];A=(xyz[:,1:]-xyz[:,:1]).transpose(0,2,1);det=np.linalg.det(A);assert np.all(det>0)
    equilateral=np.array([[1,.5,.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2/3)]])
    J=A@np.linalg.inv(equilateral);cof=np.stack([np.cross(J[:,:,1],J[:,:,2]),np.cross(J[:,:,2],J[:,:,0]),np.cross(J[:,:,0],J[:,:,1])],axis=2)
    quality=3*np.linalg.det(J)/np.sqrt(np.sum(J*J,axis=(1,2))*np.sum(cof*cof,axis=(1,2)))
    metric,counts=internal_area_determinant(g['sf'],g['owner'],g['neighbor'],len(g['cells']))
    assert min(quality)>=.1 and min(metric)>=.001
    assert abs(det.sum()/6/381.4346146072696-1)<2e-13
    mapping=json.loads((case/'wall_mapping.json').read_text());owners=np.array(mapping['nativeOwnerCellsZeroBased']);ids=np.array(mapping['nativeFaceIDsZeroBased'])
    assert len(ids)==19752 and np.array_equal(g['owner'][ids],owners)
    assert np.array_equal(g['centers'][owners],mapping['samplePositionsM']) and min(mapping['normalDistanceM'])>0
    assert abs(g['area'][ids].sum()-1.9646292565715553)<2e-13
    for name in ['constant/polyMesh/faces','constant/polyMesh/owner','constant/polyMesh/neighbour','constant/polyMesh/boundary','constant/physicalProperties','0/p','0/T','0/U']:
        assert sha(case/name)==sha(base/'rocket02'/name)
    result=dict(passed=True,originalWallExact=True,originalFarfieldExact=True,connectivityExact=True,originalWallFaces=19752,
        originalWallAreaM2=float(g['area'][ids].sum()),changedInteriorVertices=4,affectedCells=len(affected),minimumSICN=float(quality.min()),
        minimumDeterminant=float(metric.min()),badCells=0,volumeM3=float(det.sum()/6),minimumVolumeM3=float(det.min()/6),
        geometry=dict(cells=len(g['cells']),points=len(g['points']),faces=len(g['faces']),internalFaces=len(g['neighbor']),volumeM3=float(det.sum()/6)),
        inputSHA256=report['inputSHA256'],originalRocketCFDRuns=0,executionAuthorized=False,physicalQualification=False,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','minimumSICN','minimumDeterminant','volumeM3','elapsedSeconds']}))


if __name__=='__main__':main()
