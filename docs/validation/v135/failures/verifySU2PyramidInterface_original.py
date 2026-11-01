"""独立核对MATLAB混合单元、SU2原连接、中位体积及两类壁距；不授予物理解。"""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from analyzeSU2DualVolumes import dual_volumes
from su2SSTInterface import mesh,wall_geometry
from prepareSU2CheckpointContinuation import sha

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--folder',type=Path,required=True);a=parser.parse_args()
    folder=a.folder.resolve();v=loadmat(str(folder/'interface.mat'),squeeze_me=True,struct_as_record=False)['volume']
    dual,check=dual_volumes(v);points,cells,wall=mesh(folder/'engine/mesh.su2')
    assert np.array_equal(points,v.vertices)
    for kind,ids,row in zip(v.cellTypes,v.cells,cells):
        expected=np.asarray(ids,dtype=int)-1
        if kind==13:expected=expected[[0,2,1,3,5,4]]
        assert row[0]==kind and np.array_equal(row[1:1+len(expected)],expected)
    assert len(cells)==len(v.cells)==92
    assert np.array_equal(np.array([row[1:] for row in wall]),v.sourceFaces[:,[0,2,1]]-1)
    g=wall_geometry(folder/'engine/mesh.su2');assert len(g['wallIDs'])==6
    assert np.all(g['thermalDistance']>0) and np.all(g['nearestWallDistance']>0)
    log=(folder/'engine/engine_stdout.log').read_text('utf8')
    assert 'PYRAMID volume elements' not in log and 'volume elements are distorted' not in log
    report=dict(dualGeometry=check,minimumNodeDualVolume=float(min(dual)),exactExportedNodes=True,
        exactExportedConnectivity=True,wallVertices=len(g['wallIDs']),minimumThermalDistance=float(min(g['thermalDistance'])),
        minimumNearestWallDistance=float(min(g['nearestWallDistance'])),pyramidReorientation=False,
        qualityQualified=False,physicsQualified=False,newQualifiedArea=0,
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [folder/'interface.mat',folder/'engine/mesh.su2',
            folder/'engine/engine_stdout.log',Path(__file__).resolve(),ROOT/'scripts/analyzeSU2DualVolumes.py',ROOT/'scripts/su2SSTInterface.py']})
    (folder/'independent.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
