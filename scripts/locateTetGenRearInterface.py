"""将五个剩余PLC坏单元的固定点追溯到原列/层，不修改几何或求解结果。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from prepareTetGenRuntime import ROOT,sha
from cgalTetrahedralContract import faces_and_neighbors

def main():
    folder=ROOT/'outputs/v143_validation/actual'
    paths=[Path(__file__).resolve(),ROOT/'scripts/cgalTetrahedralContract.py',folder/'hybrid_volume.mat',folder/'candidate_report.json']
    locked={str(p):sha(p) for p in paths};d=loadmat(str(folder/'hybrid_volume.mat'));r=json.loads((folder/'candidate_report.json').read_text())
    P=d['vertices'];T=d['tetrahedra'].astype(np.int64)-1;V=d['sourceVertices'];nv=len(V)
    B,_=faces_and_neighbors(T);fixed=set(np.unique(B));prepared=d['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1
    rows=[]
    for row in np.array(r['badTetRowsOneBased'])-1:
        anchors=[]
        for i in T[row]:
            if i not in fixed:continue
            j=int(prepared[i]);layer,col=divmod(j,nv);assert layer==32
            anchors.append(dict(hybridNodeOneBased=int(i+1),preparedNodeOneBased=j+1,layer=layer,columnOneBased=col+1,
                coordinateM=P[i].tolist(),originalWallCoordinateM=V[col].tolist(),normalHeightM=float(d['normalHeightsByColumnM'][layer,col])))
        rows.append(dict(tetRowOneBased=int(row+1),fixedNodes=anchors))
    report=dict(passed=True,scope='geometric identity only; not a unique-cause or infeasibility proof',
        allBadTetFixedVerticesAtLayer32=True,allBadTetFixedVerticesFromAftWall=True,
        originalWallXRangeM=[min(a['originalWallCoordinateM'][0] for r in rows for a in r['fixedNodes']),
            max(a['originalWallCoordinateM'][0] for r in rows for a in r['fixedNodes'])],
        constraints=rows,sourceSHA256=locked,sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},
        physicsQualified=False,newQualifiedArea=0)
    assert report['originalWallXRangeM'][1]<-2.8 and all(report['sourceStableAfter'].values())
    dst=folder/'rear_interface.json';assert not dst.exists();dst.write_text(json.dumps(report,indent=2)+'\n')
    print(report['originalWallXRangeM'],len(rows),flush=True)

if __name__=='__main__':main()
