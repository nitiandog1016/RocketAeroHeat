"""复算已存下限点的壁列和层号，绑定原比较/初猜/网格；不推断物理原因。"""
from pathlib import Path
import json
import numpy as np
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    base=ROOT/'outputs/v130_validation';path=base/'remaining_node_layers_audit.json';assert not path.exists()
    comparison=base/'column_comparison.json';values=base/'column_final/column_values.npz'
    mesh=ROOT/'outputs/v120_validation/mesh32/mesh_report.json'
    r=read(comparison);anchors=np.load(str(values))['anchorLayers'];heights=read(mesh)['normalHeightsM'];rows={}
    assert len(anchors)==9878 and r['originalFaces']==19752
    for step,data in r['variants']['column']['snapshots'].items():
        rows[step]={}
        for key in ['originalFloorNodesRemaining','newFloorNodeIDs']:
            rows[step][key]=[dict(nodeID=i,wallParentID=int(i%9878),normalLayer=int(i//9878),
                columnAnchorLayer=int(anchors[i%9878]),withinModifiedInterior=bool(0<i//9878<anchors[i%9878]),
                normalDistanceM=heights[i//9878]) for i in data[key]]
    original=base/'remaining_node_layers.json';assert rows==read(original)
    result=dict(passed=True,identicalToFirstReadonlyReview=True,layers=rows,physicalCauseEstablished=False,
        sourceSHA256={str(p):sha(p) for p in [comparison,values,mesh,original,Path(__file__).resolve()]})
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Stored layer assignments reproduced; no solver run.')


if __name__=='__main__':main()
