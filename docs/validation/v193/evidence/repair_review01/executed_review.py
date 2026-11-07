# -*- coding: utf-8 -*-
"""只读已结束的原生refine，定位边界连接变化；不再调用网格器或CFD。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import numpy as np
from tetgenFixedInterface import read_native,tet_audit
from readEilmerNative import sha


def review(base,out):
    start=time.monotonic();base,out=map(Path,(base,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_review.py')
    src=base/'mesh01/mesh.npz';raw=base/'repair01/native/native_result.dat'
    m=np.load(src);P=m['vertices'];F=m['wall'];G=m['farfield']
    Q,T,B,refs,marks=read_native(raw);actual,geometry=tet_audit(Q,T)
    assert np.all(refs[B]>0) and np.all(marks==1)
    fixed=refs>0;delta=float(np.max(abs(Q[fixed]-P[refs[fixed]-1])))
    mapped=refs[B]-1;newkeys={tuple(sorted(f)) for f in mapped}
    oldkeys={tuple(sorted(f)) for f in np.r_[F,G]}
    removed_wall=[i+1 for i,f in enumerate(F) if tuple(sorted(f)) not in newkeys]
    removed_far=[i+1 for i,f in enumerate(G) if tuple(sorted(f)) not in newkeys]
    added=np.array([f for f in mapped if tuple(sorted(f)) not in oldkeys],dtype=int).reshape(-1,3)
    np.savetxt(out/'changed_boundary_triangles.csv',added,delimiter=',',fmt='%d',header='source_node_0,source_node_1,source_node_2',comments='')
    result=dict(passed=True,scope='read-only diagnosis of rejected native result',nativeGeometry=geometry,
        nativeNodes=len(Q),nativeTets=len(T),fixedCoordinateMaxDifferenceM=delta,
        removedOriginalWallFaceIDs=removed_wall,removedOriginalFarfieldFaceIDs=removed_far,
        addedBoundaryFaces=len(added),originalBoundaryPreserved=not removed_wall and not removed_far and len(added)==0,
        candidateAccepted=False,physicalQualification=False,sourceSHA256={str(p):sha(p) for p in [src,raw,Path(__file__)]})
    assert not result['originalBoundaryPreserved']
    # 模块二：原箭体与6个低质单元位置，固定空间比例；图仅说明几何失败位置。
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    bad=json.loads((base/'geometry01/report.json').read_text())['badCells']
    fig=plt.figure(figsize=(11,4.5));ax=fig.add_subplot(111,projection='3d')
    shell=Poly3DCollection(P[F],facecolors='#bfc7d0',edgecolors='none',alpha=.65);ax.add_collection3d(shell)
    centers=np.array([r['centerM'] for r in bad]);ax.scatter(centers[:,0],centers[:,1],centers[:,2],c='#b3261e',s=28,label='6 cells below SICN 0.1',depthshade=False)
    ax.set(xlim=(-2.95,.15),ylim=(-.22,.22),zlim=(-.22,.22),xlabel='X [m]',ylabel='Y [m]',zlabel='Z [m]')
    ax.set_box_aspect((3.1,.44,.44));ax.view_init(elev=19,azim=-64);ax.legend(loc='upper left',fontsize=8)
    ax.set_title('Original 19,752-wall rocket: Euler mesh geometry screening\n69,177 tetrahedra; candidate rejected; no new flow / heat solution',fontsize=10)
    fig.savefig(out/'rocket_mesh_screen.png',dpi=160,bbox_inches='tight');plt.close(fig)
    result['elapsedSeconds']=time.monotonic()-start
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','fixedCoordinateMaxDifferenceM','removedOriginalWallFaceIDs','addedBoundaryFaces','elapsedSeconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();review(a.base,a.out)
