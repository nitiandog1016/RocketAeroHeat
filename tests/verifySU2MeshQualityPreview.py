"""独立SVD校核真实体单元预览与来源，绘制局部科研图，不重跑网格器或MATLAB。"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    base=ROOT/'outputs/v134_validation';out=base/'quality_check.json';assert not out.exists()
    source=ROOT/'outputs/v133_validation/fan_workbench/hybrid_volume.mat'
    old=loadmat(str(source));P=old['vertices'];T=old['tetrahedra'].astype(int)-1;npw=len(old['wedges'])
    data=loadmat(str(base/'workbench02/quality_preview.mat'),squeeze_me=True,struct_as_record=False)
    info=data['info'];payload=data['payload'];ids=np.atleast_1d(info.sourceVolumeNodeIDs).astype(int)-1
    # 以奇异值定义条件数质量，独立于MATLAB的余子式实现。
    ref=np.array([[1,.5,.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2/3)]])
    p=P[T];J=np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=2)@np.linalg.inv(ref)
    singular=np.linalg.svd(J,compute_uv=False);q=3/np.sqrt(np.sum(singular**2,axis=1)*np.sum(singular**-2,axis=1))
    assert np.all(np.linalg.det(J)>0);bad=np.flatnonzero(q<.1)
    assert len(bad)==17 and np.array_equal(info.sourceTetRows.astype(int)-1,bad)
    assert np.array_equal(info.sourceCellIDs.astype(int)-1,bad+npw)
    error=float(max(abs(q[bad]-info.cellSICN)));assert error<1e-12
    assert np.array_equal(payload.vertices,P[ids])
    faces=ids[payload.faces.astype(int)];n=int(info.visibleCellFaces);owners=info.displaySourceCellIDs.astype(int)
    for tri,owner,value in zip(faces[:n],owners[:n],info.displayValues[:n]):
        row=owner-1-npw;assert set(tri)<=set(T[row]);assert abs(value-q[row])<1e-12
    original=old['sourceFaces'].astype(int)-1;context=info.contextOriginalFaceIDs.astype(int)-1
    assert np.array_equal(faces[n:],original[context]) and np.all(owners[n:]==0)
    assert np.all(np.isnan(payload.faceField.values[n:]))
    task=json.loads((base/'mesh_quality_workbench02_process.json').read_text(encoding='utf-8'))
    assert task['exitCode']==0 and all(task['sourceStableAfter'].values())
    # 科学图来自真实节点和连接，明确是质量诊断而非物理场云图或UI截图。
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    fig=plt.figure(figsize=(8.5,5.2),constrained_layout=True);ax=fig.add_subplot(111,projection='3d')
    ax.add_collection3d(Poly3DCollection(P[faces[n:]],facecolors='#d9dfe2',edgecolors='none',alpha=.45))
    coll=Poly3DCollection(P[faces[:n]],edgecolors='#273e50',linewidths=.4,cmap='viridis_r')
    coll.set_array(info.displayValues[:n]);coll.set_clim(0,.1);ax.add_collection3d(coll)
    limits=P[ids];lo=limits.min(axis=0);hi=limits.max(axis=0)
    ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_zlim(lo[2],hi[2]);ax.set_box_aspect(hi-lo)
    ax.set_xlabel('X / m');ax.set_ylabel('Y / m');ax.set_zlabel('Z / m');ax.view_init(22,-54)
    ax.set_title('鼻部外层 · 17个真实低质量四面体\n灰色为原壁面；仅网格诊断，未增加气动热资格',fontsize=11)
    fig.colorbar(coll,ax=ax,shrink=.62,pad=.09,label='SICN / 1（越小越差）')
    plot=base/'nose_quality.png';fig.savefig(str(plot),dpi=170);plt.close(fig)
    report=dict(passed=True,lowQualityCells=17,visibleCellFaces=n,contextFaces=len(context),
        minimumSICN=float(min(q)),maximumSVDQualityDifference=error,exactCoordinatesAndConnectivity=True,
        nativeFlowRuns=0,newQualifiedArea=0,plotIsNativeScreenshot=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),source,base/'workbench02/quality_preview.mat',
            base/'mesh_quality_workbench02_process.json']},imageSHA256=sha(plot))
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'},indent=2))

if __name__=='__main__':main()
