"""独立核对工作台截面的实际体节点/面来源，并画可导出的科研几何图。"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import read,sha

def main():
    base=ROOT/'outputs/v133_validation';out=base/'slice_check.json';assert not out.exists()
    r=read(base/'fan_workbench/mesh_report.json');d=loadmat(str(base/'fan_workbench/hybrid_volume.mat'))
    result=loadmat(str(base/'workbench01/slice.mat'),squeeze_me=True,struct_as_record=False)
    info=result['info'];payload=result['payload'];ids=np.atleast_1d(info.sourceVolumeNodeIDs).astype(int)-1
    tris=np.atleast_2d(info.sourceVolumeTriangles).astype(int)-1
    P=d['vertices'];W=d['wedges'].astype(int)-1
    assert int(info.sourceVertexID)==18 and len(tris)==62
    eligible=np.unique(np.r_[np.array(r['sideNodeIDs'])[:,17],np.array(r['baseNodeIDs'])[:,17]])
    assert np.array_equal(ids,eligible)
    tri=np.concatenate([W[r['originalPrisms']:,:3],W[r['originalPrisms']:,3:]])
    tri=tri[np.all(np.isin(tri,eligible),axis=1)]
    assert np.array_equal(np.unique(np.sort(tri,axis=1),axis=0),np.unique(np.sort(tris,axis=1),axis=0))
    assert np.array_equal(payload.vertices,P[ids]) and np.array_equal(ids[payload.faces.astype(int)],tris)
    areas=np.linalg.norm(np.cross(P[tris[:,1]]-P[tris[:,0]],P[tris[:,2]]-P[tris[:,0]]),axis=1)/2
    assert abs(sum(areas)-info.area)<1e-15 and min(areas)>0
    task=read(base/'crease_fan_workbench01_process.json')
    assert task['status']=='terminal' and task['exitCode']==0 and all(task['sourceStableAfter'].values())
    # 局部正交坐标仅用于几何图，显示的所有边都来自实际棱柱横截面。
    sid=np.array(r['sideNodeIDs'])[:,17];bid=np.array(r['baseNodeIDs'])[:,17]
    origin=P[sid[0]];u=P[sid[1]]-origin;u/=np.linalg.norm(u)
    v=P[bid[1]]-origin;v-=u*(v@u);v/=np.linalg.norm(v)
    planarError=float(max(abs((P[ids]-origin)@np.cross(u,v))));assert planarError<1e-12
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    fig,axs=plt.subplots(1,2,figsize=(9,4.3),constrained_layout=True)
    basis=np.column_stack([u,v]);xy=(P-origin)@basis
    edge=np.unique(np.sort(np.concatenate([tris[:,[0,1]],tris[:,[1,2]],tris[:,[2,0]]]),axis=1),axis=0)
    for ax,scale,title,unit in [(axs[0],1e3,'完整过渡截面','mm'),(axs[1],1e6,'第一层附近','µm')]:
        ax.add_collection(LineCollection(xy[edge]*scale,colors='#62717d',linewidths=.65))
        ax.plot(*(xy[sid]*scale).T,color='#176b87',lw=1.3,label='侧壁法向列')
        ax.plot(*(xy[bid]*scale).T,color='#bd7134',lw=1.3,label='底面法向列')
        ax.set_aspect('equal');ax.autoscale();ax.set_title(title,fontsize=11)
        ax.set_xlabel('局部坐标 1 / '+unit);ax.set_ylabel('局部坐标 2 / '+unit)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    axs[0].legend(frameon=False,fontsize=8);axs[1].set_xlim(-1,30);axs[1].set_ylim(-1,30)
    fig.suptitle('尾缘原壁点18 · 62个实际三角连接 · 仅几何，未获气动热资格',fontsize=12)
    image=base/'crease_fan_section.png';fig.savefig(str(image),dpi=180);plt.close(fig)
    report=dict(passed=True,faces=62,nodes=len(ids),sourceWallVertexID=18,coordinateAndFaceArraysExactlyMatched=True,
        crossSectionAreaM2=float(sum(areas)),maximumOutOfPlaneDistanceM=planarError,
        geometryOnly=True,nativeFlowRuns=0,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),base/'workbench01/slice.mat',
            base/'fan_workbench/hybrid_volume.mat',base/'fan_workbench/mesh_report.json',base/'crease_fan_workbench01_process.json']},
        imageSHA256=sha(image))
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'},indent=2))

if __name__=='__main__':main()
