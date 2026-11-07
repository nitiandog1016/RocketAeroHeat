# -*- coding: utf-8 -*-
"""原壁完整Euler PLC网格；复用成熟TetGen，旧黏性网格只读且不自动启动CFD。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import traceback
import numpy as np
from scipy.io import loadmat
from tetgenPLCInterface import execute, winding
from reviewSU2ExplicitQuality import quality
from readEilmerNative import sha

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = '22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c'


def source_surface():
    """模块一：原物理壁与旧远场直接取原身份，不用棱柱外包络或最近邻。"""
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'
    proof=ROOT/'outputs/v144_validation/actual/independent_review.json'
    baseline=ROOT/'outputs/v104_validation/baseline02/current_case.mat'
    assert sha(source)==SOURCE_SHA
    review=json.loads(proof.read_text())
    assert review['passed'] and review['geometryCandidateAccepted']
    assert review['sourceSHA256'][str(source)]==SOURCE_SHA
    d=loadmat(str(source),variable_names=['vertices','sourceVertices','sourceFaces','farfield'])
    V=d['sourceVertices'];F=d['sourceFaces'].astype(np.int64)-1
    b=loadmat(str(baseline),simplify_cells=True)['current']
    assert np.array_equal(V,b['V']) and np.array_equal(F,b['F']-1)
    assert V.shape==(9878,3) and F.shape==(19752,3)
    far=d['farfield'].astype(np.int64)-1
    nodes,inverse=np.unique(far,return_inverse=True)
    Q=d['vertices'][nodes];G=inverse.reshape(-1,3)+len(V)
    assert len(Q)==3288 and len(G)==6572
    P=np.r_[V,Q];boundary=np.r_[F[:,[0,2,1]],G]
    tri=V[F];six=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2]))
    hole=np.sum(tri.sum(axis=1)*six[:,None]/4,axis=0)/six.sum()
    assert abs(winding(V,F,hole)-1)<1e-8
    L=np.ptp(V[:,0]);witness=np.array([V[:,0].mean(),Q[:,1].min()+.1*L,0.])
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert abs(area.sum()-1.9646292565715551)<2e-13
    return P,F,G,boundary,hole,witness,area,[source,proof,baseline]


def write_su2(path,P,T,F,G):
    """模块二：原壁顺序作为输出面顺序；边界外法向指向流体域外。"""
    with path.open('w',encoding='ascii',newline='\n') as f:
        f.write('NDIME= 3\nNELEM= %d\n'%len(T))
        for row in T:f.write('10 '+' '.join(map(str,row))+'\n')
        f.write('NPOIN= %d\n'%len(P))
        for i,p in enumerate(P):f.write(' '.join('%.17g'%a for a in p)+' %d\n'%i)
        f.write('NMARK= 2\n')
        for name,faces in [('wall',F[:,[0,2,1]]),('farfield',G)]:
            f.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(name,len(faces)))
            for row in faces:f.write('5 '+' '.join(map(str,row))+'\n')


def prepare(out):
    """模块三：唯一新目录、有限原生预算、实际资源计数和准确失败记录。"""
    start=time.monotonic();out=Path(out).resolve();out.relative_to(ROOT/'outputs')
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_prepare.py')
    result=dict(schema='RAH.HybridRocketEulerMesh.v1',passed=False,
                executionAuthorized=False,physicalQualification=False,newSolvedHeatAreaM2=0)
    paths=[Path(__file__).resolve(),ROOT/'docs/validation/v193/PLAN.md']
    try:
        P,F,G,boundary,hole,witness,area,source=source_surface();paths+=source
        locked={str(p):sha(p) for p in paths}
        R,T,native=execute(out/'plc',P,boundary,hole,witness,120,paths)
        assert np.array_equal(R[:len(P)],P)
        q=quality(R[T]);quality_ok=bool(np.min(q)>=.1)
        # 全部原面都具有唯一体单元邻接关系；这是位置审计，不是外缘认证。
        faces=T[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        keys,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
        rows=np.flatnonzero(count[inv]==1)
        owners={tuple(sorted(faces[i])):int(i//4) for i in rows}
        wall_owners=np.array([owners[tuple(sorted(f))] for f in F])
        tri=R[F];center=tri.mean(1);vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        normal=vec/np.linalg.norm(vec,axis=1)[:,None]
        cell_center=R[T[wall_owners]].mean(1)
        distance=np.einsum('ij,ij->i',cell_center-center,normal)
        assert np.all(distance>0)
        np.savez_compressed(out/'mesh.npz',vertices=R,tetrahedra=T,wall=F,farfield=G,
                            originalFaceIDs=np.arange(1,len(F)+1),wallOwners=wall_owners)
        write_su2(out/'mesh.su2',R,T,F,G)
        surface=dict(vertices=R[:9878].tolist(),facesOneBased=(F+1).tolist(),faceIDs=list(range(1,len(F)+1)),
                     identity='original rocket wall; X axis, nose +X; SI; outward from solid')
        (out/'surface.json').write_text(json.dumps(surface,separators=(',',':'))+'\n')
        np.savetxt(out/'wall_sampling_geometry.csv',np.c_[np.arange(1,len(F)+1),wall_owners,area,distance,center,cell_center],
            delimiter=',',comments='',header='original_face_id,owner_zero_based,area_m2,normal_distance_m,wall_x,wall_y,wall_z,cell_x,cell_y,cell_z',fmt='%.17g')
        c=len(T)+len(boundary);payload=c*8*5*8
        result.update(passed=quality_ok,geometryPassed=True,qualityScreenPassed=quality_ok,
            originalWallExactlyPreserved=True,originalFarfieldExactlyPreserved=True,
            originalAreaM2=float(area.sum()),originalVertices=9878,originalFaces=len(F),vertices=len(R),cells=len(T),
            faces=len(keys),ghostCells=len(boundary),farfieldFaces=len(G),
            cellObjectsIncludingGhost=c,conservedVariables=5,retainedTimeLevels=4,stageArrayPayloadBytes=payload,
            minimumSICN=float(q.min()),tetsBelowSICN01=int(np.sum(q<.1)),
            qualityPercentiles=np.percentile(q,[0,1,50,99,100]).tolist(),
            wallNormalDistancePercentilesM=np.percentile(distance,[0,1,50,99,100]).tolist(),
            outerBoundsM=[P[9878:].min(0).tolist(),P[9878:].max(0).tolist()],
            native=native,sourceSHA256=locked,
            filesSHA256={n:sha(out/n) for n in ['mesh.npz','mesh.su2','surface.json','wall_sampling_geometry.csv']},
            capacityStatus='exact topology and stage payload only; full peak unknown',
            pending=['native Euler preparation and exact peak feasibility','spatial/outer-domain accuracy',
                     'wall-adjacent to attached-edge evidence','CEA/NASA7 compatibility','public thermal validation'])
        assert quality_ok,'SICN minimum below original 0.1 geometry screen'
    except Exception as exc:
        result['error']=repr(exc);(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        result['sourceStableAfter']={str(p):sha(p)==locked[str(p)] for p in paths} if 'locked' in locals() else {}
        result['elapsedSeconds']=time.monotonic()-start
        result['outputBytesBeforeReport']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        if result['outputBytesBeforeReport']>=80*1024**2:result.update(passed=False,error='output exceeds 80MiB')
        if not all(result['sourceStableAfter'].values()):result.update(passed=False,error='source changed')
        (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','cells','vertices','faces','ghostCells','minimumSICN','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    raise SystemExit(0 if prepare(p.parse_args().out)['passed'] else 1)
