# -*- coding: utf-8 -*-
"""一次成熟TetGen固定边界refine；失败PLC原样保留，不扫描参数。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import traceback
import numpy as np
from tetgenFixedInterface import execute
from reviewSU2ExplicitQuality import quality
from prepareHybridRocketEulerMesh import ROOT,write_su2
from readEilmerNative import sha


def repair(source,review,out):
    """模块一：明确原失败和可动内点依据；只在新目录作一次修复。"""
    start=time.monotonic();source,review,out=map(Path,(source,review,out));out=out.resolve()
    out.relative_to(ROOT/'outputs');out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_repair.py')
    paths=[source/'report.json',source/'mesh.npz',source/'surface.json',review,Path(__file__).resolve(),
           ROOT/'scripts/prepareHybridRocketEulerMesh.py',ROOT/'docs/validation/v193/REPAIR_PLAN.md']
    locked={str(p):sha(p) for p in paths}
    result=dict(passed=False,executionAuthorized=False,physicalQualification=False,newSolvedHeatAreaM2=0)
    try:
        old=json.loads((source/'report.json').read_text());check=json.loads(review.read_text())
        assert old['geometryPassed'] and not old['qualityScreenPassed'] and check['passed']
        assert len(check['badCells'])==6 and all(r['fixedVertexCount']<4 for r in check['badCells'])
        for n,d in old['filesSHA256'].items():assert sha(source/n)==d
        m=np.load(source/'mesh.npz');P=m['vertices'];T=m['tetrahedra'];F=m['wall'];G=m['farfield']
        Q,U,native=execute(out/'native',P,T,60,paths)
        # 接口保留旧自由点索引，压缩只移除不被新体网格使用的内部点；原边界不变。
        used=np.unique(U);assert np.array_equal(used[:13166],np.arange(13166))
        remap=np.full(len(Q),-1,dtype=int);remap[used]=np.arange(len(used));R=Q[used];U=remap[U]
        F=remap[F];G=remap[G];assert np.array_equal(R[:13166],P[:13166])
        q=quality(R[U]);good=bool(q.min()>=.1)
        faces=U[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        keys,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
        rows=np.flatnonzero(count[inv]==1);owners={tuple(sorted(faces[i])):int(i//4) for i in rows}
        wall_owners=np.array([owners[tuple(sorted(f))] for f in F])
        tri=R[F];center=tri.mean(1);vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(vec,axis=1)/2
        cell_center=R[U[wall_owners]].mean(1);distance=np.einsum('ij,ij->i',cell_center-center,vec/(2*area[:,None]))
        assert np.all(distance>0)
        np.savez_compressed(out/'mesh.npz',vertices=R,tetrahedra=U,wall=F,farfield=G,originalFaceIDs=m['originalFaceIDs'],wallOwners=wall_owners)
        write_su2(out/'mesh.su2',R,U,F,G);shutil.copyfile(source/'surface.json',out/'surface.json')
        np.savetxt(out/'wall_sampling_geometry.csv',np.c_[np.arange(1,len(F)+1),wall_owners,area,distance,center,cell_center],delimiter=',',comments='',
            header='original_face_id,owner_zero_based,area_m2,normal_distance_m,wall_x,wall_y,wall_z,cell_x,cell_y,cell_z',fmt='%.17g')
        c=len(U)+len(F)+len(G)
        result.update(schema='RAH.HybridRocketEulerMesh.v1',passed=good,geometryPassed=True,qualityScreenPassed=good,
            originalWallExactlyPreserved=True,originalFarfieldExactlyPreserved=True,originalAreaM2=float(area.sum()),
            originalVertices=9878,originalFaces=len(F),vertices=len(R),cells=len(U),faces=len(keys),ghostCells=len(F)+len(G),farfieldFaces=len(G),
            cellObjectsIncludingGhost=c,conservedVariables=5,retainedTimeLevels=4,stageArrayPayloadBytes=c*8*5*8,
            minimumSICN=float(q.min()),tetsBelowSICN01=int(np.sum(q<.1)),qualityPercentiles=np.percentile(q,[0,1,50,99,100]).tolist(),
            wallNormalDistancePercentilesM=np.percentile(distance,[0,1,50,99,100]).tolist(),outerBoundsM=old['outerBoundsM'],
            native=native,sourceSHA256=locked,sourceFailureRetained=True,
            filesSHA256={n:sha(out/n) for n in ['mesh.npz','mesh.su2','surface.json','wall_sampling_geometry.csv']},
            capacityStatus=old['capacityStatus'],pending=old['pending'])
        assert good,'Fixed-boundary refine still below original SICN 0.1 screen; stop this candidate'
    except Exception as e:
        result['error']=repr(e);(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        result['sourceStableAfter']={str(p):sha(p)==locked[str(p)] for p in paths}
        result['elapsedSeconds']=time.monotonic()-start
        result['outputBytesBeforeReport']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        if not all(result['sourceStableAfter'].values()):result.update(passed=False,error='source changed')
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','cells','vertices','minimumSICN','tetsBelowSICN01','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['source','review','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();raise SystemExit(0 if repair(a.source,a.review,a.out)['passed'] else 1)
