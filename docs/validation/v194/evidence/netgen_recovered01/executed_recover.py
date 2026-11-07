# -*- coding: utf-8 -*-
"""按准确Gmsh离散体追加语义恢复已存Netgen新批次；不重跑优化或删除低质单元。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import traceback
import numpy as np
from readEilmerNative import sha
from prepareHybridRocketEulerMesh import ROOT,write_su2
from tetgenFixedInterface import tet_audit,sorted_cells
from cgalTetrahedralContract import faces_and_neighbors
from scipy.sparse.csgraph import connected_components
from reviewSU2ExplicitQuality import quality


def recover(source,attempt,out):
    start=time.monotonic();source,attempt,out=map(Path,(source,attempt,out));out=out.resolve()
    out.relative_to(ROOT/'outputs');out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_recover.py')
    paths=[source/n for n in ['mesh.npz','report.json','surface.json']]+[
        attempt/'optimized_raw.npz',attempt/'report.json',attempt/'executed_optimize.py',
        ROOT/'outputs/v194_validation/sources/deMeshGRegion_excerpt.txt',
        ROOT/'outputs/v194_validation/sources/meshGRegionNetgen.cpp',
        ROOT/'docs/validation/v194/RECOVERY_PLAN.md',Path(__file__).resolve()]
    locked={str(p):sha(p) for p in paths}
    result=dict(schema='RAH.HybridRocketEulerMesh.v1',passed=False,executionAuthorized=False,
                physicalQualification=False,newSolvedHeatAreaM2=0)
    try:
        old=json.loads((source/'report.json').read_text());failed=json.loads((attempt/'report.json').read_text())
        assert not failed['passed'] and failed['error']=="AssertionError('duplicate tetrahedra')"
        assert all(failed['sourceStableAfter'].values())
        original_source=ROOT/'scripts/optimizeHybridRocketEulerMesh.py'
        assert sha(attempt/'executed_optimize.py')==failed['sourceSHA256'][str(original_source)]
        for n,d in old['filesSHA256'].items():assert sha(source/n)==d
        mesh=np.load(source/'mesh.npz');P,T,F,G=[mesh[k] for k in ['vertices','tetrahedra','wall','farfield']]
        raw=np.load(attempt/'optimized_raw.npz');tags=raw['tetrahedronTags'].astype(np.int64)
        nodes=raw['nodeTags'].astype(np.int64);xyz=raw['vertices'];tet_tags=raw['tetrahedronNodeTags'].astype(np.int64)
        assert len(np.unique(tags))==len(tags) and len(np.unique(nodes))==len(nodes)
        first=len(F)+len(G)+1;end=first+len(T)
        is_old=(tags>=first)&(tags<end);assert np.array_equal(tags[is_old],np.arange(first,end))
        # 允许原生内部把正四面体方向反转导入Netgen，但节点集合/原行身份必须不变。
        assert np.array_equal(np.sort(tet_tags[is_old]-1,axis=1),np.sort(T,axis=1))
        node_lookup={int(n):i for i,n in enumerate(nodes)}
        assert np.array_equal(np.array([xyz[node_lookup[i+1]] for i in range(len(P))]),P)
        assert np.array_equal(tags[~is_old],np.arange(end,end+sum(~is_old)))
        new_tags=tet_tags[~is_old];used=np.unique(new_tags)
        fixed=np.unique(np.r_[F,G]);assert np.array_equal(fixed,np.arange(13166))
        assert np.array_equal(used[:len(fixed)],fixed+1)
        R=np.array([xyz[node_lookup[int(i)]] for i in used]);lookup={int(n):i for i,n in enumerate(used)}
        U=np.array([[lookup[int(n)] for n in row] for row in new_tags],dtype=np.int64)
        assert np.array_equal(R[:len(fixed)],P[fixed])
        np.savez_compressed(out/'new_cohort.npz',vertices=R,tetrahedra=U,wall=F,farfield=G,
                            sourceElementTags=tags[~is_old],sourceNodeTags=used)
        boundary=np.r_[F[:,[0,2,1]],G];actual,ga=tet_audit(R,U)
        assert np.array_equal(sorted_cells(actual),sorted_cells(boundary))
        _,graph=faces_and_neighbors(U);assert connected_components(graph,directed=False,return_labels=False)==1
        assert abs(ga['volume']/old['native']['expectedFluidVolume']-1)<2e-10
        q=quality(R[U]);good=bool(q.min()>=.1);size_ok=len(U)<=500000 and len(R)<=100000
        faces=U[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        keys,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
        rows=np.flatnonzero(count[inv]==1);owners={tuple(sorted(faces[i])):int(i//4) for i in rows}
        wall_owner=np.array([owners[tuple(sorted(f))] for f in F]);tri=R[F];centers=tri.mean(1)
        vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(vec,axis=1)/2
        sample=R[U[wall_owner]].mean(1);distance=np.einsum('ij,ij->i',sample-centers,vec/(2*area[:,None]));assert np.all(distance>0)
        np.savez_compressed(out/'mesh.npz',vertices=R,tetrahedra=U,wall=F,farfield=G,wallOwners=wall_owner,originalFaceIDs=mesh['originalFaceIDs'])
        write_su2(out/'mesh.su2',R,U,F,G);shutil.copyfile(source/'surface.json',out/'surface.json')
        np.savetxt(out/'wall_sampling_geometry.csv',np.c_[np.arange(1,len(F)+1),wall_owner,area,distance,centers,sample],
            delimiter=',',comments='',header='original_face_id,owner_zero_based,area_m2,normal_distance_m,wall_x,wall_y,wall_z,cell_x,cell_y,cell_z',fmt='%.17g')
        c=len(U)+len(boundary)
        result.update(passed=good and size_ok,geometryPassed=True,qualityScreenPassed=good,sizeScreenPassed=size_ok,
            originalWallExactlyPreserved=True,originalFarfieldExactlyPreserved=True,originalAreaM2=float(area.sum()),
            originalVertices=9878,originalFaces=len(F),vertices=len(R),cells=len(U),faces=len(keys),ghostCells=len(boundary),farfieldFaces=len(G),
            cellObjectsIncludingGhost=c,conservedVariables=5,retainedTimeLevels=4,stageArrayPayloadBytes=c*8*5*8,
            minimumSICN=float(q.min()),tetsBelowSICN01=int(np.sum(q<.1)),qualityPercentiles=np.percentile(q,[0,1,50,99,100]).tolist(),
            wallNormalDistancePercentilesM=np.percentile(distance,[0,1,50,99,100]).tolist(),outerBoundsM=old['outerBoundsM'],
            nativeGeometry=ga,nativeSeconds=failed['nativeSeconds'],originalFailureRetained=True,
            oldCohortExactlyMatched=True,removedOldInputCohortCells=int(sum(is_old)),newCohortAllCellsPreserved=True,
            rawCells=len(tags),rawVertices=len(nodes),recovery='source-proven appended cohort; no rerun, no quality-based deletion',
            filesSHA256={n:sha(out/n) for n in ['mesh.npz','mesh.su2','surface.json','wall_sampling_geometry.csv']},sourceSHA256=locked)
        assert good and size_ok,'new Netgen cohort fails original quality or resource screen'
    except Exception as error:
        result['error']=repr(error);(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        result['sourceStableAfter']={str(p):sha(p)==locked[str(p)] for p in paths}
        if not all(result['sourceStableAfter'].values()):result.update(passed=False,error='input changed')
        result['elapsedSeconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','cells','vertices','minimumSICN','tetsBelowSICN01','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['source','attempt','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();raise SystemExit(0 if recover(a.source,a.attempt,a.out)['passed'] else 1)
