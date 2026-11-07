# -*- coding: utf-8 -*-
"""一次Gmsh/Netgen完整体网格优化；严格保留原面，原生结果先保存后审核。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
import traceback
import numpy as np
from readEilmerNative import sha
from prepareHybridRocketEulerMesh import ROOT,write_su2
from reviewSU2ExplicitQuality import quality
from tetgenFixedInterface import tet_audit,sorted_cells
from cgalTetrahedralContract import faces_and_neighbors
from scipy.sparse.csgraph import connected_components

RUNTIME=ROOT/'tools/runtime/gmsh/4.15.2'


def appended_cohort(node_tags,xyz,element_tags,connections,P,T,first):
    """完全离散体的旧单元未被Gmsh删除：按来源编号提取全部新批次，不按质量删单元。"""
    node_tags=np.asarray(node_tags);element_tags=np.asarray(element_tags)
    connections=np.asarray(connections).reshape(-1,4);end=first+len(T)
    assert len(np.unique(node_tags))==len(node_tags) and len(np.unique(element_tags))==len(element_tags)
    old=(element_tags>=first)&(element_tags<end)
    assert np.array_equal(element_tags[old],np.arange(first,end)),'old element cohort incomplete'
    assert np.array_equal(np.sort(connections[old]-1,axis=1),np.sort(T,axis=1)),'old connectivity changed'
    lookup={int(n):i for i,n in enumerate(node_tags)}
    assert np.array_equal(np.array([xyz[lookup[i+1]] for i in range(len(P))]),P),'old coordinates changed'
    assert np.any(~old) and np.array_equal(element_tags[~old],np.arange(end,end+sum(~old))),'new cohort tags invalid'
    return element_tags[~old],connections[~old]


def optimize(source,out):
    """模块一：固定输入、运行库及单次调用；不重建原PLC或覆盖失败资产。"""
    start=time.monotonic();source,out=map(Path,(source,out));out=out.resolve()
    out.relative_to(ROOT/'outputs');out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_optimize.py')
    paths=[source/n for n in ['report.json','mesh.npz','surface.json']]+[
        RUNTIME/'runtime.json',RUNTIME/'gmsh.py',RUNTIME/'gmsh-4.15.dll',Path(__file__).resolve(),
        ROOT/'docs/validation/v194/PLAN.md',ROOT/'scripts/prepareHybridRocketEulerMesh.py',
        ROOT/'scripts/reviewSU2ExplicitQuality.py',ROOT/'scripts/tetgenFixedInterface.py',
        ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'scripts/auditSU2HybridTopology.py']
    locked={str(p):sha(p) for p in paths};result=dict(schema='RAH.HybridRocketEulerMesh.v1',passed=False,
        executionAuthorized=False,physicalQualification=False,newSolvedHeatAreaM2=0,sourceSHA256=locked)
    initialized=False
    try:
        rt=json.loads((RUNTIME/'runtime.json').read_text())
        for row in rt['files']:assert sha(RUNTIME/row['path'])==row['sha256']
        old=json.loads((source/'report.json').read_text());assert old['geometryPassed'] and old['tetsBelowSICN01']==6
        for n,d in old['filesSHA256'].items():assert sha(source/n)==d
        mesh=np.load(source/'mesh.npz');P,T,F,G=[mesh[k] for k in ['vertices','tetrahedra','wall','farfield']]
        boundary=np.r_[F[:,[0,2,1]],G];fixed=np.unique(boundary)
        assert np.array_equal(fixed,np.arange(13166))
        inner=np.unique(F);far=np.unique(G);interior=np.setdiff1d(np.arange(len(P)),fixed)
        sys.path.insert(0,str(RUNTIME));import gmsh
        assert gmsh.__version__=='4.15.2'
        gmsh.initialize();initialized=True;gmsh.option.setNumber('General.NumThreads',2)
        gmsh.option.setNumber('Mesh.Renumber',0);gmsh.model.add('rocket-original-wall-euler')
        gmsh.model.addDiscreteEntity(2,1);gmsh.model.addDiscreteEntity(2,2);gmsh.model.addDiscreteEntity(3,1,[1,2])
        for dim,tag,ids in [(2,1,inner),(2,2,far),(3,1,interior)]:gmsh.model.mesh.addNodes(dim,tag,ids+1,P[ids].ravel())
        gmsh.model.mesh.addElementsByType(1,2,np.arange(1,len(F)+1),(F[:,[0,2,1]]+1).ravel())
        gmsh.model.mesh.addElementsByType(2,2,np.arange(len(F)+1,len(F)+len(G)+1),(G+1).ravel())
        tags=np.arange(len(F)+len(G)+1,len(F)+len(G)+len(T)+1)
        gmsh.model.mesh.addElementsByType(1,4,tags,(T+1).ravel())
        before=gmsh.model.mesh.getElementQualities(tags,'minSICN');local=quality(P[T])
        import_error=float(np.max(abs(before-local)));assert import_error<1e-10 and np.sum(before<.1)==6
        gmsh.option.setNumber('Mesh.OptimizeThreshold',.1)
        tick=time.monotonic();gmsh.model.mesh.optimize('Netgen',force=True,niter=1,dimTags=[(3,1)])
        result['nativeSeconds']=time.monotonic()-tick
        node_tags,xyz,_=gmsh.model.mesh.getNodes();xyz=xyz.reshape(-1,3)
        kinds,element_tags,connections=gmsh.model.mesh.getElements(3,1);assert list(kinds)==[4]
        # 模块二：任何后续审核之前先保存原生数组，查询异常不导致原生重算。
        np.savez_compressed(out/'optimized_raw.npz',nodeTags=node_tags,vertices=xyz,
                            tetrahedronTags=element_tags[0],tetrahedronNodeTags=connections[0].reshape(-1,4))
        # 固定4.15.2的deMeshGRegion对完全离散体直接返回；保留原始数组后才选择新批次。
        new_tags,new_connections=appended_cohort(node_tags,xyz,element_tags[0],connections[0],P,T,len(F)+len(G)+1)
        result.update(oldCohortExactlyMatched=True,removedOldInputCohortCells=len(T),
                      newCohortAllCellsPreserved=True,rawCells=len(element_tags[0]))
        element_tags=[new_tags];connections=[new_connections.ravel()]
        gmsh.model.mesh.rebuildNodeCache(False);gmsh.model.mesh.rebuildElementCache(False)
        order={int(n):i for i,n in enumerate(node_tags)}
        assert all(int(i+1) in order for i in fixed)
        assert np.array_equal(np.array([xyz[order[int(i+1)]] for i in fixed]),P[fixed]),'boundary node moved'
        for tag,tri in [(1,F),(2,G)]:
            kinds,_,conn=gmsh.model.mesh.getElements(2,tag);assert list(kinds)==[2]
            assert np.array_equal(sorted_cells(conn[0].reshape(-1,3)-1),sorted_cells(tri)),'declared surface triangles changed'
        free=np.flatnonzero(~np.isin(node_tags,fixed+1));R=np.r_[P[fixed],xyz[free]]
        mapping={int(i+1):int(i) for i in fixed}
        mapping.update({int(node_tags[i]):int(len(fixed)+j) for j,i in enumerate(free)})
        U=np.array([[mapping[int(v)] for v in row] for row in connections[0].reshape(-1,4)],dtype=np.int64)
        used=np.unique(U);assert np.array_equal(used[:len(fixed)],fixed)
        compact=np.full(len(R),-1,int);compact[used]=np.arange(len(used));R=R[used];U=compact[U]
        actual,geometry=tet_audit(R,U)
        assert np.array_equal(sorted_cells(actual),sorted_cells(boundary)),'actual volume boundary differs from original triangles'
        _,graph=faces_and_neighbors(U);assert connected_components(graph,directed=False,return_labels=False)==1
        assert abs(geometry['volume']/old['native']['expectedFluidVolume']-1)<2e-10
        q=quality(R[U]);engine=gmsh.model.mesh.getElementQualities(element_tags[0],'minSICN')
        quality_error=float(np.max(abs(q-engine)));assert quality_error<1e-10
        faces=U[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        keys,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
        rows=np.flatnonzero(count[inv]==1);owners={tuple(sorted(faces[i])):int(i//4) for i in rows}
        wall_owner=np.array([owners[tuple(sorted(f))] for f in F])
        tri=R[F];centers=tri.mean(1);vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(vec,axis=1)/2
        sample=R[U[wall_owner]].mean(1);distance=np.einsum('ij,ij->i',sample-centers,vec/(2*area[:,None]))
        assert np.all(distance>0)
        np.savez_compressed(out/'mesh.npz',vertices=R,tetrahedra=U,wall=F,farfield=G,wallOwners=wall_owner,originalFaceIDs=mesh['originalFaceIDs'])
        write_su2(out/'mesh.su2',R,U,F,G);shutil.copyfile(source/'surface.json',out/'surface.json')
        np.savetxt(out/'wall_sampling_geometry.csv',np.c_[np.arange(1,len(F)+1),wall_owner,area,distance,centers,sample],
                   delimiter=',',comments='',header='original_face_id,owner_zero_based,area_m2,normal_distance_m,wall_x,wall_y,wall_z,cell_x,cell_y,cell_z',fmt='%.17g')
        good=bool(q.min()>=.1);size_ok=len(U)<=500000 and len(R)<=100000;c=len(U)+len(boundary)
        result.update(passed=good and size_ok,geometryPassed=True,qualityScreenPassed=good,sizeScreenPassed=size_ok,
            originalWallExactlyPreserved=True,originalFarfieldExactlyPreserved=True,originalVertices=9878,originalFaces=len(F),
            originalAreaM2=float(area.sum()),vertices=len(R),cells=len(U),faces=len(keys),ghostCells=len(boundary),farfieldFaces=len(G),
            cellObjectsIncludingGhost=c,conservedVariables=5,retainedTimeLevels=4,stageArrayPayloadBytes=c*8*5*8,
            minimumSICN=float(q.min()),tetsBelowSICN01=int(np.sum(q<.1)),qualityPercentiles=np.percentile(q,[0,1,50,99,100]).tolist(),
            wallNormalDistancePercentilesM=np.percentile(distance,[0,1,50,99,100]).tolist(),outerBoundsM=old['outerBoundsM'],
            gmshImportMaximumQualityDifference=import_error,gmshFinalMaximumQualityDifference=quality_error,
            nativeGeometry=geometry,sourceSHA256=locked,algorithm='Gmsh4.15.2/Netgen one volume optimization invocation',
            filesSHA256={n:sha(out/n) for n in ['mesh.npz','mesh.su2','surface.json','wall_sampling_geometry.csv']})
        assert good and size_ok,'Netgen candidate rejected by original quality or stated resource screen'
    except Exception as error:
        result['error']=repr(error);(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8')
    finally:
        if initialized:gmsh.finalize()
        result['sourceStableAfter']={str(p):sha(p)==locked[str(p)] for p in paths}
        if not all(result['sourceStableAfter'].values()):result.update(passed=False,error='input changed')
        result['elapsedSeconds']=time.monotonic()-start
        result['outputBytesBeforeReport']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        if result['outputBytesBeforeReport']>80*1024**2:result.update(passed=False,error='output over 80MiB')
        (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ['passed','error','cells','vertices','minimumSICN','tetsBelowSICN01','nativeSeconds','elapsedSeconds']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();raise SystemExit(0 if optimize(a.source,a.out)['passed'] else 1)
