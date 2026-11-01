"""一次固定PLC真实完整外域：保留原棱柱和物理面，只重建纯四面体域。"""
from pathlib import Path
import json,time,traceback
import numpy as np
from scipy.io import loadmat,savemat
from prepareTetGenRuntime import ROOT,sha
from tetgenPLCInterface import execute,winding
from cgalTetrahedralContract import faces_and_neighbors
from tetgenFixedInterface import sorted_cells
from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from reviewSU2ExplicitQuality import quality

def write(p,r):p.write_text(json.dumps(r,indent=2)+'\n')

def main():
    start=time.monotonic();source=ROOT/'outputs/v141_validation/actual';folder=ROOT/'outputs/v143_validation/actual'
    small=ROOT/'docs/validation/v143/hollow_test/report.json';assert json.loads(small.read_text())['passed']
    folder.mkdir(parents=True,exist_ok=False)
    paths=[Path(__file__).resolve(),small,ROOT/'docs/validation/v143/plan.md']
    paths += [ROOT/'scripts'/n for n in ['tetgenPLCInterface.py','tetgenFixedInterface.py','cgalTetrahedralContract.py',
        'auditSU2HybridTopology.py','meshSU2HybridVolume.py','su2PyramidGeometry.py','reviewSU2ExplicitQuality.py']]
    paths += [source/n for n in ['hybrid_volume.mat','mesh_process.json','candidate_report.json','prepared_envelope.mat',
        'boundary_audit/audit_input.mat','boundary_audit/process.json','boundary_audit/intersections.json']]
    locked={str(p):sha(p) for p in paths};code=1;error=None
    try:
        oldreceipt=json.loads((source/'mesh_process.json').read_text());assert oldreceipt['sourceStable'] and oldreceipt['exitCode']==1
        assert sha(source/'hybrid_volume.mat')==oldreceipt['outputSHA256']['hybrid_volume.mat']
        cg=json.loads((source/'boundary_audit/process.json').read_text());assert cg['exitCode']==0 and cg['runtimeStableAfter']
        for n,h in cg['outputSHA256'].items():assert sha(source/'boundary_audit'/n)==h
        old=loadmat(str(source/'hybrid_volume.mat'));P=old['vertices'];T=old['tetrahedra'].astype(np.int64)-1
        W=old['wedges'].astype(np.int64)-1;far=old['farfield'].astype(np.int64)-1;V=old['sourceVertices'];F=old['sourceFaces'].astype(np.int64)-1
        boundary,_=faces_and_neighbors(T);ids=np.unique(boundary);mapping=np.full(len(P),-1,np.int64);mapping[ids]=np.arange(len(ids))
        assert len(ids)==13166 and len(boundary)==26324 and len(W)==634916
        # 模块1：将已通过CGAL的包络与本次实际四面体内边界逐点逐面对应。
        prepared=loadmat(str(source/'prepared_envelope.mat'));top=prepared['farfield'].astype(np.int64)-1
        pids=old['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1
        inverse=np.full(len(prepared['vertices']),-1,np.int64);valid=pids<len(inverse);inverse[pids[valid]]=np.flatnonzero(valid)
        assert np.all(inverse[top]>=0)
        assert np.array_equal(sorted_cells(boundary),sorted_cells(np.r_[inverse[top[:,[0,2,1]]],far]))
        cgInput=loadmat(str(source/'boundary_audit/audit_input.mat'));nodes,back=np.unique(top,return_inverse=True)
        assert np.array_equal(cgInput['V'],prepared['vertices'][nodes]) and np.array_equal(cgInput['F']-1,back.reshape(top.shape))
        assert np.array_equal(P[inverse[top]],prepared['vertices'][top])
        tri=V[F];vv=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2]))/6
        hole=np.sum(tri.sum(axis=1)/4*vv[:,None],axis=0)/sum(vv);assert abs(winding(V,F,hole)-1)<1e-8
        witness=P[T[0]].mean(axis=0)
        Q,U,inter=execute(folder/'native',P[ids],mapping[boundary],hole,witness,180,paths)
        # 模块2：边界按显式原编号恢复；新内部点单独标0，旧内部点不做最近邻迁移。
        newCount=len(Q)-len(ids);newIDs=np.arange(len(P),len(P)+newCount,dtype=np.int64)
        localToAll=np.r_[ids,newIDs];allP=np.r_[P,Q[len(ids):]];allT=localToAll[U]
        used=np.unique(np.r_[W.ravel(),allT.ravel(),far.ravel()]);remap=np.full(len(allP),-1,np.int64);remap[used]=np.arange(len(used))
        assert np.array_equal(used[:len(V)],np.arange(len(V)))
        oldMap=np.where(used<len(P),used+1,0);preparedMap=np.zeros(len(used),np.int64);wasOld=used<len(P);preparedMap[wasOld]=pids[used[wasOld]]+1
        data=dict(vertices=allP[used],wedges=remap[W]+1,tetrahedra=remap[allT]+1,farfield=remap[far]+1,
            sourceVertices=V,sourceFaces=F+1,normalHeightsByColumnM=old['normalHeightsByColumnM'],
            sourceV141NodeIDsOneBased=oldMap,sourcePreparedNodeIDsOneBased=preparedMap,
            plcToHybridNodeIDsOneBased=remap[localToAll]+1,sourcePLCBoundaryNodeIDsOneBased=ids+1)
        savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True)
        assert np.array_equal(data['vertices'][data['wedges']-1],P[W]) and np.array_equal(data['vertices'][data['farfield']-1],P[far])
        ga=audit(data);wj=float(min(wedge_jacobian_min(data['vertices'][data['wedges']-1])));assert wj>0
        write(folder/'geometry_audit.json',dict(ga,minimumPrismJacobian=wj))
        val=quality(data['vertices'][data['tetrahedra']-1]);bad=np.flatnonzero(val<.1)
        area=float(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2)
        report=dict(inter,nodes=len(used),prisms=len(W),tetrahedra=len(allT),totalElements=len(W)+len(allT),
            minimumTetSICN=float(min(val)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=val[bad].tolist(),
            badTetCentersM=data['vertices'][data['tetrahedra'][bad]-1].mean(axis=1).tolist(),originalAreaM2=area,
            originalWallExactlyPreserved=True,allPrismsAndFarfieldExactlyPreserved=True,priorCGALSurfaceExactlyMatched=True,
            priorBadTets=31,qualityGatePassed=bool(min(val)>=.1),sizeGatePassed=bool(len(used)<=400000 and len(W)+len(allT)<=900000),
            geometryPassed=True,exteriorAlgorithm='TetGen1.5.1 fixed-boundary PLC with explicit body hole',
            sourceCandidateSHA256=sha(source/'hybrid_volume.mat'),candidateAccepted=False,physicsQualified=False,newQualifiedArea=0)
        report['candidateAccepted']=report['qualityGatePassed'] and report['sizeGatePassed']
        write(folder/'candidate_report.json',report)
        assert report['candidateAccepted'],'quality or size gate failed'
        code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        r=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,error=error,sourceSHA256=locked,
            sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},
            outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},physicsQualified=False,newQualifiedArea=0)
        r['sourceStable']=all(r['sourceStableAfter'].values());write(folder/'mesh_process.json',r)
        print(json.dumps({k:v for k,v in r.items() if k not in ['sourceSHA256','sourceStableAfter','outputSHA256']},indent=2))
        assert r['sourceStable']
    raise SystemExit(code)

if __name__=='__main__':main()
