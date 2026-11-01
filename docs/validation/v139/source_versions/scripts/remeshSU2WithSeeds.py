"""v139唯一真实候选：完整内部边星中点细分，固定边界MMG保点优化，全混合域审核。"""
from pathlib import Path
import json,time,traceback
import numpy as np
from scipy.io import loadmat,savemat
from prepareCGALRemeshingRuntime import ROOT,sha
from splitClosedTetEdges import split,volumes
from mmgLocalInterface import execute,sorted_cells
from reviewMMGLocalCandidate import key_cells
from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from reviewSU2ExplicitQuality import quality

def write(path,data):path.write_text(json.dumps(data,indent=2)+'\n')

def main():
    start=time.monotonic();source=ROOT/'outputs/v138_validation/actual';folder=ROOT/'outputs/v139_validation/actual'
    oldpath=ROOT/'outputs/v136_validation/explicit/hybrid_volume.mat'
    selection=ROOT/'docs/validation/v137/real_admissibility/report.json';edgespath=source/'interior_seed_review.json'
    case=ROOT/'docs/validation/v139/seed_test/report.json';assert json.loads(case.read_text())['passed']
    folder.mkdir(parents=True,exist_ok=False)
    paths=[source/'hybrid_volume.mat',source/'mesh_report.json',source/'mesh_process.json',oldpath,selection,edgespath,case,
        Path(__file__).resolve()]+[ROOT/'scripts'/n for n in ['splitClosedTetEdges.py','mmgLocalInterface.py',
        'reviewMMGLocalCandidate.py','auditSU2HybridTopology.py','meshSU2HybridVolume.py','reviewSU2ExplicitQuality.py',
        'su2PyramidGeometry.py','cgalTetrahedralContract.py']]
    locked={p.resolve().as_posix():sha(p) for p in paths};code=1;error=None
    try:
        previous=json.loads((source/'mesh_process.json').read_text());assert previous['exitCode']==1 and previous['sourceStable']
        assert sha(source/'hybrid_volume.mat')==previous['outputSHA256']['hybrid_volume.mat']
        old=loadmat(str(source/'hybrid_volume.mat'));P=old['vertices'];T=old['tetrahedra'].astype(np.int64)-1
        W=old['wedges'].astype(np.int64)-1;Q=old['pyramids'].astype(np.int64)-1;far=old['farfield'].astype(np.int64)-1
        base=loadmat(str(oldpath));A=base['tetrahedra'].astype(np.int64)-1;selected0=np.zeros(len(A),bool)
        for part in json.loads(selection.read_text())['components']:selected0[np.array(part['sourceTetRowsOneBased'])-1]=True
        selected=~np.isin(key_cells(T),key_cells(A[~selected0]));assert sum(~selected)==170560 and sum(selected)==885
        edges=[tuple(v-1 for v in e['vertexIDsOneBased']) for e in json.loads(edgespath.read_text())['edges']];assert len(edges)==247
        # 逐边使用当前完整邻域，保存原始父单元行；优化之后不伪造父单元身份。
        S,U,mask,parent,sr=split(P,T,selected,edges)
        assert len(S)-len(P)==247 and len(U)-len(T)<=6000
        assert np.allclose(np.bincount(parent,weights=volumes(S,U),minlength=len(T)),volumes(P,T),rtol=2e-10,atol=0)
        data={k:old[k] for k in ['sourceVertices','sourceFaces','wedges','pyramids','farfield']}
        data.update(vertices=S,tetrahedra=U+1)
        pre=audit(data);wj=float(min(wedge_jacobian_min(P[W])));assert wj>0
        preparation=time.monotonic()-start;assert preparation<=60,'preparation budget exceeded'
        write(folder/'seed_geometry_audit.json',dict(pre,minimumPrismJacobian=wj,preparationSeconds=preparation))
        write(folder/'split_report.json',sr)
        savemat(str(folder/'seeded_volume.mat'),dict(data,parentTetRowsOneBased=parent+1,selectedTets=mask),do_compression=True)
        # 原纯四面体域单独传给MMG；棱柱和金字塔从未送入原生修改器。
        ids=np.unique(U);mapin=np.full(len(S),-1,np.int64);mapin[ids]=np.arange(len(ids))
        PP,UU,inter=execute(folder/'native',S[ids],mapin[U],mask,180,paths,ROOT/'tools/runtime/mmg/5.8.0-seeds/runtime.json')
        mapping=np.r_[ids,np.arange(len(S),len(S)+len(PP)-len(ids))];allP=np.r_[S,PP[len(ids):]];allT=mapping[UU]
        used=np.unique(np.r_[W.ravel(),Q.ravel(),allT.ravel(),far.ravel()]);remap=np.full(len(allP),-1,np.int64);remap[used]=np.arange(len(used))
        assert np.array_equal(used[:len(P)],np.arange(len(P)))
        data.update(vertices=allP[used],tetrahedra=remap[allT]+1,wedges=remap[W]+1,pyramids=remap[Q]+1,farfield=remap[far]+1)
        assert np.array_equal(data['vertices'][:len(P)],P)
        for name in ['sourceVertices','sourceFaces','wedges','pyramids','farfield']:assert np.array_equal(data[name],old[name])
        savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True)
        geom=audit(data);write(folder/'geometry_audit.json',dict(geom,minimumPrismJacobian=wj))
        val=quality(data['vertices'][data['tetrahedra']-1]);oldq=quality(P[T]);bad=np.flatnonzero(val<.1)
        seedCountPassed=inter['freeOutputNodes']==247 and len(used)==len(P)+247
        r=dict(inter,nodes=len(used),prisms=len(W),pyramids=len(Q),tetrahedra=len(allT),totalElements=len(W)+len(Q)+len(allT),
            minimumTetSICN=float(min(val)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=val[bad].tolist(),
            badTetCentersM=data['vertices'][data['tetrahedra'][bad]-1].mean(axis=1).tolist(),
            previousBadTets=int(sum(oldq<.1)),previousMinimumSICN=float(min(oldq)),qualityGatePassed=bool(min(val)>=.1),
            sizeGatePassed=bool(len(used)<=400000 and len(W)+len(Q)+len(allT)<=900000),seedCountPassed=bool(seedCountPassed),
            insertedNodes=247,preparationSeconds=preparation,originalCoordinatesExactlyPreserved=True,
            originalWallExactlyPreserved=True,retainedPrismsExactlyPreserved=True,retainedPyramidsExactlyPreserved=True,
            farfieldExactlyPreserved=True,geometryPassed=True,parentLineageScope='subdivision only; MMG swaps have no invented parent identity',
            exteriorAlgorithm='MMG3D 5.8.0 noinsert optim/nosurf (v138 + 247 complete edge stars)',
            algorithm='closed internal-edge midpoint subdivision + fixed-boundary keep-seeds MMG optimization',physicsQualified=False,newQualifiedArea=0)
        write(folder/'candidate_report.json',r)
        mr=json.loads((source/'mesh_report.json').read_text());mr.update(r);mr.update(sourceCandidate=str(source),sourceCandidateSHA256=sha(source/'hybrid_volume.mat'))
        write(folder/'mesh_report.json',mr)
        assert r['qualityGatePassed'] and r['sizeGatePassed'] and seedCountPassed,'candidate fails quality, size or seed count gate'
        code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        receipt=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,error=error,sourceSHA256=locked,
            sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},
            outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},physicsQualified=False,newQualifiedArea=0)
        receipt['sourceStable']=all(receipt['sourceStableAfter'].values());receipt['inputStable']=receipt['sourceStable']
        native=folder/'native/process/process.json';receipt['timedOut']=bool(native.exists() and json.loads(native.read_text()).get('stopReason')=='time')
        write(folder/'mesh_process.json',receipt);print(json.dumps(receipt,indent=2));assert receipt['sourceStable']
    raise SystemExit(code)

if __name__=='__main__':main()
