"""对v136原始完整外域执行唯一MMG局部候选，冻结混合接口并重组原W/Q/T。"""
from pathlib import Path
import json,time,traceback
import numpy as np
from scipy.io import loadmat,savemat
from prepareCGALRemeshingRuntime import ROOT,sha
from mmgLocalInterface import execute
from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from reviewSU2ExplicitQuality import quality

def write(path,data):path.write_text(json.dumps(data,indent=2)+'\n')

def main():
    start=time.monotonic();source=ROOT/'outputs/v136_validation/explicit';folder=ROOT/'outputs/v138_validation/actual'
    selection=ROOT/'docs/validation/v137/real_admissibility/report.json'
    case=ROOT/'docs/validation/v138/hollow_test02/report.json';test=json.loads(case.read_text());assert test['passed']
    assert not folder.exists();folder.mkdir(parents=True)
    paths=[source/'hybrid_volume.mat',source/'mesh_report.json',source/'mesh_process.json',selection,case,
        Path(__file__).resolve(),ROOT/'scripts/mmgLocalInterface.py',ROOT/'scripts/auditSU2HybridTopology.py',
        ROOT/'scripts/meshSU2HybridVolume.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',ROOT/'scripts/su2PyramidGeometry.py']
    locked={p.resolve().as_posix():sha(p) for p in paths};code=1;error=None
    try:
        oldprocess=json.loads((source/'mesh_process.json').read_text());assert oldprocess['exitCode']==1 and oldprocess['sourceStable']
        assert sha(source/'hybrid_volume.mat')==oldprocess['outputSHA256']['hybrid_volume.mat']
        old=loadmat(str(source/'hybrid_volume.mat'));P=old['vertices'];T=old['tetrahedra'].astype(int)-1
        W=old['wedges'].astype(int)-1;Q=old['pyramids'].astype(int)-1;far=old['farfield'].astype(int)-1;V=old['sourceVertices'];F=old['sourceFaces']
        selected=np.zeros(len(T),bool);sel=json.loads(selection.read_text())
        for comp in sel['components']:selected[np.array(comp['sourceTetRowsOneBased'])-1]=True
        assert sum(selected)==888 and len(T)==171448
        ids=np.unique(T);mapin=np.full(len(P),-1,int);mapin[ids]=np.arange(len(ids));TT=mapin[T]
        PP,UU,inter=execute(folder/'native',P[ids],TT,selected,180,paths)
        # 原W/Q节点不交给MMG改变。用显式纯四面体节点编号恢复全局身份。
        mapping=np.r_[ids,np.arange(len(P),len(P)+len(PP)-len(ids))];allP=np.r_[P,PP[len(ids):]];allT=mapping[UU]
        used=np.unique(np.r_[W.ravel(),Q.ravel(),allT.ravel(),far.ravel()]);remap=np.full(len(allP),-1,int);remap[used]=np.arange(len(used))
        assert np.array_equal(used[:len(V)],np.arange(len(V)))
        data=dict(vertices=allP[used],wedges=remap[W]+1,pyramids=remap[Q]+1,tetrahedra=remap[allT]+1,
            farfield=remap[far]+1,sourceVertices=V,sourceFaces=F)
        savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True)
        assert np.array_equal(data['vertices'][data['wedges']-1],P[W])
        assert np.array_equal(data['vertices'][data['pyramids']-1],P[Q])
        assert np.array_equal(data['vertices'][data['farfield']-1],P[far])
        geom=audit(data);wj=float(min(wedge_jacobian_min(data['vertices'][data['wedges']-1])));assert wj>0
        write(folder/'geometry_audit.json',dict(geom,minimumPrismJacobian=wj))
        val=quality(data['vertices'][data['tetrahedra']-1]);oldq=quality(P[T]);bad=np.flatnonzero(val<.1)
        report=dict(inter,nodes=len(used),prisms=len(W),pyramids=len(Q),tetrahedra=len(allT),totalElements=len(W)+len(Q)+len(allT),
            minimumTetSICN=float(min(val)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),
            badTetSICN=val[bad].tolist(),badTetCentersM=data['vertices'][data['tetrahedra'][bad]-1].mean(axis=1).tolist(),
            previousBadTets=int(sum(oldq<.1)),previousMinimumSICN=float(min(oldq)),
            qualityGatePassed=bool(min(val)>=.1),sizeGatePassed=bool(len(used)<=400000 and len(W)+len(Q)+len(allT)<=900000),
            originalWallExactlyPreserved=True,retainedPrismsExactlyPreserved=True,retainedPyramidsExactlyPreserved=True,
            farfieldExactlyPreserved=True,geometryPassed=True,exteriorAlgorithm='MMG3D 5.8.0 local optim (v136 Delaunay input)',
            algorithm='MMG3D 5.8.0 optim/nosurf + required nonselected tetrahedra',
            physicsQualified=False,newQualifiedArea=0)
        write(folder/'candidate_report.json',report)
        oldreport=json.loads((source/'mesh_report.json').read_text());oldreport.update(report)
        oldreport['sourceCandidate']=str(source);oldreport['sourceCandidateSHA256']=sha(source/'hybrid_volume.mat')
        write(folder/'mesh_report.json',oldreport)
        assert report['qualityGatePassed'] and report['sizeGatePassed'],'candidate fails quality or size gate'
        code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        receipt=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,error=error,sourceSHA256=locked,
            sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},
            outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},physicsQualified=False,newQualifiedArea=0)
        receipt['sourceStable']=all(receipt['sourceStableAfter'].values());receipt['inputStable']=receipt['sourceStable']
        native=folder/'native/process/process.json'
        receipt['timedOut']=bool(native.exists() and json.loads(native.read_text()).get('stopReason')=='time')
        write(folder/'mesh_process.json',receipt)
        print(json.dumps(receipt,indent=2));assert receipt['sourceStable']
    raise SystemExit(code)
if __name__=='__main__':main()
