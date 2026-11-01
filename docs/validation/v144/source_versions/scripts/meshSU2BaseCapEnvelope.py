"""底面内部列外层切向重分布：一次壳审核与固定PLC，保护原壁/近壁/完整尾缘扇。"""
from pathlib import Path
import json,time,traceback
import numpy as np
from scipy.io import loadmat,savemat
from prepareTetGenRuntime import ROOT,sha
from recenterBaseOuterColumn import recenter
from buildSU2CreaseFan import interface_edges
from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from runCGALSurfaceAudit import audit_surface
from tetgenPLCInterface import execute
from reviewSU2ExplicitQuality import quality

def write(p,r):p.write_text(json.dumps(r,indent=2)+'\n')

def main():
    start=time.monotonic();source=ROOT/'outputs/v143_validation/actual';folder=ROOT/'outputs/v144_validation/actual'
    small=ROOT/'docs/validation/v144/fixture/report.json';assert json.loads(small.read_text())['passed']
    folder.mkdir(parents=True,exist_ok=False);fanPath=ROOT/'outputs/v133_validation/fan_workbench/mesh_report.json'
    paths=[Path(__file__).resolve(),small,fanPath,ROOT/'docs/validation/v144/plan.md',source/'hybrid_volume.mat',source/'mesh_process.json',source/'candidate_report.json']
    paths += [ROOT/'scripts'/n for n in ['recenterBaseOuterColumn.py','buildSU2CreaseFan.py','auditSU2HybridTopology.py',
        'meshSU2HybridVolume.py','su2PyramidGeometry.py','runCGALSurfaceAudit.py','cgalSurfaceIntersections.py',
        'tetgenPLCInterface.py','tetgenFixedInterface.py','reviewSU2ExplicitQuality.py']]
    locked={str(p):sha(p) for p in paths};code=1;error=None;plcStarted=False
    try:
        receipt=json.loads((source/'mesh_process.json').read_text());assert receipt['sourceStable'] and receipt['exitCode']==1
        assert sha(source/'hybrid_volume.mat')==receipt['outputSHA256']['hybrid_volume.mat']
        old=loadmat(str(source/'hybrid_volume.mat'));P=old['vertices'];V=old['sourceVertices'];nv=len(V)
        F=old['sourceFaces'].astype(np.int64)-1;W=old['wedges'].astype(np.int64)-1;far=old['farfield'].astype(np.int64)-1
        pids=old['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1;good=pids>=0
        inverse=np.full(int(max(pids))+1,-1,np.int64);inverse[pids[good]]=np.flatnonzero(good)
        p=V[F];av=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);base=av[:,0]/np.linalg.norm(av,axis=1)<-.99
        E,rim=interface_edges(F,base);inside=np.setdiff1d(np.unique(F[base]),rim);assert len(inside)==1 and inside[0]==9876 and sum(base)==46
        col=int(inside[0]);column=inverse[np.arange(33)*nv+col];top=inverse[F+32*nv];edges=inverse[E+32*nv]
        assert min(column)>=0 and top.min()>=0 and edges.min()>=0
        fan=json.loads(fanPath.read_text());fanOld=np.unique(np.r_[np.array(fan['sideNodeIDs']).ravel(),np.array(fan['baseNodeIDs']).ravel()])
        fanIDs=inverse[fanOld];assert min(fanIDs)>=0 and not np.any(np.isin(column,fanIDs))
        h=old['normalHeightsByColumnM'][:,col];Q,offset,detail=recenter(P,column,edges,h,19)
        changed=np.flatnonzero(np.any(Q!=P,axis=1));assert np.array_equal(changed,np.sort(column[20:])) and len(changed)==13
        assert detail['maximumDisplacementM']<=.09
        near=inverse[np.arange(20*nv)];assert min(near)>=0 and np.array_equal(Q[near],P[near])
        assert np.array_equal(Q[fanIDs],P[fanIDs]) and np.array_equal(Q[W[-2852:]],P[W[-2852:]]) and np.array_equal(Q[far],P[far])
        assert np.array_equal(Q[:nv],V)
        shell=dict(vertices=Q,wedges=W+1,tetrahedra=np.empty((0,4),np.int64),farfield=top+1,exteriorFarfield=far+1,
            sourceVertices=V,sourceFaces=F+1,normalHeightsByColumnM=old['normalHeightsByColumnM'],
            baseColumnNodeIDsOneBased=column+1,baseColumnTangentialOffsetsM=offset,baseColumnSourceIDOneBased=col+1,
            columnCoordinateMode=detail['columnCoordinateMode'],sourcePreparedNodeIDsOneBased=pids+1)
        savemat(str(folder/'prepared_cap.mat'),shell,do_compression=True)
        ga=audit(shell);wj=float(min(wedge_jacobian_min(Q[W])));assert wj>0
        write(folder/'shell_geometry_audit.json',dict(ga,minimumPrismJacobian=wj))
        affected=np.any(np.isin(W,changed),axis=1);assert sum(affected)==13*46
        detail.update(baseFaces=46,rimPoints=46,sourceInteriorColumnOneBased=col+1,changedNodeIDsOneBased=(changed+1).tolist(),
            affectedPrisms=int(sum(affected)),wholeTailFanExactlyPreserved=True,protectedNearWallExactlyPreserved=True,
            allConnectivityUnchanged=True,minimumPrismJacobian=wj,preparationSeconds=time.monotonic()-start)
        write(folder/'preparation_report.json',detail)
        boundaryAudit=audit_surface(Q,top,folder/'boundary_audit',dict(source='prepared_cap.mat',sha256=sha(folder/'prepared_cap.mat')))
        assert boundaryAudit['selfIntersectionFree'] and time.monotonic()-start<=120
        ids=np.unique(np.r_[top.ravel(),far.ravel()]);mapping=np.full(len(Q),-1,np.int64);mapping[ids]=np.arange(len(ids))
        boundary=np.r_[top[:,[0,2,1]],far];baseline=json.loads((source/'candidate_report.json').read_text())
        plcStarted=True;R,U,inter=execute(folder/'native',Q[ids],mapping[boundary],np.array(baseline['holeSeed']),np.array(baseline['fluidWitness']),180,paths+[folder/'prepared_cap.mat',folder/'preparation_report.json'])
        free=len(R)-len(ids);localToAll=np.r_[ids,np.arange(len(Q),len(Q)+free,dtype=np.int64)]
        allP=np.r_[Q,R[len(ids):]];allT=localToAll[U];used=np.unique(np.r_[W.ravel(),allT.ravel(),far.ravel()])
        remap=np.full(len(allP),-1,np.int64);remap[used]=np.arange(len(used));wasOld=used<len(Q)
        sourceIDs=np.where(wasOld,used+1,0);preparedIDs=np.zeros(len(used),np.int64);preparedIDs[wasOld]=pids[used[wasOld]]+1
        data=dict(vertices=allP[used],wedges=remap[W]+1,tetrahedra=remap[allT]+1,farfield=remap[far]+1,sourceVertices=V,sourceFaces=F+1,
            normalHeightsByColumnM=old['normalHeightsByColumnM'],baseColumnTangentialOffsetsM=offset,
            baseColumnNodeIDsOneBased=remap[column]+1,baseColumnSourceIDOneBased=col+1,columnCoordinateMode=detail['columnCoordinateMode'],
            sourceV143NodeIDsOneBased=sourceIDs,sourcePreparedNodeIDsOneBased=preparedIDs,
            sourcePLCBoundaryNodeIDsOneBased=ids+1,plcToHybridNodeIDsOneBased=remap[localToAll]+1)
        savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True)
        geom=audit(data);write(folder/'geometry_audit.json',dict(geom,minimumPrismJacobian=wj))
        assert np.array_equal(data['vertices'][data['wedges']-1],Q[W]) and np.array_equal(data['vertices'][data['farfield']-1],P[far])
        values=quality(data['vertices'][data['tetrahedra']-1]);bad=np.flatnonzero(values<.1)
        report=dict(inter,nodes=len(used),prisms=len(W),tetrahedra=len(allT),totalElements=len(W)+len(allT),
            minimumTetSICN=float(min(values)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=values[bad].tolist(),
            badTetCentersM=data['vertices'][data['tetrahedra'][bad]-1].mean(axis=1).tolist(),priorBadTets=5,
            originalWallAreaM2=baseline['originalAreaM2'],originalWallAndNearWallAndTailFanPreserved=True,
            normalHeightMode=detail['columnCoordinateMode'],geometryPassed=True,qualityGatePassed=bool(min(values)>=.1),
            sizeGatePassed=bool(len(used)<=400000 and len(W)+len(allT)<=900000),exteriorAlgorithm='TetGen1.5.1 fixed PLC with tangential base interior outer-column redistribution',
            physicsQualified=False,newQualifiedArea=0)
        report['candidateAccepted']=report['qualityGatePassed'] and report['sizeGatePassed'];write(folder/'candidate_report.json',report)
        assert report['candidateAccepted'],'full exterior quality or size gate failed';code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        r=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,error=error,plcStarted=plcStarted,sourceSHA256=locked,
            sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},
            physicsQualified=False,newQualifiedArea=0);r['sourceStable']=all(r['sourceStableAfter'].values())
        write(folder/'mesh_process.json',r);assert r['sourceStable'];print(json.dumps({k:v for k,v in r.items() if k not in ['sourceSHA256','sourceStableAfter','outputSHA256']},indent=2))
    raise SystemExit(code)

if __name__=='__main__':main()
