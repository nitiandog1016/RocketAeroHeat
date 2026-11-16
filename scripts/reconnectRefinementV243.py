"""模块职责：只重接四个已具名中心八面体对角线，独立核对后生成空间对照网格。"""
from pathlib import Path
import sys
import numpy as np
from scipy.spatial import ConvexHull
from phaseEulerRegionalVolume import topology,extract_region
from tetgenFixedInterface import sorted_cells
from reviewSU2ExplicitQuality import quality
from clipEulerForebody import write_su2
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import runSU2RocketPhaseV222 as guard
BASE=ROOT/'outputs/v243_validation';DOC=ROOT/'docs/validation/v243';RUNTIME=ROOT/'tools/runtime/gmsh/4.15.2'
CLOSED=True


def boundary(t):
    faces=np.concatenate([t[:,v] for v in ([0,1,2],[0,1,3],[0,2,3],[1,2,3])]);keys,count=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True)
    assert max(count)<=2;return keys[count==1]


def fix():
    out=BASE/'reconnected01';out.mkdir(exist_ok=False)
    z=np.load(str(BASE/'refine_audit01/mesh.npz'));review=read(BASE/'diagonal01/report.json');assert review['canMeetOriginalQuality']
    x=z['vertices'];t=z['tetrahedra'].copy();support=z['nodeSourceSupport'];changed=[];proof=[]
    for record in review['badParents']:
        father=record['parent'];rows=np.flatnonzero(z['parentCells']==father)
        middle=rows[np.all(support[t[rows],0]!=support[t[rows],1],axis=1)];assert len(middle)==4
        before=t[middle].copy();best=max(record['variants'],key=lambda a:(a['minimumSICN'],tuple(-i for i in a['diagonalNodes'])))
        after=np.asarray(best['cells'],dtype=np.int64);assert np.array_equal(boundary(before),boundary(after))
        hull=ConvexHull(x[np.unique(before)]).volume
        vol=np.linalg.det(x[after[:,1:]]-x[after[:,0,None]])/6;assert np.all(vol>0) and abs(vol.sum()/hull-1)<1e-12
        t[middle]=after;changed.extend(middle.tolist());proof.append(dict(parent=father,rows=middle.tolist(),diagonal=best['diagonalNodes'],
          unchangedBoundary=True,relativeHullVolumeDifference=float(abs(vol.sum()/hull-1))))
    assert len(changed)==16;unchanged=np.ones(len(t),bool);unchanged[changed]=False;assert np.array_equal(t[unchanged],z['tetrahedra'][unchanged])
    arrays={name:z[name] for name in z.files};arrays.update(tetrahedra=t,SICN=quality(x[t]));np.savez_compressed(str(out/'mesh.npz'),**arrays)
    save(out/'local_proof.json',dict(passed=True,changedCells=changed,allCoordinatesBitEqual=True,otherCellsBitEqual=True,proof=proof))


def audit():
    out=BASE/'reconnected_audit01';out.mkdir(exist_ok=False);z=np.load(str(BASE/'reconnected01/mesh.npz'));x=z['vertices'];t=z['tetrahedra']
    original=np.load(str(BASE/'refine_audit01/mesh.npz'));proof=read(BASE/'reconnected01/local_proof.json')
    assert np.array_equal(x,original['vertices']) and len(proof['changedCells'])==16
    assert all(np.array_equal(z[k],original[k]) for k in ('wall','outlet','farfield','parentCells','parentWallRows','originalFaceIDs','childArea','parentOriginalFaceIDs','parentOriginalArea','nodeSourceSupport'))
    topo=topology(x,t,{name:z[name] for name in ('wall','outlet','farfield')});closed=extract_region(topo,np.ones(len(t),bool))['summary'];assert closed['components']==1
    oldt=original['tetrahedra'];oldvol=np.linalg.det(x[oldt[:,1:]]-x[oldt[:,0,None]])/6
    parent=z['parentCells'];vsum=np.bincount(parent,weights=topo['volume']);oldsum=np.bincount(parent,weights=oldvol);assert np.max(abs(vsum/oldsum-1))<2e-12
    # 独立成熟实现只读质量核对；不生成或优化任何节点/单元。
    sys.path.insert(0,str(RUNTIME));import gmsh
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1)
    try:
        gmsh.model.add('v243-diagonal-audit');gmsh.model.addDiscreteEntity(3,1)
        gmsh.model.mesh.addNodes(3,1,np.arange(1,len(x)+1),x.ravel());tags=np.arange(1,len(t)+1)
        gmsh.model.mesh.addElementsByType(1,4,tags,(t+1).ravel());q=gmsh.model.mesh.getElementQualities(tags,'minSICN')
    finally:gmsh.finalize()
    assert min(q)>=.1 and np.max(abs(q-z['SICN']))<1e-10
    np.savez_compressed(str(out/'independent_quality.npz'),gmshMinSICN=q)
    write_su2(out/'mesh.su2',dict(xyz=x,cells=t,markers={name:z[name] for name in ('wall','outlet','farfield')}))
    save(out/'report.json',dict(passed=True,nodes=len(x),cells=len(t),originalParents=4434,computationalWallFaces=len(z['wall']),
      originalGeometryAndParentsPreserved=True,allCoordinatesBitEqual=True,changedCentralCells=16,unchangedCells=len(t)-16,
      motherParentVolumeMaxRelativeDifference=float(np.max(abs(vsum/oldsum-1))),minimumIndependentSICN=float(min(q)),
      independentQualityMaxDifference=float(np.max(abs(q-z['SICN']))),cellsBelow01=int(sum(q<.1)),geometry=closed,
      outletNodes=len(np.unique(z['outlet'])),meshSHA256=sha(out/'mesh.su2'),physicalQualification=False,newQualifiedAreaM2=0,CFDAdmitted=False))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 reconnect closed; use NEXT_DECISION')
    if sys.argv[1]=='task_fix':fix()
    elif sys.argv[1]=='task_audit':audit()
    else:
        stage=sys.argv[1];assert stage in ('fix','audit');guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
        guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        locks=[Path(__file__).resolve(),DOC/'PLAN.md',DOC/'METHOD_DECISION.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/phaseEulerRegionalVolume.py',
          ROOT/'scripts/reviewSU2ExplicitQuality.py',ROOT/'scripts/tetgenFixedInterface.py',ROOT/'scripts/clipEulerForebody.py',BASE/'diagonal01/report.json',BASE/'refine_audit01/mesh.npz']
        if stage=='audit':locks += [BASE/'reconnected01/mesh.npz',BASE/'reconnected01/local_proof.json']+[RUNTIME/n for n in ('gmsh.py','gmsh-4.15.dll','runtime.json','LICENSE.txt')]
        # 附加方法决定每次按原字节保存，之后允许写入实际结果而不改执行回执。
        snap=BASE/('inputs_reconnect_'+stage);snap.mkdir(exist_ok=False)
        for p in locks:
            if p.suffix in ('.py','.md'):(snap/(sha(p)+p.suffix)).write_bytes(p.read_bytes())
        raise SystemExit(0 if guard.supervise('reconnect_'+stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],20 if stage=='fix' else 30,256,locks,BASE,'offline' if stage=='fix' else 'audit') else 1)
