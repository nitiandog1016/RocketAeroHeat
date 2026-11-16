"""模块职责：一次成熟Gmsh细化及原面父子审核，不改变原用户壁面数据。"""
from pathlib import Path
import sys
import numpy as np
from scipy.spatial import cKDTree
from phaseEulerRegionalVolume import read_tetra_mesh,topology,extract_region
from reviewSU2ExplicitQuality import quality
from clipEulerForebody import write_su2
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v243_validation';DOC=ROOT/'docs/validation/v243';OLD=ROOT/'outputs/v242_validation'
RUNTIME=ROOT/'tools/runtime/gmsh/4.15.2';CLOSED=False


def refine():
    """模块一：唯一refine调用；全部真实结果先落盘，不作优化或质量筛删。"""
    out=BASE/'refine01';out.mkdir(exist_ok=False)
    rt=read(RUNTIME/'runtime.json')
    for row in rt['files']:assert sha(RUNTIME/row['path'])==row['sha256']
    x,t,m=read_tetra_mesh(OLD/'run01/mesh.su2');sys.path.insert(0,str(RUNTIME));import gmsh
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('Mesh.Renumber',0)
    try:
        gmsh.model.add('fixed-physical-parent-forebody')
        tags={'wall':1,'outlet':2,'farfield':3}
        for tag in tags.values():gmsh.model.addDiscreteEntity(2,tag)
        gmsh.model.addDiscreteEntity(3,1,[1,2,3]);gmsh.model.mesh.addNodes(3,1,np.arange(1,len(x)+1),x.ravel())
        first=1
        for name,tag in tags.items():
            gmsh.model.mesh.addElementsByType(tag,2,np.arange(first,first+len(m[name])),(m[name]+1).ravel());first+=len(m[name])
        gmsh.model.mesh.addElementsByType(1,4,np.arange(first,first+len(t)),(t+1).ravel());gmsh.model.mesh.reclassifyNodes()
        gmsh.model.mesh.refine()
        nt,xx,_=gmsh.model.mesh.getNodes();kind,ct,cc=gmsh.model.mesh.getElements(3,1);assert list(kind)==[4]
        arrays=dict(nodeTags=nt,vertices=xx.reshape(-1,3),cellTags=ct[0],tetrahedronNodeTags=cc[0].reshape(-1,4),gmshMinSICN=gmsh.model.mesh.getElementQualities(ct[0],'minSICN'))
        for name,tag in tags.items():
            k,et,con=gmsh.model.mesh.getElements(2,tag);assert list(k)==[2];arrays[name+'NodeTags']=con[0].reshape(-1,3)
        np.savez_compressed(str(out/'native.npz'),**arrays)
        save(out/'report.json',dict(passed=True,nodes=len(nt),cells=len(ct[0]),wallFaces=len(arrays['wallNodeTags']),
          minimumNativeSICN=float(min(arrays['gmshMinSICN'])),nativeCellsBelow01=int(sum(arrays['gmshMinSICN']<.1)),
          refineCalls=1,optimizationCalls=0,nativeSavedBeforeGeometryAudit=True,physicalQualification=False,newQualifiedAreaM2=0))
    finally:gmsh.finalize()


def audit():
    """模块二：逐原边中点和每母面/母体积核对，保存可聚合回原面的真实计算网格。"""
    out=BASE/'refine_audit01';out.mkdir(exist_ok=False);raw=np.load(str(BASE/'refine01/native.npz'))
    oldx,oldt,oldm=read_tetra_mesh(OLD/'run01/mesh.su2');source=np.load(str(OLD/'tetgen01/mesh.npz'))
    tags=raw['nodeTags'];order=np.argsort(tags);tags=tags[order];x=raw['vertices'][order]
    assert len(np.unique(tags))==len(tags)
    def mapped(nodes):
        ids=np.searchsorted(tags,nodes);assert np.array_equal(tags[ids],nodes);return ids.astype(np.int64)
    t=mapped(raw['tetrahedronNodeTags']);m={name:mapped(raw[name+'NodeTags']) for name in oldm}
    original=mapped(np.arange(1,len(oldx)+1));assert np.array_equal(x[original],oldx)
    edges=np.unique(np.sort(np.concatenate([oldt[:,v] for v in ([0,1],[0,2],[0,3],[1,2],[1,3],[2,3])]),axis=1),axis=0)
    expected=np.r_[oldx,oldx[edges].mean(axis=1)];support=np.r_[np.column_stack([np.arange(len(oldx))]*2),edges]
    delta,match=cKDTree(expected).query(x);assert len(x)==len(expected) and len(np.unique(match))==len(expected)
    assert max(delta)<2e-14;supports=support[match]
    def parents(elements,old):
        table={tuple(sorted(row)):j for j,row in enumerate(old)}
        return np.array([table[tuple(sorted(set(supports[row].ravel())))] for row in elements],dtype=np.int64)
    parent_t=parents(t,oldt);assert np.all(np.bincount(parent_t,minlength=len(oldt))==8)
    parent_f={name:parents(f,oldm[name]) for name,f in m.items()}
    assert all(np.all(np.bincount(v,minlength=len(oldm[k]))==4) for k,v in parent_f.items())
    dx=x[t[:,1:]]-x[t[:,0,None]];vol=np.linalg.det(dx)/6;assert np.all(vol>0)
    oldvol=np.linalg.det(oldx[oldt[:,1:]]-oldx[oldt[:,0,None]])/6
    diff=abs(np.bincount(parent_t,weights=vol,minlength=len(oldt))/oldvol-1);assert max(diff)<2e-10
    area_errors={};areas={}
    for name,f in m.items():
        points=x[f];areas[name]=np.linalg.norm(np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0]),axis=1)/2
        old=oldx[oldm[name]];a=np.linalg.norm(np.cross(old[:,1]-old[:,0],old[:,2]-old[:,0]),axis=1)/2
        e=abs(np.bincount(parent_f[name],weights=areas[name],minlength=len(a))/a-1);area_errors[name]=float(max(e));assert max(e)<2e-10
    topo=topology(x,t,m);closed=extract_region(topo,np.ones(len(t),bool))['summary'];assert closed['components']==1
    assert abs(closed['volumeM3']/float(source['expectedVolume'])-1)<2e-13
    q=quality(x[t]);assert max(abs(q-raw['gmshMinSICN']))<1e-10
    sizeok=len(t)<=200000 and len(x)<=50000;qualityok=bool(q.min()>=.1)
    np.savez_compressed(str(out/'mesh.npz'),vertices=x,tetrahedra=t,wall=m['wall'],outlet=m['outlet'],farfield=m['farfield'],
      parentCells=parent_t,parentWallRows=parent_f['wall'],originalFaceIDs=source['originalFaceIDs'][parent_f['wall']],
      childArea=areas['wall'],parentOriginalFaceIDs=source['originalFaceIDs'],parentOriginalArea=source['originalArea'],SICN=q,nodeSourceSupport=supports)
    write_su2(out/'mesh.su2',dict(xyz=x,cells=t,markers=m))
    save(out/'report.json',dict(passed=sizeok and qualityok,geometryPassed=True,qualityPassed=qualityok,sizePassed=sizeok,
      nodes=len(x),cells=len(t),computationalWallFaces=len(m['wall']),originalWallParents=len(oldm['wall']),
      originalWallNodesBitEqual=True,allNewNodesFromOriginalEdgeMidpoints=True,maxMidpointCoordinateDifferenceM=float(max(delta)),
      everyParentTetHasEightChildren=True,everyParentFaceHasFourChildren=True,maxParentVolumeRelativeDifference=float(max(diff)),
      maxParentAreaRelativeDifference=area_errors,minimumIndependentSICN=float(q.min()),cellsBelow01=int(sum(q<.1)),
      outletNodes=len(np.unique(m['outlet'])),geometry=closed,originalRegionAreaM2=float(source['originalArea'].sum()),
      fullOriginalFaces=19752,fullOriginalAreaM2=1.9646292565715553,CFDAdmitted=False,physicalQualification=False,newQualifiedAreaM2=0,
      meshSHA256=sha(out/'mesh.su2')))


def main(stage):
    if CLOSED:raise SystemExit('v243 fixed refinement closed; use NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    inputs=[Path(__file__).resolve(),DOC/'PLAN.md',DOC/'METHOD_DECISION.md',ROOT/'scripts/runSU2RocketPhaseV222.py',
      ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',ROOT/'scripts/clipEulerForebody.py',
      OLD/'run01/mesh.su2',OLD/'tetgen01/mesh.npz',BASE/'inspect02/report.json']+[RUNTIME/n for n in ('runtime.json','gmsh.py','gmsh-4.15.dll','LICENSE.txt')]
    if stage=='audit':inputs += [BASE/'refine01/native.npz',BASE/'refine01/report.json']
    assert stage in ('mesh','audit')
    snap=BASE/('inputs_refine_'+stage);snap.mkdir(exist_ok=False)
    for p in inputs:
        if p.suffix in ('.md','.py'):(snap/(sha(p)+p.suffix)).write_bytes(p.read_bytes())
    return guard.supervise('refine_'+stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],60,256,inputs,BASE,'mesh' if stage=='mesh' else 'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 fixed refinement closed; use NEXT_DECISION')
    if sys.argv[1]=='task_mesh':refine()
    elif sys.argv[1]=='task_audit':audit()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
