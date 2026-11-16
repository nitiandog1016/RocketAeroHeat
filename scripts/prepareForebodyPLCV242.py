"""模块职责：定位13个坏固定面的母来源，形成保留原壁的前体单壳PLC采用依据。"""
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from phaseEulerRegionalVolume import read_tetra_mesh,topology
from reviewSU2ExplicitQuality import bound
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v242_validation';OLD=ROOT/'outputs/v240_validation'


def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False)
    x,t,m=read_tetra_mesh(OLD/'run01/mesh.su2');topo=topology(x,t,m)
    cavity=np.load(str(ROOT/'outputs/v241_validation/buffer01/cavity.npz'));provenance=np.load(str(OLD/'region01/provenance.npz'))
    assert np.array_equal(x,provenance['xyz']) and np.array_equal(t,provenance['cells'])
    fixed=cavity['fixedFaces'];upper=bound(x[fixed])[0];bad=np.flatnonzero(upper<.1);assert len(bad)==13
    lookup={tuple(sorted(f)):i for i,f in enumerate(topo['faces'])};rows=[]
    for j in bad:
        face=fixed[j];idx=lookup[tuple(sorted(face))];owners=[int(topo['owner'][idx])]
        if topo['other'][idx]>=0:owners.append(int(topo['other'][idx]))
        mothers=provenance['motherCell'][owners];partial=np.isin(mothers,provenance['partialMotherCells'])
        rows.append(dict(faceNodes=face.tolist(),singleFaceMaxSICN=float(upper[j]),owners=owners,motherCells=mothers.tolist(),
          motherPartial=partial.tolist(),generatedNodeCount=int(sum(provenance['originalNode'][face]<0)),
          boundaryMarker=next((tag for tag,v in topo['markers'].items() if idx in v),'artificial-cavity-interface')))
    wall=m['wall'];wn=np.unique(wall);xx=x[wn];mapid=np.full(len(x),-1,np.int64);mapid[wn]=np.arange(len(wn));w=mapid[wall]
    edges=np.concatenate([w[:,[0,1]],w[:,[1,2]],w[:,[2,0]]]);keys,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    assert np.all(counts<=2);ringedges=keys[counts==1];rn=np.unique(ringedges)
    assert len(rn)==46 and np.all(np.unique(ringedges,return_counts=True)[1]==2)
    adjacency={int(i):[] for i in rn}
    for a,b in ringedges:adjacency[int(a)].append(int(b));adjacency[int(b)].append(int(a))
    ring=[int(min(rn))];prev=-1
    while True:
        nex=next(v for v in sorted(adjacency[ring[-1]]) if v!=prev)
        if nex==ring[0]:break
        assert nex not in ring;prev=ring[-1];ring.append(nex)
    assert len(ring)==46;ring=np.array(ring,dtype=np.int64);cut=float(xx[ring[0],0]);assert np.all(xx[ring,0]==cut)
    # 环按yz正向排序仅改变外表面构造遍历，不改变原壁节点/三角面。
    area2=np.sum(xx[ring,1]*np.roll(xx[ring,2],-1)-xx[ring,2]*np.roll(xx[ring,1],-1))
    if area2<0:ring=ring[::-1]
    minimum=x.min(axis=0);maximum=x.max(axis=0);planes=[];covered=np.zeros(len(m['farfield']),bool)
    for axis in range(3):
        for side,value in [('low',minimum[axis]),('high',maximum[axis])]:
            mask=np.all(abs(x[m['farfield'],axis]-value)<1e-12,axis=1)
            planes.append(dict(axis=axis,side=side,coordinate=float(value),faces=int(sum(mask))))
            covered|=mask
    wallbound=bound(xx[w])[0];geom=read(OLD/'region01/report.json')
    expected=float(topo['volume'].sum());sourcearea=float(np.sum(np.linalg.norm(np.cross(xx[w[:,1]]-xx[w[:,0]],xx[w[:,2]]-xx[w[:,0]]),axis=1)/2))
    selected=bool(np.all(covered) and np.all(wallbound>=.1) and all(all(r['motherPartial']) for r in rows))
    mapping=np.loadtxt(str(OLD/'run01/original_faces.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(mapping[:,1:4].astype(int),wall)
    np.savez_compressed(str(out/'wall.npz'),vertices=xx,wall=w,ring=ring,originalVolumeNodeIDs=wn,
      originalFaceIDs=mapping[:,0].astype(np.int64),originalArea=mapping[:,4],bounds=np.array([minimum,maximum]),expectedVolume=expected)
    report=dict(passed=True,adoptSingleShell=selected,badFixedFaces=rows,allBadOwnersFromClippedMothers=all(all(r['motherPartial']) for r in rows),
      allFarfieldFacesOnBoxPlanes=bool(covered.all()),farfieldPlanes=planes,boundsM=[minimum.tolist(),maximum.tolist()],
      wallFaces=len(w),wallNodes=len(wn),wallAreaM2=sourcearea,originalWallFaceIDsExact=True,wallRingNodes=len(ring),cutXM=cut,
      wallMinSingleFaceMaxSICN=float(wallbound.min()),wallBoundsBelow01=int(sum(wallbound<.1)),
      expectedFluidVolumeM3=expected,originalRegionVolumeM3=geom.get('volumeM3'),
      physicalQualification=False,newQualifiedAreaM2=0,sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [
      OLD/'run01/mesh.su2',OLD/'run01/original_faces.csv',OLD/'region01/provenance.npz',ROOT/'outputs/v241_validation/buffer01/cavity.npz']})
    save(out/'report.json',report);print(report)


if __name__=='__main__':raise SystemExit('v242 fixed task closed; use NEXT_DECISION')
