"""模块职责：用成熟Gmsh只读质量API独立审核新前体，核对文本、原壁和出口。"""
import sys
import numpy as np
from phaseEulerRegionalVolume import read_tetra_mesh,topology,extract_region
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v242_validation'


def audit():
    out=BASE/'mesh_audit01';out.mkdir(exist_ok=False)
    z=np.load(str(BASE/'tetgen01/mesh.npz'));x,t,m=read_tetra_mesh(BASE/'tetgen01/mesh.su2')
    assert np.array_equal(x,z['vertices']) and np.array_equal(t,z['tetrahedra'])
    assert all(np.array_equal(m[k],z[k]) for k in ('wall','outlet','farfield'))
    oldx,_,oldm=read_tetra_mesh(ROOT/'outputs/v240_validation/run01/mesh.su2')
    oldmap=np.loadtxt(str(ROOT/'outputs/v240_validation/run01/original_faces.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(x[m['wall']],oldx[oldm['wall']])
    assert np.array_equal(z['originalFaceIDs'],oldmap[:,0]) and np.array_equal(z['originalArea'],oldmap[:,4])
    area=np.linalg.norm(np.cross(x[m['wall'][:,1]]-x[m['wall'][:,0]],x[m['wall'][:,2]]-x[m['wall'][:,0]]),axis=1)/2
    assert np.max(abs(area-z['originalArea']))<1e-16
    topo=topology(x,t,m);closed=extract_region(topo,np.ones(len(t),bool))['summary']
    assert closed['components']==1 and abs(closed['volumeM3']/float(z['expectedVolume'])-1)<2e-14
    outlet=np.unique(m['outlet']);cut=x[outlet[0],0];assert np.all(x[outlet,0]==cut)
    normals=topo['normal'][topo['markers']['outlet']]
    assert np.all(normals[:,0]<0) and np.all(normals[:,1:]==0)
    # Gmsh只导入已有节点/四面体；不生成、优化或改变任何一个网格量。
    runtime=ROOT/'tools/runtime/gmsh/4.15.2';sys.path.insert(0,str(runtime));import gmsh
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1)
    try:
        gmsh.model.add('v242-read-only');gmsh.model.addDiscreteEntity(3,1)
        gmsh.model.mesh.addNodes(3,1,np.arange(1,len(x)+1),x.ravel())
        tags=np.arange(1,len(t)+1);gmsh.model.mesh.addElementsByType(1,4,tags,(t+1).ravel())
        independent=gmsh.model.mesh.getElementQualities(tags,'minSICN')
    finally:gmsh.finalize()
    diff=float(np.max(abs(independent-z['SICN'])));assert diff<2e-12 and independent.min()>=.1
    np.savez_compressed(str(out/'quality.npz'),gmshMinSICN=independent)
    save(out/'report.json',dict(passed=True,nodes=len(x),cells=len(t),wallFaces=len(m['wall']),
      wallAreaM2=float(area.sum()),originalWallCoordinatesAndIDsExact=True,textRoundtripExact=True,
      originalAreaMaxDifferenceM2=float(np.max(abs(area-z['originalArea']))),outletNodes=len(outlet),
      outletFaces=len(m['outlet']),cutXM=float(cut),allOutletNormalsNegativeX=True,
      gmshVersion=gmsh.__version__,independentMinSICN=float(independent.min()),
      maxSICNDifference=diff,cellsBelow01=int(sum(independent<.1)),geometry=closed,
      physicalQualification=False,newQualifiedAreaM2=0,meshSHA256=sha(BASE/'tetgen01/mesh.su2')))


if __name__=='__main__':audit()
