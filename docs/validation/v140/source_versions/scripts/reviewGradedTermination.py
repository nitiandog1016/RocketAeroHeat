"""只读复核渐变终止的原列身份、独立单元体积/质量与失败层位置。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from prepareCGALRemeshingRuntime import ROOT,sha
from optimizeSU2FanExterior import sicn

def main():
    start=time.monotonic();folder=ROOT/'outputs/v140_validation/actual';source=ROOT/'outputs/v133_validation/fan_workbench'
    paths=[Path(__file__).resolve(),ROOT/'scripts/optimizeSU2FanExterior.py',source/'hybrid_volume.mat',source/'mesh_report.json',
        folder/'prepared_transition.mat',folder/'mesh_process.json',folder/'preparation_report.json']
    locked={str(p):sha(p) for p in paths};process=json.loads((folder/'mesh_process.json').read_text())
    assert process['exitCode']==1 and process['sourceStable'] and not process['gmshStarted']
    for n,h in process['outputSHA256'].items():assert sha(folder/n)==h
    old=loadmat(str(source/'hybrid_volume.mat'));d=loadmat(str(folder/'prepared_transition.mat'))
    P=d['vertices'];V=d['sourceVertices'];F=d['sourceFaces'].astype(np.int64)-1;K=d['columnLayers'].ravel().astype(np.int64)
    W=d['wedges'].astype(np.int64)-1;Q=d['pyramids'].astype(np.int64)-1;T=d['tetrahedra'].astype(np.int64)-1
    assert np.array_equal(P,old['vertices']) and np.array_equal(V,old['sourceVertices']) and np.array_equal(F+1,old['sourceFaces'])
    oldW=old['wedges'].astype(np.int64)-1;assert np.array_equal(W[:19*len(F)],oldW[:19*len(F)]) and np.array_equal(W[-2852:],oldW[-2852:])
    E=np.sort(F[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1);assert max(abs(K[E[:,0]]-K[E[:,1]]))==1
    # 原面逐列取层顶，检查输出包络直接来自原列，而非重新画了另一套物理表面。
    columns=np.r_[oldW[:32*len(F)].reshape(32,len(F),6)[:,:,:3],oldW[:32*len(F)].reshape(32,len(F),6)[-1:,:,3:]]
    expected=columns[K[F],np.arange(len(F))[:,None],np.arange(3)[None,:]]
    assert np.array_equal(d['farfield']-1,expected) and np.array_equal(d['exteriorFarfield'],old['farfield'])
    v=P[T];tv=np.linalg.det(np.stack([v[:,j]-v[:,0] for j in [1,2,3]],axis=2))/6;assert min(tv)>0
    p=P[Q];c=p[:,:4].mean(axis=1);pv=np.zeros(len(Q))
    for j in range(4):pv+=np.linalg.det(np.stack([p[:,j]-c,p[:,(j+1)%4]-c,p[:,4]-c],axis=2))/6
    assert min(pv)>0
    values=sicn(P,T);bad=np.flatnonzero(values<.1);r=json.loads((folder/'preparation_report.json').read_text())
    assert np.array_equal(bad+1,r['badTetRowsOneBased']);difference=float(max(abs(values[bad]-r['badTetSICN'])));assert difference<1e-12
    lineage=d['tetSourceFaceLayer'].astype(np.int64)-1;hist={str(j+1):int(sum(lineage[bad,1]==j)) for j in np.unique(lineage[bad,1])}
    rows=np.arange(len(T));base=d['sourceFaces'][lineage[:,0]]-1
    # 每个四面体底面均是该原面在本层的三个列点；第四点来自唯一继续的原列。
    expectedBase=columns[lineage[:,1],lineage[:,0]];assert np.array_equal(T[:,:3],expectedBase)
    h=np.array(json.loads((source/'mesh_report.json').read_text())['normalHeightsM'])
    lengths=np.linalg.norm(P[T[:,[1,2,0]]]-P[T[:,:3]],axis=2).max(axis=1)
    normals=np.cross(P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    altitude=np.abs(np.einsum('ij,ij->i',P[T[:,3]]-P[T[:,0]],normals));aspect=altitude/lengths
    stable={p:sha(Path(p))==v for p,v in locked.items()};assert all(stable.values())
    report=dict(passed=True,seconds=time.monotonic()-start,sourceSHA256=locked,sourceStableAfter=stable,
        allOriginalPointCoordinatesExactlyEqual=True,originalWallAndFarfieldExactlyEqual=True,protected19LayerPrisms=19*len(F),
        unchangedTailFanPrisms=2852,maximumNeighborLayerDifference=1,outerTrianglesExactlyFromOriginalColumns=True,
        positiveIndependentTetAndPyramidVolumes=True,minimumIndependentTetVolume=float(min(tv)),minimumIndependentPyramidVolume=float(min(pv)),
        independentQualityMaxDifference=difference,minimumTetSICN=float(min(values)),badTets=len(bad),badUpperLayerHistogram=hist,
        badTetAltitudeOverMaximumBaseEdgeRange=[float(min(aspect[bad])),float(max(aspect[bad]))],
        badTetAltitudeRangeM=[float(min(altitude[bad])),float(max(altitude[bad]))],
        badTetMaximumBaseEdgeRangeM=[float(min(lengths[bad])),float(max(lengths[bad]))],
        noExteriorMeshOrFlowStarted=True,scope='local termination quality failure; no inference that any smooth outer envelope must fail',
        physicsQualified=False,newQualifiedArea=0,candidateAccepted=False)
    target=folder/'independent_review.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','sourceStableAfter']},indent=2))

if __name__=='__main__':main()
