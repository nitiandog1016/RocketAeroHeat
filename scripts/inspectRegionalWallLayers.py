"""只读检查原32层中的目标区域；不构造新网格，不改原面或高度。"""
import json,time,shutil
import numpy as np
from scipy.io import loadmat
from prepareOpenFoamControl import ROOT,save,sha

def main():
    start=time.monotonic();out=ROOT/'outputs/v208_validation/layers02';out.mkdir(exist_ok=False)
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';region=ROOT/'outputs/v207_validation/region01/contract.json'
    assert sha(source)=='22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c'
    d=loadmat(str(source),variable_names=['vertices','wedges','sourceVertices','sourceFaces','normalHeightsByColumnM'])
    P=d['vertices'];V=d['sourceVertices'];F=d['sourceFaces'].astype(np.int64)-1;W=d['wedges'].astype(np.int64)-1;n=len(V)
    c=json.loads(region.read_text());ids=np.array(c['wallOriginalFaceIDs'])-1;rim=np.array(c['cutRingOriginalVertexIDs'])-1
    key=lambda f:np.sort(f,axis=1)@np.array([1,n,n*n],np.int64)
    keys=key(F[ids]);regular=(W<33*n).all(axis=1)&np.all(np.sort(W[:,:3]%n,axis=1)==np.sort(W[:,3:]%n,axis=1),axis=1)
    picked=np.flatnonzero(regular & np.isin(key(W[:,:3]%n),keys));local=W[picked];columns=np.unique(F[ids]);heights=d['normalHeightsByColumnM'][:,columns]
    counts=np.unique(key(local[:,:3]%n),return_counts=True)[1]
    assert len(picked)==7838*32 and np.all(counts==32)
    assert np.array_equal(P[:n],V) and np.all(np.diff(heights,axis=0)>0)
    rim_nodes=np.arange(33)[:,None]*n+rim[None,:];cut_error=float(np.max(abs(P[rim_nodes,0]-c['cutXM'])))
    edge=np.linalg.norm(V[F[ids]]-np.roll(V[F[ids]],1,axis=1),axis=2)
    # 用最小实际单元边长比记录各向异性；它不是y+或完整质量资格。
    q=P[local];normal=np.linalg.norm(q[:,3:]-q[:,:3],axis=2)
    tangential=np.linalg.norm(q[:,:3]-np.roll(q[:,:3],1,axis=1),axis=2)
    aspect=tangential.max(axis=1)/normal.min(axis=1)
    np.savez_compressed(out/'selection.npz',originalWedgeRows=picked,originalWallFaceIDs=ids+1,originalColumnIDs=columns+1)
    report=dict(passed=True,sourceSHA256=sha(source),regionSHA256=sha(region),
        wallFaces=len(ids),wallColumns=len(columns),fullLayers=32,protectedLayers=19,selectedWedges=len(local),
        selectedPrismPoints=len(np.unique(local)),firstHeightRangeM=[float(heights[1].min()),float(heights[1].max())],
        outerHeightRangeM=[float(heights[-1].min()),float(heights[-1].max())],
        tangentEdgeRangeM=[float(edge.min()),float(edge.max())],
        maxTangentialToNormalEdgeRatio=float(aspect.max()),medianTangentialToNormalEdgeRatio=float(np.median(aspect)),
        cutRingPlanarityMaxErrorM=cut_error,cutRingNodes=len(rim),allSourceFacesHave32Layers=True,
        normalLayerOutletPlanar=bool(cut_error<1e-12),
        needsLayerFollowingOutlet=bool(cut_error>=1e-12),
        volumeMeshGenerated=False,CFDRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    assert report['elapsedSeconds']<60
    shutil.copyfile(__file__,out/'executed_entry.py');save(out/'report.json',report);print(json.dumps(report))

if __name__=='__main__':main()
