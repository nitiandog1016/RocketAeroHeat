"""独立核对两份真实失败候选、原壁/近壁和MATLAB显示的每一个质量数值/编号。"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from reviewSU2ExplicitQuality import quality
from prepareSU2CheckpointContinuation import sha,read

def main():
    base=ROOT/'outputs/v136_validation';saved=loadmat(str(base/'inspector01/candidate_previews.mat'),squeeze_me=True,struct_as_record=False)
    prepared=loadmat(str(base/'explicit/prepared_transition.mat'));reports={}
    task=read(base/'v136_inspector01_process.json');assert task['exitCode']==0 and all(task['sourceStableAfter'].values())
    for j,name in enumerate(['explicit','hxt']):
        path=base/name/'hybrid_volume.mat';d=loadmat(str(path));P=d['vertices'];T=d['tetrahedra'].astype(int)-1
        W=d['wedges'].astype(int)-1;Q=d['pyramids'].astype(int)-1;offset=len(W)+len(Q);q=quality(P[T]);bad=np.flatnonzero(q<.1)
        info=saved['infos'][j];payload=saved['payloads'][j];ids=info.sourceVolumeNodeIDs.astype(int)-1
        assert np.array_equal(info.sourceTetRows.astype(int)-1,bad) and np.array_equal(info.sourceCellIDs.astype(int)-1,bad+offset)
        error=float(max(abs(q[bad]-info.cellSICN)));assert error<1e-12
        assert np.array_equal(payload.vertices,P[ids]);faces=ids[payload.faces.astype(int)];n=int(info.visibleCellFaces)
        for tri,owner,value in zip(faces[:n],info.displaySourceCellIDs[:n].astype(int),info.displayValues[:n]):
            k=owner-1-offset;assert set(tri)<=set(T[k]) and abs(value-q[k])<1e-12
        context=info.contextOriginalFaceIDs.astype(int)-1
        assert np.array_equal(faces[n:],d['sourceFaces'][context]-1) and np.all(np.isnan(payload.faceField.values[n:]))
        assert np.array_equal(d['sourceVertices'],prepared['sourceVertices']) and np.array_equal(d['sourceFaces'],prepared['sourceFaces'])
        for key in ['wedges','pyramids']:
            assert np.array_equal(P[d[key].astype(int)-1],prepared['vertices'][prepared[key].astype(int)-1])
        assert np.array_equal(P[d['farfield'].astype(int)-1],prepared['vertices'][prepared['exteriorFarfield'].astype(int)-1])
        assert len(np.unique(np.sort(T,axis=1),axis=0))==len(T)
        fixed=np.unique(np.r_[W.ravel(),Q.ravel(),d['farfield'].ravel()-1]);counts=np.isin(T[bad],fixed).sum(axis=1)
        reports[name]=dict(badTets=len(bad),minimumTetSICN=float(min(q)),maximumSVDQualityDifference=error,
            visibleCellFaces=n,contextFaces=len(context),originalWallAndRetainedPrismsAndPyramidsAndFarfieldExact=True,
            exactPreviewNodesFacesAndValues=True,duplicateTetrahedra=0,badTetRowsOneBased=(bad+1).tolist(),
            badTetCentersM=P[T[bad]].mean(axis=1).tolist(),fixedVertexCountHistogram=np.bincount(counts,minlength=5).tolist(),
            noseBadTets=int(sum(P[T[bad]].mean(axis=1)[:,0]>0)),sourceSHA256=sha(path))
    result=dict(candidates=reports,matlabSourceStable=True,physicsQualified=False,newQualifiedArea=0,
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [base/'inspector01/candidate_previews.mat',base/'v136_inspector01_process.json',Path(__file__).resolve(),ROOT/'scripts/reviewSU2ExplicitQuality.py']})
    target=base/'inspector01/independent.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:{n:v for n,v in row.items() if n not in ['badTetRowsOneBased','badTetCentersM','sourceSHA256']} for k,row in reports.items()},indent=2))
if __name__=='__main__':main()
