"""模块职责：全原面目标和凸锥整面角域绑定，复用供体，不重算旧波系或剖面。"""
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v255_validation';CLOSED=True


def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False)
    s=read(ROOT/'outputs/v254_validation/prepare01/source.json');old=read(ROOT/'outputs/v252_validation/matlab01/surface_base.json')
    meta=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json');a=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'))
    ids=np.flatnonzero(old['validFace']);axis=-np.array(meta['condition']['flowDirection']);axis/=np.linalg.norm(axis)
    tri=a['originalVertices'][a['originalFaces'][ids]]-np.array(meta['sphereCenter']);center=tri.mean(axis=1)
    angle=np.arctan2(np.linalg.norm(np.cross(center,axis),axis=1),center@axis)
    vertex=np.arctan2(np.linalg.norm(np.cross(tri,axis),axis=2),tri@axis)
    coneOK=np.max(vertex,axis=1)<=np.pi/3
    # 0—60°为凸圆锥；所有顶点在内即可保证每个重心组合在内。无内侧截断。
    probe=np.unique(np.r_[0,np.geomspace(1e-4,.05,65),np.deg2rad(np.arange(.25,60.01,.25))])
    s['targetAnglesRad']=np.r_[s['targetAnglesRad'][:46],angle,probe,np.pi/3].tolist()
    s['parentV254SourceSHA256']=sha(ROOT/'outputs/v254_validation/prepare01/source.json')
    s['surfaceTargetIndicesZeroBased']=list(range(46,46+len(ids)));s['surfaceFaceIDs']=(ids+1).tolist()
    save(out/'source.json',s)
    save(out/'geometry.json',dict(originalFaces=19752,originalAreaM2=float(a['originalFaceArea'].sum()),candidateFaceIDs=(ids+1).tolist(),
        candidateAnglesRad=angle.tolist(),vertexAnglesRad=vertex.tolist(),coneAdmittedFaceIDs=(ids[coneOK]+1).tolist(),
        coneRejectedFaceIDs=(ids[~coneOK]+1).tolist(),coneReason='convex circular cone theta<=60deg; nonnegative barycentric sums preserve membership; full origin-including interval',
        sourceMeshSHA256=sha(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz')))
    save(out/'report.json',dict(passed=True,targetSamples=len(s['targetAnglesRad']),candidateFaces=len(ids),coneAdmittedFaces=int(coneOK.sum()),
        coneRejectedFaces=int((~coneOK).sum()),sourceSHA256=sha(out/'source.json'),oldShockRecomputed=False,oldProfilesRecomputed=False,newQualifiedAreaM2=0))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v255 fixed preparation closed; follow NEXT_DECISION')
    prepare()
