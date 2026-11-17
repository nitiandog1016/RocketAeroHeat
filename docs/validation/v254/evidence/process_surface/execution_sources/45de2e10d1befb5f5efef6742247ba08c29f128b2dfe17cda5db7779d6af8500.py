"""模块职责：因科研图发现起步低谷，复用两档已存解筛查整面映射；不重跑数值解。"""
import json
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v254_validation';CLOSED=False


def audit():
    out=BASE/'surface_audit01';out.mkdir(exist_ok=False)
    base=read(BASE/'matlab02/base.json');fine=read(BASE/'matlab02/fine.json');t=np.array(base['angleRad'])
    errors={};bad=np.zeros(len(t),bool)
    for key in ('heatFluxWm2','massThicknessM'):
        y=np.array(base[key]);ref=np.interp(t,fine['angleRad'],fine[key]);e=abs(y/ref-1);bad|=e>.02;errors[key]=e.tolist()
    # 驻点到首超阈值站的插值区也没有完整证据；将其到最后坏站的下一站整段拒绝。
    last=int(np.flatnonzero(bad)[-1]);assert last+1<len(t);start=t[last+1]
    old=read(ROOT/'outputs/v252_validation/matlab01/surface_base.json');a=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'))
    meta=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json');candidate=np.flatnonzero(old['validFace'])
    tri=a['originalVertices'][a['originalFaces'][candidate]]-np.array(meta['sphereCenter']);center=tri.mean(axis=1)
    cn=center/np.linalg.norm(center,axis=1)[:,None];vn=tri/np.linalg.norm(tri,axis=2)[:,:,None]
    th=np.arccos(cn@(-np.array(meta['condition']['flowDirection'])))
    # atan2叉积/点积给稳定夹角；球面凸三角形包含在各顶点定义的中心球帽内。
    angularRadius=np.max(np.arctan2(np.linalg.norm(np.cross(vn,cn[:,None,:]),axis=2),np.einsum('ijk,ik->ij',vn,cn)),axis=1)
    lower=th-angularRadius;upper=th+angularRadius
    supported=(lower>=start)&(upper<=t[-1]);faceErrors={};centerBad=np.zeros(len(candidate),bool)
    for key in errors:
        y=np.interp(th,t,base[key]);ref=np.interp(th,fine['angleRad'],fine[key]);e=abs(y/ref-1);centerBad|=e>.02;faceErrors[key]=e
    assert np.all(~centerBad[supported])
    report=dict(passed=True,purpose='additional surface screening after visible startup defect; prior ring audit unchanged',
        fullCapNumericalGatePassed=False,rawCandidateFaces=len(candidate),centerFailedFaces=(candidate[centerBad]+1).tolist(),
        maximumCommonStationRelativeDifferences={k:max(v) for k,v in errors.items()},
        maximumMappedCenterRelativeDifferences={k:float(v.max()) for k,v in faceErrors.items()},
        withheldStartIntervalDeg=[0,float(np.rad2deg(start))],wholeTriangleBound='theta_center +/- max angular distance to three vertex rays; conservative enclosing spherical cap',
        retainedFaceIDs=(candidate[supported]+1).tolist(),withheldFaceIDs=(candidate[~supported]+1).tolist(),
        retainedFaces=int(supported.sum()),withheldFaces=int((~supported).sum()),
        retainedAreaM2=float(a['originalFaceArea'][candidate[supported]].sum()),withheldAreaM2=float(a['originalFaceArea'][candidate[~supported]].sum()),
        maximumRetainedCenterRelativeDifferences={k:float(v[supported].max()) for k,v in faceErrors.items()},
        originalFaces=19752,originalAreaM2=1.9646292565715553,physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,
        baseSHA256=sha(BASE/'matlab02/base.json'),fineSHA256=sha(BASE/'matlab02/fine.json'),numericalSolutionRecomputed=False,
        sourceDefect='Large first trapezoid from 0.001 rad to first uniform angle; I starts proportional to theta^4; cause consistent with observed refinement, no corrected production run yet')
    save(out/'report.json',report);save(out/'station_errors.json',dict(angleRad=t.tolist(),relativeDifferences=errors))
    print(json.dumps({k:report[k] for k in ('passed','fullCapNumericalGatePassed','centerFailedFaces','withheldFaces','retainedFaces','withheldStartIntervalDeg')}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v254 fixed surface screening closed; follow NEXT_DECISION')
    audit()
