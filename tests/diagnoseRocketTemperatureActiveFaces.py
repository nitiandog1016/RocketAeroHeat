"""审查实际Q28全部内部采样温域余量，定位下一步可能受限的区域，不更新流场。"""
import json,os
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,sha
OUT=ROOT/'outputs/v86_validation'


def main():
    state=ROOT/'outputs/v83_validation/krylov28/wall_velocity_krylov_evidence.mat'
    geometry=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat';before={str(p.relative_to(ROOT)):sha(p) for p in (state,geometry)}
    s=loadmat(state,simplify_cells=True,variable_names=['result','gas','boundary']);r=s['result'];gas=s['gas'];Q=r['Q']
    d=loadmat(geometry,simplify_cells=True,variable_names=['domain'])['domain'];a=r['info']['exactEulerFluxAudit'];f=d['openFaces']
    ids=np.asarray(a['faceID'],int)-1;T=np.asarray(a['samplePrimitive'])[:,4];upper=gas['temperatureRange'][1]-T;lower=T-gas['temperatureRange'][0]
    assert np.all(upper>=0)&np.all(lower>=0);nfaces=len(d['inputGeometry']['F']);cellT=r['info']['primitive'][:,4]
    selected=np.unique(np.r_[np.argsort(upper)[:16],19146-1 if 19146-1 in ids else [],19207-1 if 19207-1 in ids else []])
    # 诊断记录的19146/19207是全局面号，不把它们误当排序后的数组位置。
    slots={int(face):i for i,face in enumerate(ids)}
    selected=np.unique(np.r_[np.argsort(upper)[:16],[slots[k-1] for k in (19146,19207)]])
    rows=[]
    for slot in selected:
        face=int(ids[slot]);left=int(f['cellMinus'][face])-1;right=int(f['cellPlus'][face])-1
        rows.append(dict(faceID=face+1,sampleTemperatureK=float(T[slot]),upperMarginK=float(upper[slot]),
                         centroidM=f['centroid'][face].tolist(),ownerCV=left+1,neighbourCV=right+1,
                         ownerLayer=left//nfaces+1,neighbourLayer=right//nfaces+1,
                         ownerSourceFacet=left%nfaces+1,neighbourSourceFacet=right%nfaces+1,
                         ownerCellTemperatureK=float(cellT[left]),neighbourCellTemperatureK=float(cellT[right])))
    out=dict(status='observed-all-actual-Q28-interior-temperature-margins-from-accepted-NS-audit',
             sourceStateSHA256=before[str(state.relative_to(ROOT))],interiorFaceCount=len(ids),
             interfaceTemperatureRangeK=[float(T.min()),float(T.max())],
             marginBandCounts={str(v):int(np.sum(upper<=v)) for v in (.1,1.,5.,10.,50.)},
             hottestAndTrackedFaces=sorted(rows,key=lambda x:x['upperMarginK']),
             actualCellTemperatureRangeK=[float(cellT.min()),float(cellT.max())],
             prospectiveNewtonDirectionOrGlobalStepRadiusClaimed=False,
             newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,physicalSolution=False,
             inputSHA256=before,scriptSHA256=sha(__file__),processID=os.getpid())
    assert all(sha(ROOT/path)==h for path,h in before.items())
    target=OUT/'current_temperature_margins.json';assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','hottestAndTrackedFaces')}))


if __name__=='__main__':main()
