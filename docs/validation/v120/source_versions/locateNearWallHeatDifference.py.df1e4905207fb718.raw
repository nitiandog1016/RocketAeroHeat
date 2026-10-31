"""按原来流几何迎/背风和固定轴向分箱定位未收敛两网格差异，不诊断真实分离。"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v120_validation'
def load(p):return loadmat(str(p),squeeze_me=True,struct_as_record=False)
def main():
    target=OUT/'heat_difference_location.json';assert not target.exists()
    source=ROOT/'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat';old=load(source)['result']
    a=load(OUT/'review01/control_result.mat')['r'];b=load(OUT/'review01/candidate_result.mat')['r']
    assert np.array_equal(a.V,b.V) and np.array_equal(a.V,old.V) and np.array_equal(a.F,old.F) and np.array_equal(a.F,b.F)
    area=a.faceArea;qa=a.heatFluxRaw;qb=b.heatFluxRaw;absolute=area*abs(qb-qa);den=sum(area*abs(qa))
    def record(mask,label):
        return dict(label=label,faces=int(sum(mask)),areaM2=float(sum(area[mask])),
          shareOfGlobalAbsoluteHeatDifference=float(sum(absolute[mask])/sum(absolute)),
          ownAreaHeatRelativeL1=float(sum(absolute[mask])/max(sum(area[mask]*abs(qa[mask])),1e-30)),
          controlPowerW=float(area[mask]@qa[mask]),candidatePowerW=float(area[mask]@qb[mask]))
    incidence=[record(old.aerothermalRegion.flowClass==k,str(k)) for k in np.unique(old.aerothermalRegion.flowClass)]
    x=old.center[:,0];edges=np.linspace(min(old.V[:,0]),max(old.V[:,0]),9);bins=[]
    for k in range(8):
        mask=(x>=edges[k])&((x<edges[k+1]) if k<7 else (x<=edges[k+1]));item=record(mask,str(k))
        item['xRangeM']=edges[k:k+2].tolist();bins.append(item)
    assert sum(q['faces'] for q in bins)==19752
    result=dict(meaning='localization of two unconverged fields only; geometric incidence is not viscous separation or eligibility',
       referenceFrame='X axis; nose +X; original alpha20 freestream incidence',globalHeatRelativeL1=float(sum(absolute)/den),
       geometricIncidence=incidence,eightEqualAxialBins=bins,physicsQualified=False,newQualifiedArea=0,
       sourceSHA256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),source,OUT/'review01/control_result.mat',OUT/'review01/candidate_result.mat']})
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(incidence=incidence,largestAxialBin=max(bins,key=lambda q:q['shareOfGlobalAbsoluteHeatDifference'])),indent=2))
if __name__=='__main__':main()
