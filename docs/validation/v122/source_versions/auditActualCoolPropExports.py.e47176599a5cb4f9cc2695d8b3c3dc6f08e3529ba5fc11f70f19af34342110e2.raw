"""独立按原VTK三角面核对组件CSV，未知壁温仍未知，方法与转换身份一致。"""
from pathlib import Path
import csv,json,sys,hashlib
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'outputs/v122_validation/review02';folder=ROOT/'outputs/v122_validation/actual04';target=out/'independent_exports.json';assert not target.exists()
    r=loadmat(str(out/'result.mat'),squeeze_me=True,struct_as_record=False)['r'];V=r.V;F=r.F.astype(int)-1
    xyz,f=vtk_tokens(folder/'surface_flow.vtk');distance,ids=cKDTree(V).query(xyz);assert max(distance)<1e-6 and len(set(ids))==len(V)
    node=np.zeros((len(V),3));node[ids]=np.column_stack([f[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']]);face=node[F].mean(axis=1)
    table=np.genfromtxt(str(out/'fields.csv'),names=True,delimiter=',');assert np.array_equal(table['originalFaceID'],np.arange(1,19753))
    for j,key in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):assert np.allclose(table[key],face[:,j],rtol=1e-12,atol=1e-8)
    assert (table['givenWallTemperature_K']==300).all() and np.isnan(table['solvedWallTemperature_K']).all() and not table['productionQualified'].any()
    with (out/'section.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.reader(stream))[1:]
    for row in rows:
        assert abs(float(row[0])+1.5)<1e-14 and abs(float(row[10])-face[int(row[2])-1,1])<1e-8*max(abs(float(row[10])),1)
    for name in ['fields.csv.json','section.csv.json']:
        m=json.loads((out/name).read_text(encoding='utf-8-sig'));a=m['convergenceAudit']
        assert a['checkpointIteration']==19 and a['sourceCumulativeIteration']==799 and 'cumulativeCheckpointIteration' not in a
        assert m['algorithmRuntime']['flavor']=='rah-coolprop-recon-v3' and 'HEOS Air' in m['algorithmRuntime']['physicsCaption']
    report=dict(passed=True,originalFaces=19752,sectionRows=len(rows),newPhysicsLocalIteration=19,sourceCumulativeIteration=799,
        unknownWallTemperaturePreserved=True,physicsQualified=False,newQualifiedArea=0,sources={str(p):sha(p) for p in [Path(__file__).resolve(),out/'result.mat',out/'fields.csv',out/'section.csv',folder/'surface_flow.vtk']})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
