"""独立核对真实原面/突起保留、SU2壁面输出、节点负荷与UI面场/截面导出。"""
from pathlib import Path
import hashlib,json,csv,re,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(p):return loadmat(str(p),squeeze_me=True,struct_as_record=False)

def main():
    out=ROOT/'outputs/v113_validation';assert not (out/'independent_audit.json').exists()
    ui=json.loads((out/'ui01/report.json').read_text(encoding='utf-8'));feature=Path(ui['engineFolder'])
    actual=out/'actual_review02/engine';records={};paths=[Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py']
    mapped=[mat(out/'actual_review02/mapped_result.mat')['r'],mat(feature/'mapped_research.mat')['result']]
    for (name,folder),r in zip([('actual19752',actual),('fairing3432',feature)],mapped):
        V=r.V;F=r.F.astype(int)-1;x,field=vtk_tokens(folder/'surface_flow.vtk')
        distance,ids=cKDTree(V).query(x);assert max(distance)<1e-6 and len(np.unique(ids))==len(V)
        values=np.zeros((len(V),3));values[ids]=np.column_stack([field[k][:,0] for k in ('Pressure','Heat_Flux','Temperature')])
        average=values[F].mean(axis=1)
        for i,key in enumerate(['pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw']):
            assert np.allclose(average[:,i],getattr(r,key),rtol=3e-12,atol=1e-8)
        av=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;area=np.linalg.norm(av,axis=1)
        nav=np.zeros_like(V)
        for j in range(3):np.add.at(nav,F[:,j],av/3)
        native=float(np.linalg.norm(nav,axis=1)@values[:,1]);physical=float(area@average[:,1])
        history=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),delimiter=',',names=True))
        assert np.array_equal(history['Inner_Iter'],[0,1]);assert abs(native/history['HF'][-1]-1)<2e-6
        assert abs(physical/r.targetHeatPowerRaw-1)<1e-12
        rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1);rho=rest[:,4]
        e=rest[:,8]/rho-.5*np.sum((rest[:,5:8]/rho[:,None])**2,axis=1)
        assert np.isfinite(rest).all() and min(rho)>0 and min(e)>0
        receipt=json.loads((folder/'engine_process.json').read_text(encoding='utf-8'))
        assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut']
        for file,h in receipt['outputSHA256'].items():assert sha(folder/file)==h,file
        text=(folder/'engine_stdout.log').read_text();assert 'All volume elements are correctly oriented.' in text and 'All surface elements are correctly oriented.' in text
        orth=re.search(r'Orthogonality Angle \(deg.\)\|\s*([0-9.]+)\|',text)
        records[name]=dict(originalFaces=len(F),bodyArea=float(area.sum()),engineSeconds=receipt['elapsedSeconds'],
                          nativePower=native,physicalAreaPower=physical,minimumDensity=float(min(rho)),minimumInternalEnergy=float(min(e)),
                          minimumNativeOrthogonalityDeg=float(orth.group(1)),rhoLog=float(history['rmsRho'][-1]),
                          energyLog=float(history['rmsRhoE'][-1]),checkpointIteration=int(r.convergenceAudit.checkpointIteration),
                          originalOrientationAccepted=True,physicalSolutionQualified=False)
        paths += [folder/f for f in ('engine_process.json','case.cfg','mesh.su2','history.csv','surface_flow.vtk','restart_flow.csv')]
    original=mat(out/'actual_input.mat');assert np.array_equal(mapped[0].V,original['V']) and np.array_equal(mapped[0].F,original['F'])
    before=mat(out/'ui01/feature_geometry.mat')['before'];assert np.array_equal(mapped[1].V,before.orientedMesh.V)
    peak=float(np.max(before.modelSource.truth.nodeRadialOffset));assert abs(peak-.0192)<1e-15
    table=np.genfromtxt(str(out/'ui01/fields.csv'),delimiter=',',names=True);r=mapped[1]
    assert np.array_equal(table['originalFaceID'],np.arange(1,len(r.F)+1)) and not table['productionQualified'].any()
    assert np.isnan(table['solvedWallTemperature_K']).all() and np.all(table['givenWallTemperature_K']==300)
    for c,k in [('pressureRaw_Pa','pressureRaw'),('heatFluxRaw_W_m2','heatFluxRaw'),('fluidNodeTemperatureMean_K','fluidNodeTemperatureMeanRaw')]:
        assert np.allclose(table[c],getattr(r,k),rtol=1e-12,atol=1e-8)
    with (out/'ui01/section.csv').open(encoding='utf-8-sig',newline='') as stream:sections=list(csv.reader(stream))[1:]
    featureX=float(before.modelSource.truth.features.centerX)
    for row in sections:
        assert abs(float(row[0])-featureX)<1e-12 and abs(float(row[10])-r.heatFluxRaw[int(row[2])-1])<1e-8
    # 开口/倒向原表面应在实际网格器前拒绝；不再运行成功网格任务。
    sys.path.insert(0,str(ROOT/'scripts'));from meshSU2HybridVolume import prepare_layers
    options=json.loads(str(original['optionsJSON']));V=original['V'];F=original['F'].astype(int)-1
    rejects=[]
    for label,bad in [('open',F[:-1]),('inward',F[:,[0,2,1]])]:
        try:prepare_layers(V,bad,options)
        except AssertionError:rejects.append(label)
    assert rejects==['open','inward']
    paths += [out/'actual_input.mat',out/'ui01/feature_geometry.mat',out/'ui01/fields.csv',out/'ui01/section.csv',ROOT/'scripts/meshSU2HybridVolume.py']
    report=dict(cases=records,featurePeakRetainedM=peak,sectionAtFeatureCenterX=featureX,sectionRows=len(sections),
                fullCSVRows=len(table),geometryRejections=rejects,pressureHeatTemperatureIndependentlyMatched=True,
                newQualifiedArea=0,physicsQualified=False,inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    (out/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:report[k] for k in ['cases','featurePeakRetainedM','sectionRows','fullCSVRows']}))

if __name__=='__main__':main()
