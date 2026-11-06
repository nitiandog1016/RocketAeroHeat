"""模块用途。一次主工作台只读验收及独立导出，按实际终态选取完整场或小例。"""
import csv
import sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from auditSU2SSTState import vtk_tokens
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import run,save


def main():
    base=ROOT/'outputs/v167_validation'
    assert read(base/'review01/report.json')['passedChecks']
    # WebGL并发文件只作本地准确来源快照，不混入公开提交。
    for n in ('web/viewer.js','web/section.js','tests/verifyWebSectionExport.py'):
        p=base/'collaborator_sources'/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((ROOT/n).read_bytes())
    command=[sys.executable,str(ROOT/'scripts/runMatlabResearch.py'),'main01',
        "addpath(pwd); addpath(fullfile(pwd,'tests')); testSU2SSTEnergyWorkbench(fullfile(pwd,'outputs','v167_validation','main01'));",
        '--folder','outputs/v167_validation/matlab']
    for n in ('outputs/v167_validation/real01/engine_process.json','outputs/v167_validation/small01/engine_process.json',
              'outputs/v167_validation/small01/case_manifest.json','outputs/v167_validation/small01/audit002/report.json',
              'outputs/v157_validation/fixture_workbench/mesh_input.mat','outputs/v145_validation/mesh/hybrid_volume.mat'):
        assert (ROOT/n).is_file();command+=['--input',n]
    run(base,'main01_parent',command,300)
    receipt=read(base/'matlab/main01_process.json');assert receipt['exitCode']==0 and all(receipt['sourceStableAfter'].values())
    result=read(base/'main01/report.json')
    if result['fullDomainDisplayed']:
        sys.path.insert(0,str(ROOT/'tests'))
        from verifySU2IdealAirFullDomain import main as verify
        verify(base,'main01',199,200,('OTH','k/ω共同松弛','舍入余量向内收紧'))
        return
    source=ROOT/'outputs/v157_validation/fixture_workbench/mesh_input.mat';mesh=loadmat(str(source),squeeze_me=True)
    V,F=mesh['V'],mesh['F'].astype(int)-1
    xyz,fields=vtk_tokens(base/'small01/surface_flow.vtk');distance,ids=cKDTree(V).query(xyz)
    assert max(distance)<1e-12 and len(set(ids))==len(V)
    with (base/'main01/fields.csv').open(encoding='utf8') as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==len(F)==8;errors={}
    for name,native in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),('fluidNodeTemperatureMean_K','Temperature')]:
        values=np.empty(len(V));values[ids]=fields[native][:,0];expected=values[F].mean(axis=1)
        actual=np.array([float(row[name]) for row in rows]);errors[name]=float(max(abs(actual-expected)/np.maximum(1,abs(expected))))
        assert errors[name]<1e-12
    assert all(float(r['givenWallTemperature_K'])==300 and np.isnan(float(r['solvedWallTemperature_K'])) for r in rows)
    with (base/'main01/section.csv').open(encoding='utf-8-sig') as stream:section=list(csv.reader(stream))
    assert len(section)==5 and read(base/'main01/section.csv.json')['field']=='heatFluxRaw'
    for row in section[1:]:
        face=int(row[2])-1;assert abs(float(row[10])-float(rows[face]['heatFluxRaw_W_m2']))<1e-12*max(1,abs(float(row[10])))
        assert all(w in row[12] for w in ('k/ω共同松弛','SST约束候选冷重启','算子改变','非无中断等价'))
    for n in ('fields.csv.json','section.csv.json'):
        m=read(base/'main01'/n);a=m['convergenceAudit']
        assert m['algorithmRuntime']['flavor']=='rah-sst-energy-v1' and not m['productionQualified']
        assert a['checkpointIteration']==2 and a['completedUpdatesBeforeRun']==4 and a['cumulativeCheckpointIteration']==6
    save(base/'export_check.json',dict(passed=True,originalFaces=8,sectionSegments=4,fieldRelativeErrors=errors,
        methodAndUnknownTwVerified=True,inputSHA256={str(p):sha(p) for p in (source,base/'main01/fields.csv',base/'main01/section.csv')},
        physicsQualified=False,newQualifiedArea=0))
    print(errors)


if __name__=='__main__':main()
