"""准备固定源的一次原生更新记录；运行包变化只为观察，数值配置保持。"""
from pathlib import Path
import argparse,copy,json,os
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import prepare as checkpoint,OUTPUT_CONTROL,validate_source_artifacts
from prepareSU2SSTCase import set_value
from su2SSTInterface import SPATIAL_OUTPUT

BASE=ROOT/'outputs/v147_validation'
PLAN=ROOT/'docs/validation/v147/PLAN.md'
FLAVOR='rah-update-trace-v1'
STUDY='v147-native-implicit-update-observation'
VARIANTS=('small_official','small_off','small_on','actual')

def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def source_path(variant):
    assert variant in VARIANTS
    return ROOT/('outputs/v146_validation/control' if variant=='actual' else 'outputs/v145_validation/fixture02/case')

def expected_settings(variant):
    c=config(source_path(variant)/'case.cfg')
    c.update(ITER='1',OUTPUT_WRT_FREQ='1',SCREEN_WRT_FREQ_INNER='1',
        VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL, TIMESTEP)',
        HISTORY_OUTPUT='(ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)')
    return c

def validate_case(folder,data,c):
    """只允许三个事前小例及一个整域诊断，不开放记录包作普通默认求解。"""
    spec=data['sst']['updateTraceStudy'];variant=spec['variant'];assert variant in VARIANTS
    assert spec['studyID']==STUDY and spec['planSHA256']==sha(PLAN)
    assert spec['enabled']==(variant in ('small_on','actual'))
    assert data.get('requiredRuntimeFlavor','official')==('official' if variant=='small_official' else FLAVOR)
    assert c==expected_settings(variant)
    init=data['initialization'];assert Path(init['sourceDirectory']).resolve()==source_path(variant).resolve()
    assert init['sourceLocalIteration']==(99 if variant=='actual' else 1)
    assert init['diagnosticOnly'] and init['diagnosticPurpose']=='spatial-implicit-RHS-and-pseudotime'
    locks=validate_source_artifacts(folder,init)
    if variant=='actual':
        proof=BASE/'small_comparison.json';assert sha(proof)==spec['compatibilitySHA256'] and read(proof)['passed']
        locks[str(proof)]=sha(proof)
    for p in [PLAN,Path(__file__).resolve()]:locks[str(p)]=sha(p)
    return locks

def validate_runtime(folder,data,c,runtime):
    assert runtime['flavor']==FLAVOR and runtime['sourceCommit']=='7048fd35f56d4bf52e586429d6746a36a10255f8'
    assert runtime['allCFDObjectsRebuilt'] and runtime['priorThermalAndFrozenHooksRemoved']
    assert runtime['buildReceipt']['exitCode']==0 and not runtime['buildReceipt']['timedOut']
    proof=read(BASE/'preflight/report.json')
    assert runtime['patchSHA256']==proof['patchSHA256']==sha(ROOT/'third_party/su2_update_trace/update_trace.diff')
    locks=validate_case(folder,data,c)
    locks[str(BASE/'preflight/report.json')]=sha(BASE/'preflight/report.json')
    return locks

def prepare(variant):
    source=source_path(variant);folder=BASE/variant;assert not folder.exists()
    if variant=='actual':assert read(BASE/'small_comparison.json')['passed']
    folder.mkdir(parents=True)
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    for key,value in expected_settings(variant).items():text=set_value(text,key,value)
    os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'))
    contract=copy.deepcopy(read(source/'sst_contract.json'))
    for key in ['numericalStudy','initializationStudy']:contract.pop(key,None)
    contract.update(spatialOutputContract=SPATIAL_OUTPUT,spatialDiagnosticOnly=True)
    manifest=copy.deepcopy(read(source/'case_manifest.json'))
    manifest['flight'].update(iterations=1,outputEvery=1,namedCheckpoints=True)
    checkpoint(folder,source,text,contract,manifest)
    init=read(folder/'sst_initialization.json')
    init.update(diagnosticOnly=True,diagnosticPurpose='spatial-implicit-RHS-and-pseudotime',runtimeObservationOnly=True)
    write(folder/'sst_initialization.json',init)
    spec=dict(studyID=STUDY,variant=variant,enabled=variant in ('small_on','actual'),planSHA256=sha(PLAN))
    if variant=='actual':spec['compatibilitySHA256']=sha(BASE/'small_comparison.json')
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),updateTraceStudy=spec)
    write(folder/'sst_contract.json',contract)
    manifest.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'),
                    requiredRuntimeFlavor='official' if variant=='small_official' else FLAVOR)
    validate_case(folder,manifest,config(folder/'case.cfg'))
    write(folder/'case_manifest.json',manifest)
    write(folder/'engine_request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        threads=2,maximumSeconds=180 if variant=='actual' else 60,dryRun=False,
        sourceFiles=['scripts/prepareSU2UpdateTraceCase.py','docs/validation/v147/PLAN.md',
                     'docs/validation/v147/NUMERICAL_PREFLIGHT.md']))
    print(json.dumps(dict(prepared=str(folder),sourceLocalIteration=init['sourceLocalIteration'],newLocalIteration=0)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant',choices=VARIANTS);prepare(p.parse_args().variant)
