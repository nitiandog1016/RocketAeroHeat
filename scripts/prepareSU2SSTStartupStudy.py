"""同源完整域SST启动对照：只切换流动MUSCL，恢复必须先通过事前筛选。"""
from pathlib import Path
import argparse, copy, json, os
from prepareSU2CheckpointContinuation import ROOT, read, sha, config
from prepareSU2SSTCheckpoint import validate_source_artifacts
from prepareSU2SSTCase import set_value

BASE = ROOT/'outputs/v146_validation'
PLAN = ROOT/'docs/validation/v146/PLAN.md'
KIND = 'same-mesh-SST-reconstruction-startup-comparison'
STUDY = 'v146-common-SST149-first-order-startup'
HISTORY = '(ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)'

def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def original_source():
    return Path(read(ROOT/'outputs/v145_validation/main01/engine_folder.json')['folder']).resolve()

def settings(source, variant):
    """共用预算/输出控制；唯一数值差异是流动一阶与二阶重构。"""
    assert variant in ('control', 'first_order', 'restore_second_order')
    expected = dict(source)
    expected.update(ITER='100', OUTPUT_WRT_FREQ='50', SCREEN_WRT_FREQ_INNER='25', HISTORY_OUTPUT=HISTORY,
                    MUSCL_FLOW='NO' if variant=='first_order' else 'YES')
    return expected

def changes(source, variant):
    target = settings(source, variant)['MUSCL_FLOW']
    return {} if source['MUSCL_FLOW']==target else {'MUSCL_FLOW':dict(source=source['MUSCL_FLOW'], target=target)}

def validate_comparison(source, target, init):
    """禁止混入CFL/物性/下限/网格变化，局部迭代与实际阶数必须一致。"""
    variant = init['variant']
    assert init['kind']==KIND and init['studyID']==STUDY
    assert init['sameMesh'] and init['samePhysics'] and init['localIterationStartsAt']==0
    assert not init['physicsQualified'] and not init['uninterruptedTrajectoryEquivalent']
    assert init['sourceStateHasNativePrimitive'] and init['restartCFLReset']
    assert init['sourceLocalIteration']==(99 if variant=='restore_second_order' else 149)
    assert source['MUSCL_FLOW']==('NO' if variant=='restore_second_order' else 'YES')
    assert target==settings(source, variant)
    assert init['numericalChanges']==changes(source, variant)
    assert init['samePhysicsAndNumerics']==(variant=='control')
    assert init['flowSpatialOrder']==(1 if variant=='first_order' else 2)
    assert source['MUSCL_TURB']=='NO' and source['CONV_NUM_METHOD_TURB']=='SCALAR_UPWIND'
    assert source['CONV_RESIDUAL_MINVAL']=='-8' and source['NUM_METHOD_GRAD']=='WEIGHTED_LEAST_SQUARES'

def validate_source(folder, target, init):
    source = Path(init['sourceDirectory']).resolve()
    assert source==(BASE/'first_order' if init['variant']=='restore_second_order' else original_source())
    validate_comparison(config(source/'case.cfg'), target, init)
    locked = validate_source_artifacts(folder, init)
    assert sha(PLAN)==init['planSHA256']
    if init['variant']=='restore_second_order':
        gate = BASE/'startup_comparison.json'
        assert sha(gate)==init['restorationDecisionSHA256'] and read(gate)['restoreAllowed']
        locked[str(gate)]=sha(gate)
    for path in [PLAN, Path(__file__).resolve()]:locked[str(path)]=sha(path)
    return locked

def prepare(variant):
    """只读共享网格和原生七量，独立写输入清单；拒绝覆盖真实运行。"""
    source = BASE/'first_order' if variant=='restore_second_order' else original_source()
    if variant=='restore_second_order':assert read(BASE/'startup_comparison.json')['restoreAllowed']
    folder = BASE/variant
    assert not folder.exists()
    expected = settings(config(source/'case.cfg'), variant)
    folder.mkdir(parents=True)
    text = (source/'case.cfg').read_text(encoding='utf-8-sig')
    for key, value in expected.items():text=set_value(text,key,value)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    for old,new in [('mesh.su2','mesh.su2'),('restart_flow.csv','initial_flow.csv')]:
        os.link(str(source/old),str(folder/new))
    names=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json',
           'case_manifest.json','restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    init=dict(kind=KIND,studyID=STUDY,variant=variant,numericalChanges=changes(config(source/'case.cfg'),variant),
        flowSpatialOrder=1 if variant=='first_order' else 2,sourceDirectory=str(source),
        sourceLocalIteration=99 if variant=='restore_second_order' else 149,sourceWasDiagnostic=False,
        upstreamInitialization=read(source/'sst_initialization.json'),localIterationStartsAt=0,
        sameMesh=True,samePhysics=True,samePhysicsAndNumerics=variant=='control',sourceStateHasNativePrimitive=True,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(source/n):sha(source/n) for n in names},
        planSHA256=sha(PLAN),startupOnly=variant=='first_order')
    if variant=='restore_second_order':init['restorationDecisionSHA256']=sha(BASE/'startup_comparison.json')
    validate_source(folder,config(folder/'case.cfg'),init)
    write(folder/'sst_initialization.json',init)
    manifest=copy.deepcopy(read(source/'case_manifest.json'));contract=copy.deepcopy(manifest['sst'])
    contract.pop('initializationStudy',None)
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        numericalStudy={k:init[k] for k in ['studyID','variant','numericalChanges','flowSpatialOrder']})
    write(folder/'sst_contract.json',contract)
    manifest['sourceRunControls']={k:manifest['flight'][k] for k in ['iterations','outputEvery','namedCheckpoints']}
    manifest['flight'].update(iterations=100,outputEvery=50,namedCheckpoints=True)
    manifest.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    write(folder/'case_manifest.json',manifest)
    write(folder/'engine_request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        threads=2,maximumSeconds=480,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTStartupStudy.py','docs/validation/v146/PLAN.md']))
    print(json.dumps(dict(prepared=str(folder),sourceLocalIteration=init['sourceLocalIteration'],localFinalIteration=99)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('variant',choices=['control','first_order','restore_second_order'])
    prepare(parser.parse_args().variant)
