"""有界SST受控包的显式比较合同；复用控制初态，不把控制终态当新源场。"""
from pathlib import Path
import copy,json,os
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCase import set_value
FLAVOR='rah-bounded-conservative-v1'
KIND='same-mesh-SST-bounded-conservative-comparison'
STUDY='v131-SST-bounded-conservative-chain-rule'
BASE=ROOT/'outputs/v131_validation';REFERENCE=ROOT/'outputs/v130_validation/column_final'
PLAN=ROOT/'docs/validation/v131/JACOBIAN_PREFLIGHT.md'
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def validate_runtime(folder,data,c,runtime):
    """只允许已定义的两步包兼容/接口和唯一整箭对照，不开放通用默认选项。"""
    assert runtime['flavor']==FLAVOR and runtime['sourceCommit']=='7048fd35f56d4bf52e586429d6746a36a10255f8'
    assert runtime['patchSHA256']=='59b4a400777c60b823ef8d2206a099a1195533a6805cdd81afbd373026501096'
    assert runtime['allCFDObjectsRebuilt'] and runtime['priorThermalAndFrozenHooksRemoved']
    assert runtime['buildReceipt']['exitCode']==0 and not runtime['buildReceipt']['timedOut']
    assert c['SOLVER']=='RANS' and c['SST_OPTIONS']=='V2003m' and c['FLUID_MODEL']=='IDEAL_GAS'
    assert c['MGLEVEL']=='0' and c.get('TIME_DOMAIN','NO')=='NO' and 'MARKER_ISOTHERMAL' in c
    assert 'MARKER_PYTHON_CUSTOM' not in c
    s=data['sst']
    if 'boundedPackageTest' in s:
        t=s['boundedPackageTest'];assert t['studyID']==STUDY and t['variant'] in ('compat','bounded')
        assert c['ITER']=='2' and c['RESTART_SOL']=='NO' and c['OUTPUT_WRT_FREQ']=='1'
        source=ROOT/'outputs/v126_validation/format01'
        expected=config(source/'case.cfg');expected['ITER']='2'
        if t['variant']=='bounded':expected['CONV_NUM_METHOD_TURB']='BOUNDED_SCALAR'
        assert c==expected and sha(folder/'mesh.su2')==sha(source/'mesh.su2')
        assert t['planSHA256']==sha(PLAN)
    else:
        assert data['initialization']['kind']==KIND and 'initializationStudy' not in s
        validate_source(folder,c,data['initialization'])
    return {str(PLAN):sha(PLAN),str(Path(__file__).resolve()):sha(Path(__file__))}

def validate_source(folder,c,init):
    """核对原回执、初场、保存场、BVP来源及唯一配置差异。"""
    assert init['kind']==KIND and init['studyID']==STUDY and init['variant']=='bounded'
    assert init['sameMesh'] and init['samePhysics'] and not init['samePhysicsAndNumerics']
    assert init['localIterationStartsAt']==0 and init['sourceLocalIteration']==149
    assert init['referenceRunLocalIteration']==99 and init['reusedReferenceInitialState']
    assert not init['physicsQualified'] and not init['uninterruptedTrajectoryEquivalent']
    assert Path(init['sourceDirectory']).resolve()==REFERENCE.resolve()
    source=config(REFERENCE/'case.cfg');expected=dict(source);expected['CONV_NUM_METHOD_TURB']='BOUNDED_SCALAR'
    assert c==expected and c['ITER']=='100' and c['OUTPUT_WRT_FREQ']=='50'
    assert init['numericalChanges']=={'CONV_NUM_METHOD_TURB':dict(source='SCALAR_UPWIND',target='BOUNDED_SCALAR')}
    assert init['runtimeChange']['target']==FLAVOR and init['runtimeChange']['source']=='official'
    assert sha(folder/'mesh.su2')==sha(REFERENCE/'mesh.su2')
    assert sha(folder/'initial_flow.csv')==sha(REFERENCE/'initial_flow.csv')==init['initialStateSHA256']
    assert sha(folder/'initial_flow.csv')!=sha(REFERENCE/'restart_flow.csv')
    r=read(REFERENCE/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
    for name in ['case.cfg','mesh.su2','sst_initialization.json','sst_contract.json','case_manifest.json','initial_flow.csv']:
        assert sha(REFERENCE/name)==r['inputSHA256'][str(REFERENCE/name)]
    for name in ['restart_flow_000050.csv','restart_flow_000099.csv','surface_flow.vtk','history.csv']:
        assert sha(REFERENCE/name)==r['outputSHA256'][name]
    assert read(REFERENCE/'sst_state_audit.json')['interfacePassed']
    original=read(REFERENCE/'sst_initialization.json')
    assert init['upstreamInitialization']==original and original['sourceLocalIteration']==149
    for p,h in original['sourceSHA256'].items():assert sha(Path(p))==h
    assert sha(REFERENCE/'column_initialization.json')==original['columnReportSHA256']
    for p,h in init['sourceSHA256'].items():assert sha(Path(p))==h
    assert init['planSHA256']==sha(PLAN)
    compatibility=BASE/'small_comparison.json'
    assert sha(compatibility)==init['compatibilitySHA256'] and read(compatibility)['passed']
    runtime=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'
    assert sha(runtime)==init['runtimeSHA256']
    return dict(init['sourceSHA256'],**{str(PLAN):sha(PLAN),str(runtime):sha(runtime),
        str(compatibility):sha(compatibility),str(Path(__file__).resolve()):sha(Path(__file__))})

def prepare(variant):
    assert variant in ['official','compat','bounded','actual']
    source=REFERENCE if variant=='actual' else ROOT/'outputs/v126_validation/format01'
    folder=BASE/('actual' if variant=='actual' else 'small_'+variant);assert not folder.exists();folder.mkdir()
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    text=set_value(text,'ITER',100 if variant=='actual' else 2)
    if variant in ['actual','bounded']:text=set_value(text,'CONV_NUM_METHOD_TURB','BOUNDED_SCALAR')
    (folder/'case.cfg').write_text(text,encoding='utf-8');os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'))
    m=copy.deepcopy(read(source/'case_manifest.json'));contract=copy.deepcopy(m['sst'])
    contract['configurationSHA256']=sha(folder/'case.cfg')
    if variant in ['actual','bounded']:contract['turbulentConvection']='native bounded scalar; controlled conservative divergence Jacobian /rho'
    m['flight'].update(iterations=100 if variant=='actual' else 2,outputEvery=50 if variant=='actual' else 1,namedCheckpoints=True)
    if variant!='official':m['requiredRuntimeFlavor']=FLAVOR
    if variant=='actual':
        os.link(str(source/'initial_flow.csv'),str(folder/'initial_flow.csv'))
        names=['case.cfg','mesh.su2','sst_initialization.json','sst_contract.json','case_manifest.json','initial_flow.csv',
               'engine_process.json','sst_state_audit.json','column_initialization.json','history.csv','surface_flow.vtk',
               'restart_flow_000050.csv','restart_flow_000099.csv']
        init=dict(kind=KIND,studyID=STUDY,variant='bounded',sourceDirectory=str(source),sourceLocalIteration=149,
            referenceRunLocalIteration=99,reusedReferenceInitialState=True,upstreamInitialization=read(source/'sst_initialization.json'),
            localIterationStartsAt=0,sameMesh=True,samePhysics=True,samePhysicsAndNumerics=False,restartCFLReset=True,
            sourceStateHasNativePrimitive=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False,
            numericalChanges={'CONV_NUM_METHOD_TURB':dict(source='SCALAR_UPWIND',target='BOUNDED_SCALAR')},
            runtimeChange=dict(source='official',target=FLAVOR),initialStateSHA256=sha(folder/'initial_flow.csv'),
            sourceSHA256={str(source/n):sha(source/n) for n in names},planSHA256=sha(PLAN),
            compatibilitySHA256=sha(BASE/'small_comparison.json'),
            runtimeSHA256=sha(ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'))
        validate_source(folder,config(folder/'case.cfg'),init)
        write(folder/'sst_initialization.json',init);m['initialization']=init
        contract.pop('initializationStudy',None)
        contract['initializationSHA256']=sha(folder/'sst_initialization.json')
        contract['numericalStudy']={k:init[k] for k in ['studyID','variant','numericalChanges','runtimeChange']}
    else:contract['boundedPackageTest']=dict(studyID=STUDY,variant=variant,planSHA256=sha(PLAN))
    write(folder/'sst_contract.json',contract)
    m.update(sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'));write(folder/'case_manifest.json',m)
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,
        maximumSeconds=360 if variant=='actual' else 30,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTBoundedStudy.py','docs/validation/v131/JACOBIAN_PREFLIGHT.md'])
    write(folder/'engine_request.json',request)
    print(json.dumps(dict(prepared=str(folder),variant=variant)))
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('variant');prepare(p.parse_args().variant)

