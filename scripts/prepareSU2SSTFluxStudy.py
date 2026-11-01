"""准备固定官方通量的小例与整箭同源对照；不混入更新保护/物理模型改动。"""
from pathlib import Path
import argparse,copy,json,os
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import validate_source_artifacts
from prepareSU2SSTCase import set_value

BASE=ROOT/'outputs/v148_validation';PLAN=ROOT/'docs/validation/v148/PLAN.md'
KIND='same-mesh-SST-flow-flux-comparison';STUDY='v148-common-SST99-AUSM-vs-AUSMPLUSUP2'
VARIANTS=('small_up2','control','up2')
def write(p,v):p.write_bytes((json.dumps(v,indent=2,allow_nan=False)+'\n').encode('utf-8'))
def source_path(variant):
    assert variant in VARIANTS
    return ROOT/('outputs/v145_validation/fixture02/case' if variant=='small_up2' else 'outputs/v146_validation/control')
def settings(source,variant):
    c=dict(source);small=variant=='small_up2'
    c.update(ITER='1' if small else '50',OUTPUT_WRT_FREQ='1' if small else '25',SCREEN_WRT_FREQ_INNER='1' if small else '25',
        HISTORY_OUTPUT='(ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)',
        CONV_NUM_METHOD_FLOW='AUSM' if variant=='control' else 'AUSMPLUSUP2')
    return c
def changes(variant):
    return {} if variant=='control' else {'CONV_NUM_METHOD_FLOW':dict(source='AUSM',target='AUSMPLUSUP2')}
def validate_source(folder,c,init):
    """逐项固定源、唯一通量变化及前置证据，拒绝默认方案静默变化。"""
    variant=init['variant'];source=source_path(variant);src=config(source/'case.cfg')
    assert Path(init['sourceDirectory']).resolve()==source.resolve()
    assert init['kind']==KIND and init['studyID']==STUDY and init['localIterationStartsAt']==0
    assert init['sameMesh'] and init['samePhysics'] and not init['physicsQualified'] and not init['uninterruptedTrajectoryEquivalent']
    assert init['sourceLocalIteration']==(1 if variant=='small_up2' else 99)
    assert init['samePhysicsAndNumerics']==(variant=='control') and init['numericalChanges']==changes(variant)
    assert src['CONV_NUM_METHOD_FLOW']=='AUSM' and src['MUSCL_FLOW']=='YES'
    assert c==settings(src,variant) and init['flowFlux']==c['CONV_NUM_METHOD_FLOW']
    assert c.get('USE_ACCURATE_FLUX_JACOBIANS','NO')=='NO' and init['accurateFluxJacobianEnabled'] is False
    locks=validate_source_artifacts(folder,init)
    for p in [PLAN,Path(__file__).resolve(),BASE/'flux_preflight02/report.json']:
        locks[str(p)]=sha(p)
    assert sha(PLAN)==init['planSHA256'] and sha(BASE/'flux_preflight02/report.json')==init['fluxPreflightSHA256']
    assert read(BASE/'flux_preflight02/report.json')['passed']
    if variant!='small_up2':
        p=BASE/'small_up2/sst_state_audit.json';r=read(BASE/'small_up2/engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        a=read(p);assert a['interfacePassed'] and a['model']['numericalStudy']['flowFlux']=='AUSMPLUSUP2'
        assert a['sourceSHA256'][str(BASE/'small_up2/engine_process.json')]==sha(BASE/'small_up2/engine_process.json')
        assert sha(p)==init['smallInterfaceSHA256'];locks[str(p)]=sha(p)
    return locks
def prepare(variant):
    source=source_path(variant);folder=BASE/variant;assert not folder.exists()
    assert read(BASE/'flux_preflight02/report.json')['passed']
    if variant!='small_up2':assert read(BASE/'small_up2/sst_state_audit.json')['interfacePassed']
    folder.mkdir(parents=True);c=settings(config(source/'case.cfg'),variant)
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    for k,v in c.items():text=set_value(text,k,v)
    (folder/'case.cfg').write_bytes(text.encode('utf-8'))
    for a,b in [('mesh.su2','mesh.su2'),('restart_flow.csv','initial_flow.csv')]:os.link(str(source/a),str(folder/b))
    names=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json','case_manifest.json',
        'restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    init=dict(kind=KIND,studyID=STUDY,variant=variant,sourceDirectory=str(source),sourceLocalIteration=1 if variant=='small_up2' else 99,
        sourceWasDiagnostic=False,upstreamInitialization=read(source/'sst_initialization.json'),localIterationStartsAt=0,
        sameMesh=True,samePhysics=True,samePhysicsAndNumerics=variant=='control',sourceStateHasNativePrimitive=True,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False,numericalChanges=changes(variant),
        flowFlux=c['CONV_NUM_METHOD_FLOW'],accurateFluxJacobianEnabled=False,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(source/n):sha(source/n) for n in names},
        planSHA256=sha(PLAN),fluxPreflightSHA256=sha(BASE/'flux_preflight02/report.json'))
    if variant!='small_up2':init['smallInterfaceSHA256']=sha(BASE/'small_up2/sst_state_audit.json')
    validate_source(folder,c,init);write(folder/'sst_initialization.json',init)
    manifest=copy.deepcopy(read(source/'case_manifest.json'));contract=copy.deepcopy(manifest['sst'])
    for key in ['initializationStudy','numericalStudy','updateTraceStudy']:contract.pop(key,None)
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        numericalStudy={k:init[k] for k in ['studyID','variant','numericalChanges','flowFlux','accurateFluxJacobianEnabled']})
    write(folder/'sst_contract.json',contract);manifest['flight'].update(iterations=int(c['ITER']),outputEvery=int(c['OUTPUT_WRT_FREQ']),namedCheckpoints=True)
    manifest.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'),requiredRuntimeFlavor='official')
    write(folder/'case_manifest.json',manifest)
    write(folder/'engine_request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,
        maximumSeconds=60 if variant=='small_up2' else 240,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTFluxStudy.py','docs/validation/v148/PLAN.md']))
    print(json.dumps(dict(prepared=str(folder),sourceLocalIteration=init['sourceLocalIteration'],updates=int(c['ITER']),flowFlux=init['flowFlux'])))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant',choices=VARIANTS);prepare(p.parse_args().variant)
