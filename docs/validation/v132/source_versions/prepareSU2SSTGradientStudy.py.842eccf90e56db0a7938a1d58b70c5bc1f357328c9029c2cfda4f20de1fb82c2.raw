"""v132官方梯度的有限原生对照：同终态、同重构、完整原域，不改默认。"""
from pathlib import Path
import argparse,copy,json,os
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import validate_source_artifacts
from prepareSU2SSTCase import set_value

BASE=ROOT/'outputs/v132_validation'
SOURCE=ROOT/'outputs/v130_validation/column_final'
PLAN=ROOT/'docs/validation/v132/NATIVE_PREFLIGHT.md'
KIND='same-mesh-SST-gradient-comparison'
STUDY='v132-tail-rim-viscous-gradient'

def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def changes(variant):
    assert variant in ('wls','gg')
    return {} if variant=='wls' else {'NUM_METHOD_GRAD':dict(source='WEIGHTED_LEAST_SQUARES',target='GREEN_GAUSS')}

def settings(source,variant):
    """两组显式固定原先继承的WLS重构，唯一算法差异为黏性/源项梯度。"""
    result=dict(source)
    result.update(ITER='10',OUTPUT_WRT_FREQ='10',SCREEN_WRT_FREQ_INNER='1',
                  NUM_METHOD_GRAD_RECON='WEIGHTED_LEAST_SQUARES')
    for key,value in changes(variant).items():result[key]=value['target']
    return result

def validate_comparison(source,target,init):
    assert init['kind']==KIND and init['studyID']==STUDY
    assert init['sameMesh'] and init['samePhysics'] and init['localIterationStartsAt']==0
    assert init['sourceLocalIteration']==99 and not init['physicsQualified']
    assert not init['uninterruptedTrajectoryEquivalent'] and init['sourceStateHasNativePrimitive']
    assert source['NUM_METHOD_GRAD']=='WEIGHTED_LEAST_SQUARES' and 'NUM_METHOD_GRAD_RECON' not in source
    assert init['numericalChanges']==changes(init['variant'])
    assert init['samePhysicsAndNumerics']==(init['variant']=='wls')
    assert target==settings(source,init['variant'])
    assert target['CONV_NUM_METHOD_TURB']=='SCALAR_UPWIND' and target['CONV_RESIDUAL_MINVAL']=='-8'

def validate_source(folder,c,init):
    assert Path(init['sourceDirectory']).resolve()==SOURCE.resolve()
    validate_comparison(config(SOURCE/'case.cfg'),c,init)
    locked=validate_source_artifacts(folder,init)
    assert sha(PLAN)==init['planSHA256'] and sha(BASE/'local_diffusion.json')==init['localReviewSHA256']
    locked.update({str(p):sha(p) for p in [PLAN,BASE/'local_diffusion.json',Path(__file__).resolve()]})
    return locked

def prepare(variant):
    folder=BASE/variant;assert not folder.exists();folder.mkdir()
    expected=settings(config(SOURCE/'case.cfg'),variant)
    text=(SOURCE/'case.cfg').read_text(encoding='utf-8-sig')
    for key,value in expected.items():text=set_value(text,key,value)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    for old,new in [('mesh.su2','mesh.su2'),('restart_flow.csv','initial_flow.csv')]:
        os.link(str(SOURCE/old),str(folder/new))
    locks=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json',
           'case_manifest.json','restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    init=dict(kind=KIND,studyID=STUDY,variant=variant,numericalChanges=changes(variant),
        sourceDirectory=str(SOURCE),sourceLocalIteration=99,sourceWasDiagnostic=False,
        upstreamInitialization=read(SOURCE/'sst_initialization.json'),localIterationStartsAt=0,
        sameMesh=True,samePhysics=True,samePhysicsAndNumerics=variant=='wls',
        sourceStateHasNativePrimitive=True,restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(SOURCE/n):sha(SOURCE/n) for n in locks},
        physicsQualified=False,planSHA256=sha(PLAN),localReviewSHA256=sha(BASE/'local_diffusion.json'),
        reconstructionGradient='WEIGHTED_LEAST_SQUARES',gradientScope='flow-and-turbulence-viscous-and-source')
    validate_source(folder,config(folder/'case.cfg'),init)
    write(folder/'sst_initialization.json',init)
    m=copy.deepcopy(read(SOURCE/'case_manifest.json'));contract=copy.deepcopy(m['sst'])
    contract.pop('initializationStudy',None)
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
                    numericalStudy={k:init[k] for k in ['studyID','variant','numericalChanges']})
    write(folder/'sst_contract.json',contract)
    m['sourceRunControls']={k:m['flight'][k] for k in ['iterations','outputEvery','namedCheckpoints']}
    m['flight'].update(iterations=10,outputEvery=10,namedCheckpoints=True)
    m.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    write(folder/'case_manifest.json',m)
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=180,
        dryRun=False,sourceFiles=['scripts/prepareSU2SSTGradientStudy.py','docs/validation/v132/NATIVE_PREFLIGHT.md'])
    write(folder/'engine_request.json',request)
    print(json.dumps(dict(prepared=str(folder),sourceIteration=99,localFinalIteration=9)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant',choices=['wls','gg']);prepare(p.parse_args().variant)
