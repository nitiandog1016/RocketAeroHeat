"""完整域SST的显式同源推进对照；仅允许事前声明的CFL上限改变。"""
from pathlib import Path
import argparse,json,os,copy
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import OUTPUT_CONTROL,validate_source_artifacts
from prepareSU2SSTCase import set_value
from su2SSTInterface import PHASE_OUTPUT

KIND='same-mesh-SST-CFL-cap-comparison'
STUDY='v127-common-SST49-CFL5-vs-0.5'
HISTORY='(ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)'


def validate_comparison(source,target,init):
    """除输出控制外只允许一个具名参数变化；控制、候选标签与实际配置必须一致。"""
    assert init['kind']==KIND and init['studyID']==STUDY and init['variant'] in ['control','cap05']
    assert init['sameMesh'] and init['samePhysics'] and init['localIterationStartsAt']==0
    assert not init['physicsQualified'] and not init['uninterruptedTrajectoryEquivalent']
    assert init['sourceLocalIteration']==49 and init['sourceStateHasNativePrimitive']
    assert source['CFL_ADAPT']=='YES' and source['CFL_NUMBER']=='0.1'
    assert source['CFL_ADAPT_PARAM']=='(0.5, 1.1, 0.05, 5.0)'
    expected=dict(source);changes={}
    if init['variant']=='cap05':
        changes={'CFL_ADAPT_PARAM':dict(source=source['CFL_ADAPT_PARAM'],target='(0.5, 1.1, 0.05, 0.5)')}
        expected['CFL_ADAPT_PARAM']=changes['CFL_ADAPT_PARAM']['target']
    assert init['numericalChanges']==changes and init['samePhysicsAndNumerics']==(not changes)
    assert {k:v for k,v in target.items() if k not in OUTPUT_CONTROL}=={k:v for k,v in expected.items() if k not in OUTPUT_CONTROL}
    assert target['ITER']=='150' and target['OUTPUT_WRT_FREQ']=='75' and target['HISTORY_OUTPUT']==HISTORY
    assert target['WRT_RESTART_COMPACT']=='NO' and target['VOLUME_OUTPUT']=='(COORDINATES, SOLUTION, PRIMITIVE)'
    assert target['RESTART_SOL']=='YES' and target['SOLUTION_FILENAME']=='initial_flow' and target['READ_BINARY_RESTART']=='NO'


def validate_source(folder,c,init):
    source=Path(init['sourceDirectory']).resolve();source.relative_to(ROOT/'outputs')
    validate_comparison(config(source/'case.cfg'),c,init)
    locked=validate_source_artifacts(folder,init)
    locked[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
    assert sha(ROOT/'docs/validation/v127/PLAN.md')==init['planSHA256']
    locked[str(ROOT/'docs/validation/v127/PLAN.md')]=init['planSHA256']
    return locked


def current_manifest(source_manifest,contract,init,contract_hash):
    """源工况包含旧输出控制时单列保留；新任务清单的步数/保存频率必须与新配置一致。"""
    m=copy.deepcopy(source_manifest)
    m['sourceRunControls']={k:m['flight'][k] for k in ['iterations','outputEvery','namedCheckpoints']}
    m['flight'].update(iterations=150,outputEvery=75,namedCheckpoints=True)
    m.update(initialization=init,sst=contract,sstContractSHA256=contract_hash)
    return m


def prepare(variant):
    assert variant in ['control','cap05']
    source=Path(read(ROOT/'outputs/v126_validation/main01/engine_folder.json')['folder']).resolve()
    folder=ROOT/'outputs/v127_validation'/variant;assert not folder.exists();folder.mkdir(parents=True)
    r=read(source/'engine_process.json');audit=read(source/'sst_state_audit.json')
    assert audit['checkpointLocalIteration']==49 and r['physicsModel']['outputStateContract']==PHASE_OUTPUT
    text=(source/'case.cfg').read_text(encoding='utf-8-sig');changes={}
    if variant=='cap05':changes={'CFL_ADAPT_PARAM':dict(source=config(source/'case.cfg')['CFL_ADAPT_PARAM'],target='(0.5, 1.1, 0.05, 0.5)')}
    updates=dict(ITER=150,OUTPUT_WRT_FREQ=75,SCREEN_WRT_FREQ_INNER=25,HISTORY_OUTPUT=HISTORY)
    for key,values in changes.items():updates[key]=values['target']
    for key,value in updates.items():text=set_value(text,key,value)
    for name in ['mesh.su2','restart_flow.csv']:
        os.link(str(source/name),str(folder/('initial_flow.csv' if name=='restart_flow.csv' else name)))
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    locks=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json',
           'case_manifest.json','restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    init=dict(kind=KIND,studyID=STUDY,variant=variant,numericalChanges=changes,
        sourceDirectory=str(source),sourceLocalIteration=49,sourceWasDiagnostic=False,
        upstreamInitialization=read(source/'sst_initialization.json'),localIterationStartsAt=0,
        sameMesh=True,samePhysics=True,samePhysicsAndNumerics=not changes,
        sourceStateHasNativePrimitive=True,restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(source/n):sha(source/n) for n in locks},
        physicsQualified=False,planSHA256=sha(ROOT/'docs/validation/v127/PLAN.md'))
    validate_source(folder,config(folder/'case.cfg'),init)
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n',encoding='utf-8')
    contract=copy.deepcopy(r['physicsModel'])
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        numericalStudy=dict(studyID=STUDY,variant=variant,numericalChanges=changes))
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    m=current_manifest(read(source/'case_manifest.json'),contract,init,sha(folder/'sst_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=600,
        dryRun=False,sourceFiles=['scripts/prepareSU2SSTNumericalStudy.py','docs/validation/v127/PLAN.md'])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(folder=str(folder),variant=variant,updates=150)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant',choices=['control','cap05']);prepare(p.parse_args().variant)
