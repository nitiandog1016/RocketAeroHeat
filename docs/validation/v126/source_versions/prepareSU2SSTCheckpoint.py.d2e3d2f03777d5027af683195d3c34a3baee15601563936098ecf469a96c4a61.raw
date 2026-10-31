"""复用已通过同阶段接口审核的SST研究场；不把重启动称为已收敛或不中断等价。"""
from pathlib import Path
import json,os,re
from prepareSU2CheckpointContinuation import ROOT,read,sha,config,CONTROL
from su2SSTInterface import state,PHASE_OUTPUT

KIND='same-mesh-SST-checkpoint-research'
OUTPUT_CONTROL=CONTROL|{'WRT_RESTART_COMPACT'}

def settings(text):
    pairs=re.findall(r'^\s*([A-Z0-9_]+)\s*=\s*([^%\r\n]+)',text,re.M)
    assert len(pairs)==len(dict(pairs));return {k:v.strip() for k,v in pairs}


def validate_source(folder,c,init):
    """检查原生终态、阶段审核、具名末场和物理/离散配置；字段缺失不能回退。"""
    assert init['kind']==KIND and init['sameMesh'] and init['samePhysicsAndNumerics']
    assert init['localIterationStartsAt']==0 and not init['physicsQualified'] and not init['uninterruptedTrajectoryEquivalent']
    source=Path(init['sourceDirectory']).resolve();source.relative_to(ROOT/'outputs');assert source!=folder
    r=read(source/'engine_process.json');a=config(source/'case.cfg');audit=read(source/'sst_state_audit.json')
    assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled']
    assert all(r['inputsStable'].values()) and r['runtimeFlavor']=='official' and r['physicsModel']['model']=='SST-2003m'
    assert audit['interfacePassed'] and audit['model']==r['physicsModel'] and not audit['physicsQualified']
    assert audit['heatGradientRelativeError']<2e-6 and audit['EOSRelativeError']<2e-6 and audit['yPlusRelativeError']<2e-6
    assert audit['sourceSHA256'][str(source/'engine_process.json')]==sha(source/'engine_process.json')
    assert {k:v for k,v in c.items() if k not in OUTPUT_CONTROL}=={k:v for k,v in a.items() if k not in OUTPUT_CONTROL},'SST源与目标的物理/离散配置不一致'
    for name in ['mesh.su2','case.cfg','sst_contract.json']:
        assert sha(source/name)==r['inputSHA256'][str(source/name)]
    assert sha(folder/'mesh.su2')==sha(source/'mesh.su2')
    for name in ['restart_flow.csv','surface_flow.vtk','history.csv']:
        assert sha(source/name)==r['outputSHA256'][name]
    step=init['sourceLocalIteration'];assert step==audit['checkpointLocalIteration']
    for stem,ext in [('restart_flow','csv'),('surface_flow','vtk')]:
        name=stem+'_%06d.'%step+ext
        assert sha(source/name)==r['outputSHA256'][name]==r['outputSHA256'][stem+'.'+ext]
    assert sha(folder/'initial_flow.csv')==sha(source/'restart_flow.csv')==init['initialStateSHA256']
    extended=r['physicsModel'].get('outputStateContract')==PHASE_OUTPUT
    assert init['sourceStateHasNativePrimitive']==extended
    state(folder/'initial_flow.csv',True,extended)
    for p,h in init['sourceSHA256'].items():assert sha(Path(p))==h
    return dict(init['sourceSHA256'],**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())})


def prepare(folder,source,text,contract,manifest):
    """目标目录已有新网格/参考配置；只读共享源终态，写独立的新输入清单。"""
    from prepareSU2SSTCase import set_value
    folder=Path(folder).resolve();source=Path(source).resolve()
    folder.relative_to(ROOT/'outputs');source.relative_to(ROOT/'outputs')
    assert folder!=source and not (folder/'engine_process.json').exists() and not (folder/'initial_flow.csv').exists()
    r=read(source/'engine_process.json');audit=read(source/'sst_state_audit.json')
    for key in ['intensity','viscosityRatio','Prt']:assert contract[key]==r['physicsModel'][key]
    for key in ['initialization','diagnosticOnly','phaseEvidence']:contract.pop(key,None)
    init=dict(kind=KIND,sourceDirectory=str(source),sourceLocalIteration=audit['checkpointLocalIteration'],
        sourceWasDiagnostic=bool(r['physicsModel'].get('diagnosticOnly',False)),localIterationStartsAt=0,
        sameMesh=True,samePhysicsAndNumerics=True,initialStateSHA256=sha(source/'restart_flow.csv'),
        sourceStateHasNativePrimitive=r['physicsModel'].get('outputStateContract')==PHASE_OUTPUT,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False)
    locks=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json','case_manifest.json',
        'restart_flow.csv','surface_flow.vtk','history.csv']
    if (source/'sst_initialization.json').exists():
        locks.append('sst_initialization.json');prior=read(source/'sst_initialization.json')
        init['upstreamInitialization']=prior
    init['sourceSHA256']={str(source/n):sha(source/n) for n in locks}
    for key,value in dict(RESTART_SOL='YES',SOLUTION_FILENAME='initial_flow',READ_BINARY_RESTART='NO').items():text=set_value(text,key,value)
    os.link(str(source/'restart_flow.csv'),str(folder/'initial_flow.csv'))
    validate_source(folder,settings(text),init)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n',encoding='utf-8')
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        meshSHA256=sha(folder/'mesh.su2'))
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    manifest.update(sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'),initialization=init)
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return contract
