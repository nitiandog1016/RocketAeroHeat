"""一次SST阶段诊断：保留失败来源身份，不伪装生产续算或层流模型转换。"""
from pathlib import Path
import argparse,json,os
import numpy as np
from prepareSU2CheckpointContinuation import ROOT,config,read,sha,CONTROL
from prepareSU2SSTCase import set_value
from su2SSTInterface import state

KIND='same-mesh-SST-single-step-phase-diagnostic'


def validate_settings(c,old,init):
    """只放行一个同物理更新；输出增强与迭代控制以外的变更拒绝。"""
    assert init['kind']==KIND and init['diagnosticOnly'] is True
    assert init['sameMesh'] and init['samePhysicsAndNumerics']
    assert init['localIterationStartsAt']==0 and init['requestedUpdates']==1
    assert init['sourceNativeSucceeded'] and init['sourcePostAuditPassed'] is False
    assert init['physicsQualified'] is False
    assert c['ITER']=='1' and c['OUTPUT_WRT_FREQ']=='1'
    assert c['OUTPUT_FILES']=='(SURFACE_PARAVIEW_ASCII, RESTART_ASCII, CSV)'
    assert c['VOLUME_OUTPUT']=='(COORDINATES, SOLUTION, PRIMITIVE)'
    assert c['CONV_RESIDUAL_MINVAL']=='-8'
    assert {k:v for k,v in c.items() if k not in CONTROL}=={k:v for k,v in old.items() if k not in CONTROL}


def validate_source(init,folder,c):
    """审核实际失败的上游回执、具名检查点及每个输入；不得偷换源身份。"""
    source=Path(init['sourceDirectory']).resolve();source.relative_to(ROOT/'outputs')
    assert source!=folder
    r=read(source/'engine_process.json');p=read(source/'workflow/pipeline_process.json')
    assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled']
    assert all(r['inputsStable'].values()) and r['runtimeFlavor']=='official'
    assert r['physicsModel']['model']=='SST-2003m'
    assert p['status']=='terminal' and p['exitCode']==1
    validate_settings(c,config(source/'case.cfg'),init)
    for name in ('case.cfg','mesh.su2','case_manifest.json','sst_contract.json'):
        assert sha(source/name)==r['inputSHA256'][str(source/name)]
    step=init['sourceLocalIteration'];assert isinstance(step,int) and step>=0
    checkpoint=source/('restart_flow_%06d.csv'%step)
    assert sha(checkpoint)==r['outputSHA256'][checkpoint.name]==r['outputSHA256']['restart_flow.csv']
    assert sha(folder/'initial_flow.csv')==sha(checkpoint)==init['initialStateSHA256']
    assert sha(folder/'mesh.su2')==sha(source/'mesh.su2')
    hist=np.atleast_1d(np.genfromtxt(str(source/'history.csv'),names=True,delimiter=','))
    assert sha(source/'history.csv')==r['outputSHA256']['history.csv']
    assert np.array_equal(hist['Inner_Iter'],np.arange(len(hist))) and step==int(hist['Inner_Iter'][-1])
    return {str(Path(__file__).resolve()):sha(Path(__file__).resolve())}


def prepare(source,target):
    """仅创建新目录，历史网格/输入作只读共享；重启动输出使用另一基名。"""
    source=Path(source).resolve();target=Path(target).resolve()
    source.relative_to(ROOT/'outputs');target.relative_to(ROOT/'outputs')
    assert not target.exists()
    r=read(source/'engine_process.json');contract=read(source/'sst_contract.json')
    assert r['physicsModel']==contract
    a,keys,nd,rho,u,e=state(source/'restart_flow.csv',True)
    hist=np.atleast_1d(np.genfromtxt(str(source/'history.csv'),names=True,delimiter=','))
    step=int(hist['Inner_Iter'][-1]);checkpoint=source/('restart_flow_%06d.csv'%step)
    assert sha(checkpoint)==sha(source/'restart_flow.csv')
    lineage=dict(kind=KIND,diagnosticOnly=True,sourceDirectory=str(source),sourceLocalIteration=step,
        localIterationStartsAt=0,requestedUpdates=1,sameMesh=True,samePhysicsAndNumerics=True,
        sourceNativeSucceeded=True,sourcePostAuditPassed=False,physicsQualified=False,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,initialStateSHA256=sha(checkpoint))
    locks=[source/n for n in ['engine_process.json','workflow/pipeline_process.json','case.cfg','case_manifest.json',
        'mesh.su2','sst_contract.json','sst_initialization.json','history.csv',checkpoint.name]]
    lineage['sourceSHA256']={str(p):sha(p) for p in locks}
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    for k,v in dict(ITER=1,OUTPUT_WRT_FREQ=1,SCREEN_WRT_FREQ_INNER=1,
        OUTPUT_FILES='(SURFACE_PARAVIEW_ASCII, RESTART_ASCII, CSV)').items():text=set_value(text,k,v)
    target.mkdir(parents=True)
    os.link(str(source/'mesh.su2'),str(target/'mesh.su2'));os.link(str(checkpoint),str(target/'initial_flow.csv'))
    (target/'case.cfg').write_text(text,encoding='utf-8')
    validate_source(lineage,target,config(target/'case.cfg'))
    (target/'sst_initialization.json').write_text(json.dumps(lineage,indent=2)+'\n',encoding='utf-8')
    contract.update(initializationSHA256=sha(target/'sst_initialization.json'),configurationSHA256=sha(target/'case.cfg'),
        phaseEvidence='native-volume-CSV; flow primitive before current turbulence update',diagnosticOnly=True)
    (target/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    manifest=read(source/'case_manifest.json');manifest['sst']=contract;manifest['initialization']=lineage
    manifest['sstContractSHA256']=sha(target/'sst_contract.json')
    manifest['flight'].update(iterations=1,outputEvery=1)
    (target/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    request=dict(workDirectory=str(target),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=180,
        dryRun=False,sourceFiles=['scripts/prepareSU2SSTPhaseProbe.py','docs/validation/v125/PLAN.md'])
    (target/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    return dict(directory=str(target),nodes=len(a),sourceIteration=step,updates=1,maximumSeconds=180)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--target',required=True);a=p.parse_args()
    print(json.dumps(prepare(a.source,a.target)))
