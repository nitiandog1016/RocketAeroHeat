"""真实299步源场的继续输入和拒绝反例；不运行引擎、不改原文件。"""
from pathlib import Path
import json,os,sys,uuid
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import prepare,sha,read


def main():
    out=ROOT/'outputs/v119_validation';out.mkdir(exist_ok=True)
    source=ROOT/'outputs/v118_validation/actual01/budget_run'
    target=out/'continuation01';assert not target.exists()
    def setup(folder):
        folder.mkdir()
        for name in ['mesh.su2','projection.mat']:os.link(str(source/name),str(folder/name))
        (folder/'case.cfg').write_bytes((source/'case.cfg').read_bytes())
        m=read(source/'case_manifest.json');m.pop('budgetRetry',None)
        (folder/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    setup(target)
    # 在源场仍锁定的情况下，准备器验证所有输入并复用精确原生守恒场。
    lineage=prepare(source,target)
    assert lineage['completedUpdatesBeforeRun']==300 and lineage['sourceSavedLocalIteration']==299
    assert sha(target/'initial_flow.csv')==sha(source/'restart_flow_000299.csv')
    assert (target/'initial_flow.csv').stat().st_ino==(source/'restart_flow.csv').stat().st_ino
    rejected=[]
    for mode in ['wrong_mach','wrong_mesh','existing_initial']:
        folder=out/mode;setup(folder)
        if mode=='wrong_mach':
            p=folder/'case.cfg';p.write_text(p.read_text().replace('MACH_NUMBER= 5','MACH_NUMBER= 4'),encoding='utf-8')
        if mode=='wrong_mesh':
            # 不改共享原网格；创建不同名字的小网格并交换目标名称。
            (folder/'mesh.su2').rename(folder/'original_mesh.su2')
            (folder/'mesh.su2').write_text('NDIME= 2\n')
        if mode=='existing_initial':(folder/'initial_flow.csv').write_text('preexisting; must not overwrite\n')
        try:prepare(source,folder)
        except AssertionError:rejected.append(mode)
        else:raise AssertionError(mode+' unexpectedly accepted')
    assert len(rejected)==3
    request=dict(workDirectory=str(target),configFile='case.cfg',meshFile='mesh.su2',threads=2,
                 maximumSeconds=720,dryRun=False,jobToken=str(uuid.uuid4()),sourceFiles=[
                     str(target/'continuation.json'),str(target/'case_manifest.json'),str(target/'projection.mat'),
                     str(source/'engine_process.json'),'scripts/prepareSU2CheckpointContinuation.py',
                     'tests/testSU2ContinuationPreparation.py','docs/validation/v119/PLAN.md'])
    (target/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    report=dict(passed=True,rejected=rejected,sourceAndInitialExact=True,inputIsHardlink=True,
                sourceLocalIteration=299,completedUpdatesBeforeRun=300,physicsQualified=False,
                inputSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/prepareSU2CheckpointContinuation.py',source/'engine_process.json',target/'continuation.json']})
    (out/'preparation_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
