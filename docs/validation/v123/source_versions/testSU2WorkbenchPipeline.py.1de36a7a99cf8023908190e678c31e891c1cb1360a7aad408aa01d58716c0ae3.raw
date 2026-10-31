"""后台监督单元测试：真实短子进程取消/超时，合成阶段不冒充数值求解。"""
from pathlib import Path
import json, os, sys, threading
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runSU2WorkbenchPipeline import Pipeline
from runEmbeddedSU2 import atomic_write_json
def main():
    out=ROOT/'outputs/v123_validation/supervision_contracts';out.mkdir(exist_ok=False)
    checks=[]
    for mode in ['timeout','cancel','wrongToken']:
        folder=out/mode;folder.mkdir();request=folder/'engine_request.json'
        request.write_text(json.dumps(dict(workDirectory=str(folder),jobToken='owned-token',maximumSeconds=1,
            propertyInitialization=dict(model='HEOS-Air',runtimeFlavor='rah-coolprop-recon-v3',sourceDirectory='unused'))))
        p=Pipeline(request);timer=None
        def prepare(*args):
            nonlocal timer
            if mode in ['cancel','wrongToken']:
                token='owned-token' if mode=='cancel' else 'foreign-token'
                timer=threading.Timer(.3,lambda:atomic_write_json(folder/'engine_cancel.json',{'jobToken':token}));timer.start()
            seconds=.4 if mode=='timeout' else 3
            duration=20 if mode in ['timeout','cancel'] else .6
            p.run_property(p.out,'synthetic_child',[sys.executable,'-c',f'import time;time.sleep({duration})'],ROOT,dict(os.environ),seconds)
        # 只有独立监督合同：合成引擎/审核为空；报告明确没有SU2数值运行。
        with patch('prepareSU2CoolPropCase.prepare',prepare),patch.object(p,'run_engine',lambda:None),patch('auditSU2CoolPropState.audit',lambda *a:None),patch('runSU2WorkbenchPipeline.sha',side_effect=lambda path:'0'*64 if Path(path).name=='audit.json' else __import__('runEmbeddedSU2').sha(path)):
            code=p.execute()
        if timer:timer.join()
        expected={'timeout':124,'cancel':125,'wrongToken':0}[mode];assert code==expected
        q=json.loads((p.out/'synthetic_child_process.json').read_text());assert q['status']=='terminal'
        assert not (folder/'engine_process.json').exists() and not p.state['engineStarted']
        if mode=='wrongToken':assert q['exitCode']==0 and not q['cancelled'] and not q['timedOut']
        else:assert q['exitCode']!=0 and (q['cancelled'] if mode=='cancel' else q['timedOut'])
        checks.append(dict(mode=mode,pipelineExit=code,childExit=q['exitCode'],cancelled=q['cancelled'],timedOut=q['timedOut']))
    report=dict(passed=True,checks=checks,realSyntheticChildren=3,actualSU2Runs=0,numericalSolverClaimed=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
