# -*- coding: utf-8 -*-
"""生成唯一流向/法向加密控制，复用准确气体，只准备网格与初场。"""
from pathlib import Path
import argparse
import json
import shutil
import time
from prepareHybridEulerCase import prepare
from hybridEulerEvolution import audit_evolution_inputs
from runEilmerProcess import run,sha


def execute(case,record,baseline,coarse_audit):
    start=time.monotonic();case,record,baseline,coarse_audit=map(Path,(case,record,baseline,coarse_audit))
    c=json.loads(coarse_audit.read_text())
    assert c['numericalConvergencePassed'] and not c['planAcceptancePassed']
    assert not case.exists() and not record.exists()
    design=prepare(case,(120,2,40));record.mkdir(parents=True)
    shutil.copyfile(baseline/'rah-air.gas',case/'rah-air.gas')
    job=(case/'job.lua').read_text().replace('config.max_step=400','config.max_step=2400')
    job=job.replace('config.write_flow_solution_at_step=200','config.write_flow_solution_at_step=1200')
    job=job.replace('config.dt_plot_value=1.0','config.dt_plot_value=0.0007')
    (case/'job.lua').write_text(job,encoding='utf8')
    design['nativeConfig'].update(max_step=2400,write_flow_solution_at_step=1200)
    design.update(candidate='refined-120-2-40',executionAuthorized=False)
    design['filesSHA256']={n:sha(case/n) for n in design['filesSHA256']}
    (case/'manifest.json').write_text(json.dumps(design,indent=2)+'\n')
    sources=[Path(__file__),Path(__file__).with_name('prepareHybridEulerCase.py'),Path(__file__).with_name('hybridEulerEvolution.py'),
             baseline/'rah-air.gas',baseline/'lmrsim/config',coarse_audit]
    for p in sources[:3]:shutil.copyfile(p,record/p.name)
    locks=sources+[p for p in case.iterdir() if p.is_file()];commands=[]
    for i,cmd in enumerate([['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]):
        receipt=run(cmd,case,record/('native%02d'%i),60-(time.monotonic()-start),64,locks,4096,4096)
        commands.append(receipt);assert receipt['passed']
    # 用明确模板构造控制期望，不由读回结果自证本应采用的终止/保存值。
    ctl=json.loads((baseline/'lmrsim/control').read_text());ctl.update(max_step=2400,write_flow_solution_at_step=1200)
    result=dict(schema='RAH.HybridEulerEvolutionPreparation.v1',profile='refined-120-2-40',passed=True,
        case=str(case.resolve()),geometry=design['geometry'],expectedControl=ctl,
        configOutputOverrides={'dt_plot_schedule_values':[.0007]},
        baselineConfig=str((baseline/'lmrsim/config').resolve()),
        sourceSHA256={str(p.resolve()):sha(p) for p in sources[3:]},
        nativeInputSHA256={str(p.resolve()):sha(p) for p in case.rglob('*') if p.is_file()},
        sourceStep=0,sourceSnapshot=0,physicalQualification=False,commands=commands)
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    audit,grid,mapping=audit_evolution_inputs(case,record/'report.json')
    (record/'geometry_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(passed=True,elapsedSeconds=time.monotonic()-start,cells=audit['cells'],ghosts=audit['ghostCells'])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','record','baseline','coarse-audit']:p.add_argument('--'+k,required=True)
    a=p.parse_args();execute(a.case,a.record,a.baseline,a.coarse_audit)
