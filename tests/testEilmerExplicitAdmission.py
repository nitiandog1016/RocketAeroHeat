"""完整域唯一候选的来源/目录/资源拒绝检查；不启动CFD或重复小例。"""
from pathlib import Path
import argparse,copy,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from eilmerExplicitAdmission import audit_admission,sha
from eilmerExplicitProfile import settings
from readEilmerNative import audit_run_scope
from auditEilmerCompletion import audit_completion
from prepareEilmerRocketCase import make_job,validate_condition


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'outputs/v184_validation';c=json.loads((base/'condition.json').read_text())
    m=json.loads((ROOT/'outputs/v178_validation/real01/case_manifest.json').read_text());m['condition']=c
    a,_=audit_admission(m,base/'real01')
    original=json.loads((ROOT/'outputs/v183_validation/small01/case_manifest.json').read_text())
    old=settings(original['condition']['explicit']);new=settings(c['explicit'])
    changed={k for k in old if old[k]!=new[k]}
    assert changed=={'max_step','write_loads_at_step','write_flow_solution_at_step'}
    assert 'NewtonKrylovGlobalConfig' not in make_job(m['mesh']['nativeSourcePath'],c)
    rejected=[]
    def reject(label,fn):
        try:fn()
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    reject('closed-full-domain-profile',lambda:audit_run_scope(m,base/'real01'))
    reject('unregistered-directory',lambda:audit_admission(m,base/'other'))
    for label,mutate in [('changed-Mach',lambda v:v['condition'].update(Mach=6)),
                         ('changed-mesh',lambda v:v['mesh'].update(sourceSHA256='0'*64)),
                         ('changed-wall-count',lambda v:v['mesh'].update(wallFaces=19751)),
                         ('changed-area',lambda v:v['mesh'].update(wallAreaM2=1.9)),
                         ('wrong-candidate-SHA',lambda v:v['condition']['explicit'].update(candidateSHA256='0'*64))]:
        bad=copy.deepcopy(m);mutate(bad)
        reject(label,lambda bad=bad:audit_admission(bad,base/'real01'))
    for label,budget in [('memory-increase',dict(memoryMiB=16384)),('auto-continue',dict(automaticContinuation=1)),('time-increase',dict(nativeWallClockSeconds=1200))]:
        changed_candidate=copy.deepcopy(a);changed_candidate['budget'].update(budget)
        path=out/(label+'.json');path.write_text(json.dumps(changed_candidate))
        bad=copy.deepcopy(m);bad['condition']['explicit'].update(candidate=str(path),candidateSHA256=sha(path))
        reject(label,lambda bad=bad:audit_admission(bad,base/'real01'))
    bad=copy.deepcopy(original);bad['mesh']=m['mesh']
    reject('interface-profile-full-domain',lambda:audit_run_scope(bad,base/'real01'))
    # 新保存分派不能破坏已经交付的真实六步元数据读取；只读、不重算。
    result=audit_completion(ROOT/'outputs/v183_validation/small01',ROOT/'outputs/v183_validation/small01_run/process00/stdout.txt')
    assert result['savedSteps']==[0,3,6] and result['completedRequestedSteps']
    r=dict(passed=True,rejected=rejected,numericalChangeKeys=sorted(changed),legacyCompletionPassed=True,
           sourceSHA256={name:sha(ROOT/name) for name in ['scripts/eilmerExplicitAdmission.py','scripts/eilmerExplicitProfile.py','scripts/readEilmerNative.py','tests/testEilmerExplicitAdmission.py']},
           nativeRuns=0,physicalQualification=False,newQualifiedAreaM2=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);test(p.parse_args().out)
