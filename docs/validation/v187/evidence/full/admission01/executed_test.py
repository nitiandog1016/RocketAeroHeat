"""完整域新存储准入的来源/预算拒绝；只读复用小例，不启动原生计算。"""
from pathlib import Path
import argparse
import copy
import json
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from eilmerExplicitAdmission import audit_admission,sha
from eilmerExplicitProfile import settings
from eilmerWorkspaceRuntime import audit_workspace_execution as audit_dense_execution
from readEilmerNative import audit_run_scope


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    result=dict(passed=False,nativeRuns=0,physicalQualification=False)
    try:
        base=ROOT/'outputs/v187_validation/full';case=base/'real01'
        c=json.loads((base/'condition.json').read_text())
        m=json.loads((ROOT/'outputs/v178_validation/real01/case_manifest.json').read_text());m['condition']=c
        a,_=audit_admission(m,case)
        prior=ROOT/'outputs/v187_validation'
        old=json.loads((prior/'small01/case_manifest.json').read_text())
        oldrun=json.loads((prior/'small01_run/report.json').read_text())
        assert audit_dense_execution(old['condition'],oldrun)
        oldsettings=settings(old['condition']['explicit']);newsettings=settings(c['explicit'])
        delta={k for k in oldsettings if oldsettings[k]!=newsettings[k]}
        assert delta=={'max_step','write_loads_at_step','write_flow_solution_at_step'}
        rejected=[]
        def reject(label,fn):
            try:fn()
            except ValueError:rejected.append(label)
            else:raise AssertionError(label)
        audit_run_scope(m,case)
        reject('unregistered-directory',lambda:audit_admission(m,base/'other'))
        for label,mutate in [('changed-Mach',lambda v:v['condition'].update(Mach=6)),
            ('changed-mesh',lambda v:v['mesh'].update(sourceSHA256='0'*64)),
            ('changed-wall-count',lambda v:v['mesh'].update(wallFaces=19751)),
            ('changed-area',lambda v:v['mesh'].update(wallAreaM2=1.9)),
            ('wrong-candidate-SHA',lambda v:v['condition']['explicit'].update(candidateSHA256='0'*64)),
            ('wrong-runtime-SHA',lambda v:v['condition']['explicit'].update(runtimePackageSHA256='0'*64))]:
            bad=copy.deepcopy(m);mutate(bad)
            reject(label,lambda bad=bad:audit_admission(bad,case))
        for label,key,value in [('memory-increase','memoryMiB',16384),('auto-continue','automaticContinuation',1),('time-increase','nativeWallClockSeconds',1200)]:
            changed=copy.deepcopy(a);changed['budget'][key]=value
            path=out/(label+'.json');path.write_text(json.dumps(changed))
            bad=copy.deepcopy(m);bad['condition']['explicit'].update(candidate=str(path),candidateSHA256=sha(path))
            reject(label,lambda bad=bad:audit_admission(bad,case))
        changed=copy.deepcopy(a);changed['evidence'].pop('storageVerification')
        path=out/'missing-comparison.json';path.write_text(json.dumps(changed))
        bad=copy.deepcopy(m);bad['condition']['explicit'].update(candidate=str(path),candidateSHA256=sha(path))
        reject('missing-storage-comparison',lambda:audit_admission(bad,case))
        badrun=copy.deepcopy(oldrun);badrun['commands'][0]['inputSHA256'][str(ROOT/'scripts/eilmerWorkspaceRuntime.py')]='0'*64
        reject('forged-historical-source',lambda:audit_dense_execution(old['condition'],badrun))
        v186=ROOT/'outputs/v186_validation/real01'
        reject('v186-still-closed',lambda:audit_run_scope(json.loads((v186/'case_manifest.json').read_text()),v186))
        legacy=json.loads((ROOT/'outputs/v184_validation/real01/case_manifest.json').read_text())
        reject('v184-still-closed',lambda:audit_run_scope(legacy,ROOT/'outputs/v184_validation/real01'))
        small=copy.deepcopy(old);small['mesh']=m['mesh']
        reject('small-profile-full-domain',lambda:audit_run_scope(small,case))
        result.update(passed=True,rejected=rejected,numericalChangeKeys=sorted(delta),oldExactSourceAccepted=True,
            sourceSHA256={n:sha(ROOT/n) for n in ['scripts/eilmerExplicitAdmission.py','scripts/eilmerExplicitProfile.py',
                'scripts/eilmerWorkspaceRuntime.py','scripts/readEilmerNative.py','tests/testEilmerDenseAdmission.py']})
    except Exception as error:
        result['error']=repr(error);raise
    finally:
        result['elapsedSeconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);test(p.parse_args().out)
