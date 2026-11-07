"""离散残差读取拒绝检查：损坏边界、场身份及中断步不能冒充有效数据。"""
from pathlib import Path
import json,sys,time,tempfile,subprocess,shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from openFoamConservation import array_entry,boundary_values
from openFoamFields import audit_state,internal_field
from openFoamAirThermo import evaluate
from auditOpenFoamDevelopment import log_progress


def main():
    start=time.monotonic();out=ROOT/'outputs/v199_validation/reader_tests02';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py');rows=[]
    for name in ['openFoamConservation.py','auditOpenFoamDevelopment.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    def reject(name,action):
        try:action()
        except (ValueError,AssertionError) as error:rows.append(dict(name=name,rejected=True,error=str(error)));return
        raise AssertionError('应拒绝但通过：'+name)

    # 模块一：实际value的完整性，不接受freestreamValue代替写出边界。
    assert np.array_equal(array_entry('value uniform (1 2 3);','value',2,True),[[1,2,3],[1,2,3]])
    assert array_entry('freestreamValue uniform 300;','value',1) is None
    assert np.array_equal(array_entry('value nonuniform List<scalar> 2 (1 2);','value',2),[1,2])
    for name,text,count,vector in [
        ('uniform_dimension','value uniform (1 2);',1,True),
        ('count','value nonuniform List<scalar> 2 (1 2);',3,False),
        ('type','value nonuniform List<vector> 1 ((1 2 3));',1,False),
        ('truncated','value nonuniform List<scalar> 2 (1);',2,False),
        ('nonfinite','value uniform nan;',1,False),
        ('trailing_corruption','value uniform 300 broken;',1,False)]:
        reject(name,lambda t=text,c=count,v=vector:array_entry(t,'value',c,v))
    with tempfile.TemporaryDirectory(dir=out) as temp:
        d=Path(temp);g=dict(neighbor=np.array([],int),owner=np.array([0]),
            patches=[dict(name='wall',startFace=0,nFaces=1)],sf=np.array([[1.,0,0]]),area=np.ones(1))
        def boundary(name,body,internal):
            (d/name).write_text('boundaryField { wall { '+body+' } }')
            return boundary_values(d,name,np.array(internal),g)
        assert np.array_equal(boundary('U','type slip;',[[1.,2,3]]),[[0,2,3]])
        assert np.array_equal(boundary('T','type zeroGradient;',[300]),[300])
        reject('missing_actual',lambda:boundary('T','type freestream; freestreamValue uniform 300;',[300]))
        reject('missing_type',lambda:boundary('T','value uniform 300;',[300]))
        reject('surface_owner_substitution',lambda:boundary('phi','type zeroGradient;',[0]))
        (d/'T').write_text('boundaryField { other { type zeroGradient; } }')
        reject('missing_patch',lambda:boundary_values(d,'T',np.array([300]),g))
        g['patches']=[]
        reject('uncovered_boundary',lambda:boundary_values(d,'T',np.array([300]),g))
        (d/'T').write_text('format ascii; internalField nonuniform List<scalar> 2 (300);')
        reject('internal_count',lambda:internal_field(d/'T',2))
    gas=evaluate(5000,300)
    fields={k:np.array([v]) for k,v in [('p',5000.),('T',300.),('rho',float(gas['densityKgm3'])),('e',float(gas['internalEnergyJkg']))]}
    fields['U']=np.array([[1.,2,3]]);audit_state(fields)
    changed={k:v.copy() for k,v in fields.items()};changed['e']+=100
    reject('altered_internal_energy',lambda:audit_state(changed))
    changed={k:v.copy() for k,v in fields.items()};changed['rho']*=1.1
    reject('altered_density',lambda:audit_state(changed))

    # 模块二：最后只有Time行时不计完成；关闭候选必须在调用原生进程前拒绝。
    progress=log_progress('Time = 0.1s\nExecutionTime = 1 s\nTime = 0.2s\n')
    assert progress['announcedSteps']==2 and progress['completedLogSteps']==1
    assert progress['lastCompletedLogTimeSeconds']==.1
    reject('orphan_completion',lambda:log_progress('ExecutionTime = 1 s\n'))
    reject('unfinished_before_next_step',lambda:log_progress('Time = 0.1s\nTime = 0.2s\n'))
    reject('decreasing_time',lambda:log_progress('Time = 0.2s\nExecutionTime = 1 s\nTime = 0.1s\n'))
    gaps=log_progress('Time = 0.1s\nExecutionTime = 1 s\nExecutionTime = 2 s\nTime = 0.2s\n',inspect_gaps=True)
    assert gaps['completedLogSteps']==1 and gaps['completionMarkerCount']==2 and gaps['unpairedCompletionLines']==[3]
    assert not gaps['logTimePairingComplete']
    for flags in [[],['--development'],['--development','--recover-observation']]:
        p=subprocess.run([sys.executable,str(ROOT/'scripts/openFoamSynchronousControl.py'),*flags],capture_output=True,timeout=10)
        # 字节匹配避免Windows终端编码影响；拒绝发生在导入WSL监督器之前。
        assert p.returncode!=0 and b'ValueError' in p.stderr and b'ModuleNotFoundError' not in p.stderr
        rows.append(dict(name='closed_cli_'+'_'.join(flags),rejected=True))
    result=dict(passed=True,rejections=rows,rejectionCount=len(rows),positiveChecks=8,
        CFDRuns=0,engineBuilds=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['passed','rejectionCount','positiveChecks','elapsedSeconds']}))


if __name__=='__main__':main()
