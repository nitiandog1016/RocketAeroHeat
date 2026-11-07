"""显式热壁接口的独立原面复核及故障拒绝；不启动求解器。"""
from pathlib import Path
import argparse,copy,csv,hashlib,json,shutil,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import audit_prepared,export_loads,read_field,audit_run_scope
from auditEilmerCompletion import audit_completion
from eilmerExplicitProfile import audit_explicit_config
from prepareEilmerRocketCase import validate_condition,make_job


def verify(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'outputs/v183_validation';case=base/'small01';receipt=base/'small01_recovery/report.json'
    log=base/'small01_run/process00/stdout.txt'
    before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in case.rglob('*') if p.is_file()}
    m=json.loads((case/'case_manifest.json').read_text());c=m['condition']
    geometry,_,_=audit_prepared(case);completion=audit_completion(case,log)
    r=json.loads(receipt.read_text());assert r['passed'] and r['newNativeRuns']==0
    assert r['commands'][0]['argv'][0]=='lmr-run'
    assert not r['commands'][0]['liveGroupAfterCleanup']
    assert json.loads((base/'small01_run/report.json').read_text())['passed'] is False
    # 模块一：从原生有符号热分项独立复核原面CSV，不用导出结果自身当真值。
    load=case/'lmrsim/loads/0001/blk-0000-bndry-0-wall.dat'
    wall=export_loads(case,load,out/'wall.csv',receipt)
    a=np.atleast_1d(np.genfromtxt(load,names=True,deletechars=''))
    rows=list(csv.DictReader((out/'wall.csv').open()))
    err=0.0
    for row in rows:
        j=int(row['nativeBoundaryRow'])
        expected=-a['outsign'][j]*(a['q_cond'][j]+a['q_diff'][j])
        err=max(err,abs(float(row['heatFluxIntoWall_W_m2'])-expected)/max(1,abs(expected)))
        assert float(row['givenWallTemperature_K'])==300 and row['status']=='research'
        assert row['physicsQualified']=='False' and 'local pseudo-time' in row['method']
    assert err<2e-12 and wall['wallFaces']==7 and wall['newQualifiedAreaM2']==0
    fields=[read_field(case,i,4) for i in [0,1,2]]
    pressure_change=float(np.max(np.abs(fields[-1]['p']-fields[0]['p'])))
    assert pressure_change>0
    assert all(np.all(f['tq-tke']>=0) and np.all(f['tq-omega']>0) for f in fields)
    # 模块二：原始二进制数据缺失/非有限/零步长与伪成功日志均拒绝。
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (ValueError,FileNotFoundError):rejected.append(label)
        else:raise AssertionError('未拒绝: '+label)
    text=log.read_text()
    for label,changed in [
        ('native-integration-failed',text+'\nNote that integrateInTime failed.\n'),
        ('caught-exception',text+'\nAn exception was caught during the time step.\n'),
        ('wrong-number-type',text.replace('Number-type: real','Number-type: complex')),
        ('missing-final-marker',text.replace('FINAL-STEP: 6','NO-FINAL-STEP: 6')),
        ('incomplete-step-window',text.replace('FINAL-STEP: 6','FINAL-STEP: 5')),
        ('missing-time-step-history',text.replace('Step=      2 t=','Other=      2 t=')),
        ('wrong-pseudo-time',text.replace('FINAL-TIME: 1.41987e-05','FINAL-TIME: 0.5'))]:
        path=out/(label+'.txt');path.write_text(changed)
        reject(label,lambda path=path:audit_completion(case,path))
    cfg=json.loads((case/'lmrsim/config').read_text());ctrl=json.loads((case/'lmrsim/control').read_text())
    for key,val in [('with_local_time_stepping',False),('gasdynamic_update_scheme','euler'),('solver_mode','steady')]:
        bad=copy.deepcopy(cfg);bad[key]=val
        reject(key,lambda bad=bad:audit_explicit_config(bad,c,ctrl))
    for key,val in [('viscous_signal_factor',0),('turbulent_signal_factor',0),('stringent_cfl',False)]:
        bad=copy.deepcopy(ctrl);bad[key]=val
        reject(key,lambda bad=bad:audit_explicit_config(cfg,c,bad))
    for label in ['missing-field','zero-local-dt','truncated-load','wrong-load-step']:
        dest=out/label;shutil.copytree(case,dest)
        if label=='missing-field':(dest/'lmrsim/snapshots/0002/fluid-0000').rename(dest/'saved-field')
        elif label=='zero-local-dt':
            data=np.fromfile(dest/'lmrsim/snapshots/0002/fluid-0000',dtype='<f8').reshape(-1,4)
            data[list(fields[-1]).index('dt_local'),0]=0
            data.tofile(dest/'lmrsim/snapshots/0002/fluid-0000')
        elif label=='truncated-load':
            path=dest/'lmrsim/loads/0001/blk-0000-bndry-0-wall.dat'
            path.write_text('\n'.join(path.read_text().splitlines()[:-1])+'\n')
        else:
            path=dest/'lmrsim/loads/loads-metadata'
            path.write_text(path.read_text().replace(' 6\n',' 5\n'))
        reject(label,lambda dest=dest:audit_completion(dest,log))
    bad=copy.deepcopy(m);bad['mesh']['cells']=705490
    reject('unregistered-full-domain',lambda:audit_run_scope(bad))
    bad=copy.deepcopy(c);bad['explicit']['maxSteps']=100
    reject('unregistered-long-run',lambda:validate_condition(bad))
    bad=copy.deepcopy(c);bad['maxNewtonSteps']=6
    reject('mixed-NK-explicit-input',lambda:validate_condition(bad))
    # 模块三：共享生成/完成分派的旧方案回归仅只读复用，不重跑NK。
    legacy=ROOT/'outputs/v181_validation/small01'
    lm=json.loads((legacy/'case_manifest.json').read_text())
    assert make_job(lm['mesh']['nativeSourcePath'],lm['condition'])==(legacy/'job.lua').read_text()
    old_completion=audit_completion(legacy,ROOT/'outputs/v181_validation/small01_run/process00/stdout.txt')
    assert old_completion['finalStep']==8
    after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in case.rglob('*') if p.is_file()}
    assert before==after
    h=np.loadtxt(case/'lmrsim/transient-residuals.dat')
    result=dict(passed=True,originalFailedReceiptPreserved=True,nativeCompletion=completion,wall=wall,
                independentSignedHeatCSVScaledError=err,pressureChangePa=pressure_change,
                lastOverFirstMassResidual=float(h[-1,3]/h[0,3]),lastOverFirstEnergyResidual=float(h[-1,7]/h[0,7]),
                rejected=rejected,sourceCaseStable=True,legacyNKTextAndCompletionPreserved=True,
                newNativeRuns=0,physicalQualification=False,newQualifiedAreaM2=0,
                note='Six-step interface only; residual increase retained, no convergence or full-domain/heat-gradient qualification.',
                elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result));return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
