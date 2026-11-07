"""入口测试：超域、错误场/警告必须拒绝，旧CEA参照默认不变。"""
from pathlib import Path
import sys,json,tempfile,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from openFoamAirThermo import evaluate,contract
from openFoamFields import internal_field,audit_state,SavedGuard
from hybridEulerReference import control_reference


def main():
    start=time.monotonic();checks=[]
    t=np.array([200.,300.,1000.,1500.,2500.]);g=evaluate(5000,t)
    state=dict(p=np.full(5,5000.),T=t,rho=g['densityKgm3'],e=g['internalEnergyJkg'],U=np.zeros((5,3)))
    audit_state(state);checks.append('valid-negative-reference-internal-energy')
    for key,values in [('p',[-1]*5),('rho',[0]*5),('T',[199.]*5),('T',[2501.]*5),
        ('e',g['internalEnergyJkg']+100),('rho',g['densityKgm3']*1.01),('U',np.full((5,3),np.nan))]:
        bad={k:np.array(v).copy() for k,v in state.items()};bad[key]=np.asarray(values)
        try:audit_state(bad)
        except ValueError:checks.append('rejected-'+key)
        else:raise AssertionError('invalid '+key+' accepted')
    with tempfile.TemporaryDirectory(dir=ROOT/'outputs/v198_validation') as tmp:
        p=Path(tmp)/'field';p.write_text('FoamFile { format ascii; }\ninternalField nonuniform List<vector> 2 ((1 2 3) (4 5 6));')
        assert np.array_equal(internal_field(p,2,True),[[1,2,3],[4,5,6]]);checks.append('vector-parser')
        try:internal_field(p,3,True)
        except ValueError:checks.append('rejected-count')
        else:raise AssertionError('wrong count accepted')
        log=Path(tmp)/'log';log.write_text('attempt to use janafThermo out of temperature range')
        assert SavedGuard(tmp,5,Path(tmp)/'guard.json')(log,False)['reason']=='native-warning-or-fatal';checks.append('warning-stop')
    source=subprocess.check_output(['git','show','06bfe219ab0a5d1e9da66f6d82ad7851a3256942:scripts/hybridEulerReference.py'],cwd=ROOT,text=True,encoding='utf-8')
    ns={};exec(compile(source,'historical-reference','exec'),ns)
    assert ns['control_reference']()==control_reference();checks.append('old-CEA-default-exact')
    from runOpenFoamControl import main as run_control
    try:run_control()
    except ValueError as error:
        assert '已关闭' in str(error);checks.append('closed-candidate-before-launch')
    else:raise AssertionError('failed candidate can still launch')
    result=dict(passed=True,checks=checks,contract=contract(),elapsedSeconds=time.monotonic()-start)
    out=ROOT/'outputs/v198_validation/entry_tests.json';out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,checks=len(checks),elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
