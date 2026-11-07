"""只读确认新入口真实切换精度并保存完整帧；不以此授予小域物理资格。"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import read_field,audit_numerics


def verify(case,receipt,out):
    start=time.monotonic();case=Path(case);out=Path(out);assert not out.exists()
    run=json.loads(Path(receipt).read_text());assert run['passed'] and not run['commands'][0]['liveGroupAfterCleanup']
    m=json.loads((case/'case_manifest.json').read_text());nk=json.loads((case/'lmrsim/nk-config').read_text());audit_numerics(nk,m['condition'])
    h=np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
    assert np.array_equal(h['step'],np.arange(1,9)) and np.array_equal(h['phase'],[0,1,1,1,1,1,1,1])
    assert h['CFL'][0]==.1 and h['CFL'][1]==1
    loads=np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1)
    assert np.array_equal(loads[:,2],[2,4,6,8])
    snapshots=np.atleast_2d(np.loadtxt(case/'lmrsim/snapshots/restart-metadata'))
    assert np.array_equal(snapshots[:,0],[4,8])
    for slot in range(3):read_field(case,slot,4)
    rejected=[]
    for key,val in [('number_of_phases',1),('max_cfl',.1),('steps_between_snapshots',99)]:
        bad=dict(nk);bad[key]=val
        try:audit_numerics(bad,m['condition'])
        except ValueError:rejected.append(key)
        else:raise AssertionError(key)
    r=dict(passed=True,nativeSteps=8,diagnosticPhase1StartsAt=2,loadsSteps=[2,4,6,8],restartSteps=[4,8],
        phaseMeaning='phase as recorded by native diagnostics; not an assumption about within-step update ordering',
        rejected=rejected,lastGlobalRelativeResidual=float(h['global-residual-rel'][-1]),physicsQualified=False,
        elapsedSeconds=time.monotonic()-start)
    out.write_text(json.dumps(r,indent=2)+'\n');return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--receipt',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(verify(a.case,a.receipt,a.out)))
