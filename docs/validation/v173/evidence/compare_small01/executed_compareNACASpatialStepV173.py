"""模块用途：只读比较新旧同源小场，复核阶段输出登记和真实文件；不重跑CFD。"""
import argparse
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import read,sha
from su2IdealAirInterface import PHASE_FLAVORS
from su2IdealAirCheckpoints import inventory,output_length
from prepareNACASpatialStep import FLAVOR
from auditSU2SSTState import vtk_tokens
from validateIdealAir import save


def compare(output):
    """模块一：新分派真实清点；模块二：全部体场/壁场与旧已运行对照逐值核对。"""
    started=time.monotonic();out=Path(output).resolve();out.relative_to(ROOT/'outputs')
    out.mkdir(parents=True,exist_ok=False)
    new=ROOT/'outputs/v173_validation/small01';old=ROOT/'outputs/v172_validation/small02'
    paths=[Path(__file__).resolve(),ROOT/'scripts/su2IdealAirInterface.py',ROOT/'scripts/su2IdealAirCheckpoints.py',
        ROOT/'scripts/su2NACASpatialStepInterface.py',new/'engine_process.json',old/'engine_process.json',
        new/'case_manifest.json',old/'case_manifest.json',new/'case.cfg',old/'case.cfg']
    locks={str(p):sha(p) for p in paths};save(out/'input_manifest.json',locks)
    for p in paths:
        if p.suffix=='.py':(out/('executed_'+p.name)).write_bytes(p.read_bytes())
    nr=read(new/'engine_process.json');orr=read(old/'engine_process.json')
    assert nr['status']==orr['status']=='terminal' and nr['exitCode']==orr['exitCode']==0
    assert all(nr['inputsStable'].values()) and all(orr['inputsStable'].values())
    assert FLAVOR in PHASE_FLAVORS
    cp=inventory(new,nr)
    assert cp['normalScheduleComplete'] and cp['normalFinalPairPresent'] and not cp['rejected']
    assert [r['localIteration'] for r in cp['complete']]==[0,2,4,6,7]
    assert cp['complete']==orr['checkpoints']['complete'], '同源合法小场的全部双帧SHA应相同'
    assert (new/'case.cfg').read_bytes()==(old/'case.cfg').read_bytes()
    measured=output_length(new,nr['inputSHA256'])
    limit=read(new/'case_manifest.json')['idealAir']['checkpointContract']['maximumOutputBytes']
    assert measured<=limit==32*1024**2
    rows=[];maximum=0.
    for frame in cp['complete']:
        for name in (frame['volume'],frame['surface']):
            a=new/name;b=old/name
            for p,r in ((a,nr),(b,orr)):
                locks[str(p)]=sha(p);assert locks[str(p)]==r['outputSHA256'][name]
            assert a.read_bytes()==b.read_bytes(), '合法路径文件出现变化: '+name
            if name.endswith('.csv'):
                av=np.loadtxt(str(a),delimiter=',',skiprows=1);bv=np.loadtxt(str(b),delimiter=',',skiprows=1)
                assert np.array_equal(av,bv);fields=av.shape[1];points=len(av)
            else:
                ax,aa=vtk_tokens(a);bx,bb=vtk_tokens(b)
                assert np.array_equal(ax,bx) and set(aa)==set(bb)
                assert all(np.array_equal(aa[k],bb[k]) for k in aa)
                fields=len(aa);points=len(ax)
            rows.append(dict(file=name,byteIdentical=True,fields=fields,points=points))
    assert not (new/'conservative_steps.csv').exists(),'本小例若触发限步须独立分析，不能称合法路径对照'
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    save(out/'report.json',dict(passed=True,comparisons=rows,allSavedFieldsExactlyEqual=True,
        inputSHA256=locks,inputsStable=stable,offlineCheckpoints=cp,
        newFlavorInNativePhaseRegistry=True,nativeCheckpointsRecorded='checkpoints' in nr,
        nativeOutputBudgetRecorded='generatedOutputBytes' in nr,
        offlineGeneratedOutputBytes=measured,maximumOutputBytes=limit,
        originalSmallNativeOutputCapEnforcementClaimed=False,actualStepEventsObserved=False,
        newCFDTasks=0,physicalQualified=False,elapsedSeconds=time.monotonic()-started))
    print(dict(passed=True,identicalFiles=len(rows),offlineOutputBytes=measured))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    compare(parser.parse_args().output)
