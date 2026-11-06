"""模块用途：只读比较新旧同源小场，复核阶段输出登记和真实文件；不重跑CFD。"""
import argparse
import csv
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
    old_frames=orr['checkpoints']['complete']
    for a,b in zip(cp['complete'],old_frames):
        assert {k:v for k,v in a.items() if not k.endswith('SHA256')}=={k:v for k,v in b.items() if not k.endswith('SHA256')}
    assert (new/'case.cfg').read_bytes()==(old/'case.cfg').read_bytes()
    measured=output_length(new,nr['inputSHA256'])
    limit=read(new/'case_manifest.json')['idealAir']['checkpointContract']['maximumOutputBytes']
    assert measured<=limit==32*1024**2
    rows=[];maximum=0.;all_equal=True
    for frame in cp['complete']:
        for name in (frame['volume'],frame['surface']):
            a=new/name;b=old/name
            for p,r in ((a,nr),(b,orr)):
                locks[str(p)]=sha(p);assert locks[str(p)]==r['outputSHA256'][name]
            equal=a.read_bytes()==b.read_bytes();all_equal=all_equal and equal
            differences={}
            if name.endswith('.csv'):
                av=np.loadtxt(str(a),delimiter=',',skiprows=1);bv=np.loadtxt(str(b),delimiter=',',skiprows=1)
                with a.open(encoding='ascii') as stream:columns=next(csv.reader(stream))
                with b.open(encoding='ascii') as stream:assert columns==next(csv.reader(stream))
                assert av.shape==bv.shape and np.array_equal(av[:,:4],bv[:,:4])
                assert np.isfinite(av).all() and np.isfinite(bv).all()
                fields=av.shape[1];points=len(av)
                for j,key in enumerate(columns):
                    differences[key]=dict(maximumAbsoluteDifference=float(np.max(abs(av[:,j]-bv[:,j]))),
                        maximumScaledDifference=float(np.max(abs(av[:,j]-bv[:,j])/np.maximum(1,abs(bv[:,j])))))
            else:
                ax,aa=vtk_tokens(a);bx,bb=vtk_tokens(b)
                assert np.array_equal(ax,bx) and set(aa)==set(bb)
                for key in aa:
                    assert aa[key].shape==bb[key].shape and np.isfinite(aa[key]).all() and np.isfinite(bb[key]).all()
                    differences[key]=dict(maximumAbsoluteDifference=float(np.max(abs(aa[key]-bb[key]))),
                        maximumScaledDifference=float(np.max(abs(aa[key]-bb[key])/np.maximum(1,abs(bb[key])))))
                fields=len(aa);points=len(ax)
            maximum=max(maximum,max(v['maximumScaledDifference'] for v in differences.values()))
            rows.append(dict(file=name,byteIdentical=equal,fields=fields,points=points,differences=differences))
    assert not (new/'conservative_steps.csv').exists(),'本小例若触发限步须独立分析，不能称合法路径对照'
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    save(out/'report.json',dict(passed=True,comparisonCompleted=True,comparisons=rows,allSavedFieldsExactlyEqual=all_equal,
        exactEqualityGatePassed=all_equal,maximumScaledDifference=maximum,roundoffAsSoleCauseProven=False,
        comparisonMeaning='measured differences only; exact-equality failure retained; no new tolerance gate',
        inputSHA256=locks,inputsStable=stable,offlineCheckpoints=cp,
        newFlavorInNativePhaseRegistry=True,nativeCheckpointsRecorded='checkpoints' in nr,
        nativeOutputBudgetRecorded='generatedOutputBytes' in nr,
        offlineGeneratedOutputBytes=measured,maximumOutputBytes=limit,
        originalSmallNativeOutputCapEnforcementClaimed=False,actualStepEventsObserved=False,
        newCFDTasks=0,physicalQualified=False,elapsedSeconds=time.monotonic()-started))
    print(dict(comparisonCompleted=True,identicalFiles=sum(r['byteIdentical'] for r in rows),
        files=len(rows),maximumScaledDifference=maximum,offlineOutputBytes=measured))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    compare(parser.parse_args().output)
