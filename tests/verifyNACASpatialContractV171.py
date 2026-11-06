"""模块用途：新空间Tw组合的启动前拒绝、预算与旧分派隔离审核，不运行CFD。"""
import argparse
import copy
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2NACASpatialInterface import validate,read_profile,resource_budget
from su2ThermalRuntime import read,sha,resolve
from prepareNACASpatialCase import FLAVOR,SPATIAL,KIND,specification
from su2IdealAirInterface import settings
from validateIdealAir import save


def verify(folders,out):
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    rejected=[];checked=[];locks={str(Path(__file__).resolve()):sha(Path(__file__).resolve())}
    for folder in folders:
        m=read(folder/'case_manifest.json');text=(folder/'case.cfg').read_text()
        exe,runtime,flavor,profile,inputs=resolve(folder,text);locks.update(inputs)
        request=read(folder/'request.json');resource_budget(request,m)
        variant=m['specification']['variant'];checked.append(dict(variant=variant,flavor=flavor,hasProfile=profile is not None))
        for name,wrong in [('axis',text.replace('AXISYMMETRIC= NO','AXISYMMETRIC= YES')),
                           ('SST',text.replace('SOLVER= NAVIER_STOKES','SOLVER= RANS')),
                           ('residual',text.replace('CONV_RESIDUAL_MINVAL= -8','CONV_RESIDUAL_MINVAL= -3')),
                           ('unplanned-updates',text.replace('ITER= 4\n','ITER= 5\n')),
                           ('wrong-Mach-vector',text.replace(', -1.97, 0, 0, outlet',', -1, 0, 0, outlet'))]:
            assert wrong!=text,name
            try:validate(folder,m,wrong,runtime)
            except AssertionError:rejected.append(variant+':'+name)
            else:raise AssertionError('未拒绝 '+name)
        for name,change in [('timeout',dict(maximumSeconds=31)),('threads',dict(threads=3)),
                            ('wrong-profile',dict(resourceProfile='v171-naca-spatial-validation'))]:
            bad=dict(request,**change)
            try:resource_budget(bad,m)
            except AssertionError:rejected.append(variant+':'+name)
            else:raise AssertionError('未拒绝 '+name)
        bad=copy.deepcopy(m);bad['specification']['updates']=5
        try:validate(folder,bad,text,runtime)
        except AssertionError:rejected.append(variant+':altered-specification')
        else:raise AssertionError('未拒绝步数合同')
        bad_runtime=dict(runtime,flavor='rah-ideal-air-v1')
        try:validate(folder,m,text,bad_runtime)
        except AssertionError:rejected.append(variant+':wrong-runtime')
        else:raise AssertionError('未拒绝运行包')
    # 专用表解析拒绝：不碰原任务输入，所有错误表独占本测试目录。
    examples=dict(header='WRONG\n-1,335\n0,335\n',duplicate=SPATIAL['schema']+'\n0,335\n0,336\n',
        nonfinite=SPATIAL['schema']+'\n-1,nan\n0,335\n',low=SPATIAL['schema']+'\n-1,59\n0,335\n',
        high=SPATIAL['schema']+'\n-1,335\n0,2001\n',coverage=SPATIAL['schema']+'\n-1,335\n-0.1,335\n')
    for name,content in examples.items():
        file=out/(name+'.csv');file.write_text(content,encoding='ascii')
        try:read_profile(file,np.array([-.127,0]))
        except (AssertionError,ValueError):rejected.append('profile:'+name)
        else:raise AssertionError('未拒绝错误表 '+name)
    # 实际网格预算尚未启动；独立容量公式与允许的唯一登记输入核对。
    from prepareSU2NACALaminar import capacity
    plan=capacity(190334,3226,3096,2000,500,1024*1024**2)
    assert plan['maximumNamedPairs']==5 and plan['estimatedMaximumBytes']==1015107072
    real=dict(caseKind=KIND,specification=specification('validation-profile'),requiredRuntimeFlavor=FLAVOR,nodes=190334,
        planSHA256=sha(ROOT/'docs/validation/v171/PLAN.md'),idealAir=dict(checkpointContract=plan))
    good=dict(maximumSeconds=1800,threads=2,resourceProfile='v171-naca-spatial-validation');resource_budget(good,real)
    for name,change in [('longer',dict(maximumSeconds=1801)),('more-threads',dict(threads=4))]:
        try:resource_budget(dict(good,**change),real)
        except AssertionError:rejected.append('actual-budget:'+name)
        else:raise AssertionError('未拒绝实际预算')
    assert all(sha(Path(p))==h for p,h in locks.items())
    report=dict(passed=True,variants=checked,rejected=rejected,actualEstimatedOutputBytes=plan['estimatedMaximumBytes'],
        inputSHA256=locks,inputStable=True,nativeTasks=0,physicalQualified=False)
    save(out/'report.json',report);print(dict(passed=True,variants=len(checked),rejected=len(rejected)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folders',nargs='+',required=True)
    parser.add_argument('--output',required=True);args=parser.parse_args()
    verify([Path(p).resolve() for p in args.folders],Path(args.output).resolve())
