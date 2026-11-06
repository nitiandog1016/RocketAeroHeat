"""模块用途：只读绑定准确SU2运行源，核验自适应CFL语义；不启动或重放CFD。"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import time

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
BUILD=ROOT/'outputs/v171_validation/build01'
PACKAGE=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
PARAMETERS=(.1,1.2,.4,100.,.001,0.)
SOURCES={
    'Common/src/CConfig.cpp':[(1922,1944),(2065,2069),(3793,3806),(5945,5973)],
    'SU2_CFD/src/solvers/CSolver.cpp':[(113,120),(1894,2130)],
    'SU2_CFD/src/iteration/CFluidIteration.cpp':[(81,99),(412,458)],
    'SU2_CFD/src/solvers/CEulerSolver.cpp':[(1772,1824),(2567,2602)],
    'SU2_CFD/src/variables/CNSVariable.cpp':[(140,187)],
    'SU2_CFD/include/fluid/CSU2IdealAir.hpp':[(1,40)],
    'SU2_CFD/include/fluid/IdealAir.hpp':[(13,36),(120,165)],
    'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp':[(425,443),(552,575)],
    'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl':[(266,280),(587,605),(854,922),(993,1017)],
    'SU2_CFD/src/variables/CFlowVariable.cpp':[(80,89)],
    'SU2_CFD/src/integration/CMultiGridIntegration.cpp':[(747,771)],
    'SU2_CFD/src/output/CFlowCompOutput.cpp':[(435,446)],
    'Common/src/geometry/CGeometry.cpp':[],
    'config_template.cfg':[(1575,1592)],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def scalar_update(cfl,relaxation,linear_residual,iteration=0,reduce=False,reset=False,parameters=PARAMETERS):
    """模块一：源语义的有限标量算例，不称实际C++执行或真实局部CFL预测。"""
    down,up,lower,upper,acceptable,start=parameters
    tol=max(acceptable,1e-5)
    reduce=reduce or (linear_residual>1.2*tol and iteration>=start)
    reset=reset or linear_residual>.99
    increase=linear_residual<tol and iteration>=start
    factor=down if relaxation<.1 or reduce else (1. if relaxation<1 or not increase else up)
    result=max(lower,min(upper,cfl*factor))
    return lower if reset else result


def historical_feedback(history,start):
    """模块二：只重放既有历史的全局反馈谓词；真实松弛因子缺失且场不会随此重放变化。"""
    rows=history[start:];length=min(100,len(rows)-1);series=[0.]*length;counter=0;new=0.
    reduced=[];reset=[];linear_ok=0
    for iteration,row in enumerate(rows):
        old=new;new=sum(row[key] for key in ('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]'))/5
        if iteration==1 and new-old>10:old=new
        series[counter]=new-old;counter=(counter+1)%length
        lin=10**row['LinSolRes'];linear_ok+=lin<.001
        do_reduce=lin>.0012;do_reset=lin>.99
        if iteration>=length:
            total=sum(series);previous=series[0];changes=0
            for value in series:
                changes+=(previous>0)!=(value>0);previous=value
            do_reduce|=changes>length//4 and total>-.5
            if total>2:
                do_reset=True;counter=0;series=[0.]*length
        if do_reduce:reduced.append(int(row['Inner_Iter']))
        if do_reset:reset.append(int(row['Inner_Iter']))
    return dict(sourceFirstIteration=int(rows[0]['Inner_Iter']),sourceLastIteration=int(rows[-1]['Inner_Iter']),
        rows=len(rows),circularHistoryLength=length,linearAllowsIncreaseRows=int(linear_ok),
        globalReductionSourceIterations=reduced,globalResetSourceIterations=reset,
        localUnderRelaxationAvailable=False,actualAdaptiveFlowPredicted=False,
        caution='fixed-field history predicates only; missing local relaxation and different future solution prevent actual-CFL prediction')


def main(output):
    """模块三：来源SHA、准确选源、有限分支/解析增长核验与只读回执。"""
    started=time.monotonic();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    runtime_path=PACKAGE/'runtime.json';tree_path=BUILD/'source_tree.json'
    runtime=json.loads(runtime_path.read_text());tree=json.loads(tree_path.read_text())
    assert sha(tree_path)==runtime['sourceTreeSHA256']=='c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8'
    binary=PACKAGE/'bin/SU2_CFD.exe';assert sha(binary)==next(x['sha256'] for x in runtime['files'] if x['path']=='bin/SU2_CFD.exe')
    locks={str(p):sha(p) for p in (runtime_path,tree_path,binary,Path(__file__).resolve())}
    records=[];texts={};excerpts=[]
    for relative,intervals in SOURCES.items():
        original=SOURCE/relative;expected=tree[str(original)]
        candidates=[BUILD/'executed_sources'/relative,ROOT/'outputs/v163_validation/build03/executed_sources'/relative,original]
        selected=next((p for p in candidates if p.is_file() and sha(p)==expected),None)
        assert selected is not None,'缺准确源：'+relative
        locks[str(selected)]=expected;text=selected.read_text('utf8');texts[relative]=text
        records.append(dict(relativePath=relative,selectedSource=str(selected),sha256=expected,
            sharedSourceMatches=original.is_file() and sha(original)==expected))
        lines=text.splitlines()
        for first,last in intervals:
            excerpts.append('\n'+relative+':'+str(first)+'-'+str(last)+' SHA256='+expected+'\n')
            excerpts.extend(str(i+1)+': '+lines[i]+'\n' for i in range(first-1,min(last,len(lines))))
    (out/'exact_source_excerpts.txt').write_text(''.join(excerpts),encoding='utf8')
    cfg=texts['Common/src/CConfig.cpp'];defaults=[float(v) for _,v in re.findall(r'default_cfl_adapt\[(\d)\] = ([\d.]+);',cfg)]
    assert defaults==[.1,1.2,10.,100.,.001,0.]
    assert PARAMETERS[:2]==tuple(defaults[:2]) and PARAMETERS[3:]==tuple(defaults[3:])
    assert 'GetnMGLevels() != 0' in texts['Common/src/geometry/CGeometry.cpp']
    assert 'addDoubleOption("MAX_UPDATE_FLOW", MaxUpdateFlow, 0.2)' in cfg
    state=texts['SU2_CFD/src/variables/CNSVariable.cpp']
    assert state.index('FluidModel->SetTDState_rhoe(density, staticEnergy);')<state.index('if (check_dens || check_press')
    assert 'catch' not in state
    cases=[('growth',.4,1.,1e-7,False,False,.48),('partial-relaxation',2.,.5,1e-7,False,False,2.),
        ('strong-relaxation',20.,.05,1e-7,False,False,2.),('linear-reduction',20.,1.,.0013,False,False,2.),
        ('linear-deadband',20.,1.,.0011,False,False,20.),('linear-reset',20.,1.,1.,False,False,.4),
        ('nonlinear-reduction',20.,1.,1e-7,True,False,2.),('nonlinear-reset',20.,1.,1e-7,False,True,.4),
        ('upper-bound',99.,1.,1e-7,False,False,100.),('lower-bound',.4,.05,1e-7,False,False,.4)]
    checked=[]
    for name,cfl,relax,lin,reduce,reset,expected in cases:
        actual=scalar_update(cfl,relax,lin,reduce=reduce,reset=reset)
        assert abs(actual-expected)<1e-14
        checked.append(dict(name=name,actual=actual,expected=expected))
    first_default=scalar_update(.4,1.,1e-7,parameters=tuple(defaults));assert first_default==10.
    ramp=.4;analytic_error=0.;first_cap=None
    for count in range(1,41):
        ramp=scalar_update(ramp,1.,1e-7)
        analytic_error=max(analytic_error,abs(ramp-min(100.,.4*1.2**count)))
        if ramp==100. and first_cap is None:first_cap=count
    assert first_cap==math.ceil(math.log(100/.4)/math.log(1.2))==31 and analytic_error<1e-12
    history_path=ROOT/'outputs/v171_validation/real01/history.csv';locks[str(history_path)]=sha(history_path)
    with history_path.open() as f:raw=list(csv.DictReader(f,skipinitialspace=True))
    history=[{k.strip():float(v) for k,v in row.items()} for row in raw]
    assert [int(row['Inner_Iter']) for row in history]==list(range(1716))
    stable={path:sha(Path(path))==value for path,value in locks.items()};assert all(stable.values())
    report=dict(sourceReviewPassed=True,sourceCommit=runtime['sourceCommit'],repositoryBaseline='2865300feb83df5b093e8b14cd32b0f1d4b8cc5c',
        runtimeFlavor=runtime['flavor'],sources=records,nativeDefaults=defaults,proposedSingleParameters=list(PARAMETERS),
        soleDefaultParameterAdaptation='CFL lower bound 10 -> .4 to preserve existing start and avoid immediate floor jump',
        firstDefaultUpdateFromPoint4=first_default,scalarSemanticChecks=checked,rampAnalyticMaximumError=analytic_error,
        idealizedFullRelaxationUpdatesToUpperBound=first_cap,semanticChecksAreNotNativeExecution=True,
        fixedHistoryFeedback=historical_feedback(history,0),coldReferenceTailFeedback=historical_feedback(history,1501),
        nativeTemperatureIntervalPreservationProven=False,temperatureOutside60To2000KRemainsFatal=True,
        nativeCFLFeedbackCompatibleWithCurrentModel=True,coldRestartContinuousEquivalenceProven=False,
        sourceSHA256=sha(Path(__file__).resolve()),inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started,
        nativeCFDRuns=0,builds=0,newPhysicalQualification=False)
    save(out/'report.json',report)
    print(json.dumps(dict(passed=True,elapsedSeconds=report['elapsedSeconds'],sources=len(records),
        actualCFDRuns=0,coldReferenceTail=report['coldReferenceTailFeedback']),ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();main(args.output)
