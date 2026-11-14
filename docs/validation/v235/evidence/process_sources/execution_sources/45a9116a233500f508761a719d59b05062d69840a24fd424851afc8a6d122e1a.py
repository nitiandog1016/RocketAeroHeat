"""模块职责：绑定单连通生产循环准确源码，生成只读观测适配与固定构建输入。"""
from pathlib import Path
import json,re,shutil,subprocess,sys,difflib
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v235_validation';DOC=ROOT/'docs/validation/v235';GAS=ROOT/'outputs/v220_validation'
OLD=ROOT/'tools/runtime/su2-thermal-build/source';BUILDOLD=OLD.parent/'build-v117'
CLOSED=False
SPECS=['SU2_CFD/src/drivers/CSinglezoneDriver.cpp','SU2_CFD/src/iteration/CFluidIteration.cpp',
 'SU2_CFD/src/integration/CIntegration.cpp','SU2_CFD/src/integration/CMultiGridIntegration.cpp',
 'SU2_CFD/src/integration/CSingleGridIntegration.cpp','SU2_CFD/src/variables/CNSVariable.cpp',
 'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
 'SU2_CFD/src/solvers/CTurbSolver.cpp','SU2_CFD/src/solvers/CTurbSASolver.cpp',
 'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl']


def sources():
    target=BASE/'source01';target.mkdir(exist_ok=False);known=read(GAS/'build02/dependencies.json')
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');bound=[]
    for rel in SPECS:
        expected=known.get(str(GAS/'source'/rel),tree.get(str(OLD/rel)));assert expected,rel
        candidates=[GAS/'source'/rel,ROOT/'outputs/v222_validation/source01'/rel,
          ROOT/'outputs/v167_validation/build01/executed_sources'/rel,
          ROOT/'outputs/v219_validation/source01/source'/rel,OLD/rel]
        p=next((p for p in candidates if p.is_file() and sha(p)==expected),None)
        assert p,('EXACT_EXECUTED_SOURCE_MISSING',rel)
        dest=target/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        bound.append(dict(relative=rel,source=str(p),sha256=expected,snapshot=str(dest)))
    save(target/'bindings.json',bound)
    save(target/'report.json',dict(passed=True,sourceFiles=len(bound),flowUpdates=0,nativeExecution=False))


def main(stage):
    if CLOSED:raise SystemExit('v235 preparation closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      GAS/'build02/dependencies.json',ROOT/'outputs/v167_validation/build01/source_tree.json']
    if stage=='sources':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_sources'],120,128,locks,BASE)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v235 preparation closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_sources':sources()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
