"""模块职责：绑定成熟源/对象与固定PETSc版本，输出一次真实NK整体控制的采用依据。"""
from pathlib import Path
import hashlib,sys
import prepareOuterStencilV249 as task
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=task.BASE;DOC=task.DOC;CLOSED=True
SRC=ROOT/'tools/runtime/su2-thermal-build/source'

SPECS=[(SRC/'SU2_CFD/src/integration/CNewtonIntegration.cpp',[(66,117),(120,161),(186,229),(267,330),(365,394)]),
 (ROOT/'outputs/v220_validation/source/Common/include/CConfig.hpp',[(439,444)]),
 (ROOT/'outputs/v220_validation/source/Common/src/CConfig.cpp',[(1904,1911),(3807,3820)]),
 (SRC/'SU2_CFD/src/solvers/CSolverFactory.cpp',[(272,286)]),
 (SRC/'SU2_CFD/src/integration/CIntegrationFactory.cpp',[(50,69)]),
 (SRC/'SU2_CFD/include/variables/CVariable.hpp',[(335,346)]),
 (ROOT/'outputs/v235_validation/source01/SU2_CFD/src/iteration/CFluidIteration.cpp',[(60,101),(118,139)])]


def review():
    out=BASE/'method01';out.mkdir(exist_ok=False);manifest=read(BASE/'external_sources/manifest.json');external=[]
    assert manifest['ref']=='v3.26.1'
    for entry in manifest['sources']:
        p=BASE/'external_sources'/entry['path'].replace('/','__');raw=p.read_bytes()
        blob=hashlib.sha1(('blob %d\0'%len(raw)).encode()+raw).hexdigest();assert blob==entry['sha'],entry['path']
        external.append(dict(path=entry['path'],localPath=str(p),url=entry['display_url'],gitBlob=blob,sha256=sha(p)))
    excerpts=[]
    for n,(p,ranges) in enumerate(SPECS):
        lines=p.read_text(encoding='utf8').splitlines();h=sha(p)
        for first,last in ranges:
            name='%02d_%s_%d.txt'%(n,p.stem,first);content='\n'.join('%6d %s'%(i+1,lines[i]) for i in range(first-1,last))+'\n'
            (out/name).write_text(content,encoding='utf8');excerpts.append(dict(path=str(p),sha256=h,firstLine=first,lastLine=last,excerpt=name))
    members=read(ROOT/'outputs/v245_validation/build01/members.json');objects=[]
    for rows in members.values():
        for item in rows:
            if any(word in item['path'] for word in ('CNewtonIntegration','CSolverFactory','CIntegrationFactory')):
                assert sha(item['path'])==item['sha256'];objects.append(item)
    assert len(objects)==3
    n=(SPECS[0][0]).read_text(encoding='utf8')
    assert 'TODO: Clever back-tracking and CFL adaptation based on residual reduction.' in n
    assert 'CompleteImplicitIteration' in n and 'FGMRES_LinSolver' in n and 'ResEvalType::EXPLICIT' in n
    version=(BASE/'external_sources/include__petscversion.h').read_text();assert 'PETSC_VERSION_MINOR      26' in version
    report=dict(passed=True,adoptedForNextBoundedControl='existing SU2 8.6 native Newton-Krylov + unchanged staggered SA',
      productionAdopted=False,proofOfImprovedConvergence=False,fullResidualMeritLineSearch=False,
      PETSc=dict(version='3.26.1',refObject=manifest['tagObject'],sources=external,
        capabilities=['separate configurable level smoothers','optional coarse correction line search','function-domain error handling'],
        currentDecision='defer full replacement: no existing six-state SU2/PETSc per-level function and preconditioner bridge',
        variableBoundsSupportedByFAS=False,newInstallation=False,allSourcesBlobVerified=True),
      nativeSourceExcerpts=excerpts,nativeObjectsReused=objects,
      nativeSemantics=dict(flowVariables=5,SA='separate native auxiliary solve',pseudoTimeAddedToJv=True,
        perturbation='Solution_Old + scaled direction, not accumulated perturbations',
        finiteDifferenceStep='native norm of internal conservative coordinates; not claimed zero-offset invariant',
        defaultStartupWouldHideNK=True,controlledStartupIterations=0,controlledStartupResidual=0,
        domainRecoveryDuringProductsMustStopControl=True,nativeBacktrackingImplemented=False),
      control=dict(input='v248 accepted saved state200 as research-state method control',
        notContinuationOfOldFAS=True,initialSixStateMustMatch=True,resetCFL=.1,updates=1,
        newCFDNotYetAdmitted=True,newThermalQualification=False),physicalQualification=False)
    save(out/'report.json',report);print(dict(passed=True,adopted=report['adoptedForNextBoundedControl'],PETScVersion='3.26.1',nativeBacktracking=False))


def main():
    if CLOSED:raise SystemExit('v249 source review closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(task.__file__),DOC/'PLAN.md',ROOT/'outputs/v245_validation/build01/members.json']
    locks += [p for p,_ in SPECS]+[p for p in (BASE/'external_sources').iterdir() if p.is_file()]
    return task.setup_guard().supervise('method',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v249 source review closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':review()
    else:raise SystemExit(0 if main() else 1)
