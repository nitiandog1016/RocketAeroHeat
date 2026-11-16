"""模块职责：绑定真实拒绝与准确成熟源码，比较原生恢复和带线搜索的粗层修正。"""
from pathlib import Path
import hashlib,sys
import runFASDevelopmentV247 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;OUT=BASE/'method_review';CLOSED=True


def review():
    event=read(BASE/'audit03/actual_rejection.json')['events'][0]
    assert event['nativeRecoverySucceeded'] and event['returnedFiveStateBitEqualOld']
    specs=[('outputs/v235_validation/stage01/CNSVariable.cpp','outputs/v235_validation/stage01/dependencies.json',
        'if (rah_single_observer::rejection) rah_single_observer::rejection(this,iPoint,trial.data(),restored.data(),true,rejected_reason.c_str());'),
      ('outputs/v245_validation/stage01/CMultiGridIntegration_observed.cpp','outputs/v245_validation/stage01/dependencies.json',
        'SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);')]
    source=[]
    for name,manifest,needle in specs:
        p=ROOT/name;expected=read(ROOT/manifest)[str(p)];assert sha(p)==expected
        lines=p.read_text(encoding='utf8').splitlines();indices=[i for i,s in enumerate(lines) if needle in s];assert len(indices)==1
        i=indices[0];source.append(dict(path=name,sha256=expected,startLine=i-3,endLine=i+7,excerpt='\n'.join(lines[max(0,i-4):i+7])))
    external=[]
    for name,meta in read(OUT/'github_sources.json').items():
        p=OUT/name;data=p.read_bytes();blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        assert blob==meta['gitBlobSHA1'],(name,blob)
        s=data.decode('utf8');anchors={'petsc_fas.c':['SNESLineSearchApply(fas->coarseCorrectionLineSearch','does not support bounds'],
          'petsc_linesearchbt.c':['MatMult(jac, Y, W)','lambda *= .5'],
          'petsc_linesearchsecant.c':['lambda     = .5 * (lambda + lambda_old)','SNES_LINESEARCH_FAILED_REDUCT']}[name]
        for token in anchors:assert token in s
        lines=s.splitlines();excerpts=[]
        for needle in anchors:
            i=next(i for i,line in enumerate(lines) if needle in line);excerpts.append(dict(line=i+1,text='\n'.join(lines[max(0,i-3):i+4])))
        external.append(dict(file=name,sha256=sha(p),excerpts=excerpts,**meta))
    save(OUT/'report.json',dict(passed=True,actualNativeRecoverySucceededBeforeObserverAbort=True,
      sourceBindings=source,externalSourceReview=external,
      adoptedNext='retain verified native legal-old-state recovery; passive recovered=true observation and whole-cycle accounting before any qualification',
      notAdoptedNow='installing PETSc or adding custom prolongation damping/line search; not needed to remove observer-only abort',
      conditionalAlternative='PETSc FAS coarse correction line search with SECANT; requires pinned build, side-effect-safe native residual adapter and domain-error semantics',
      nativeFallbackNotConservativeProof=True,actualPreProlongationStateMissing=True,postCallbackFullFlowSACompletionNotTested=True,
      terminalResidualCodeCompiledNotReached=True,newCFDAdmitted=False,newQualifiedAreaM2=0,physicalQualification=False))
    print(dict(passed=True,nativeSources=len(source),externalFiles=len(external),newCFDAdmitted=False))


def main():
    if CLOSED:raise SystemExit('v247 method decision closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',BASE/'audit03/actual_rejection.json',
      ROOT/'outputs/v235_validation/stage01/CNSVariable.cpp',ROOT/'outputs/v235_validation/stage01/dependencies.json',
      ROOT/'outputs/v245_validation/stage01/CMultiGridIntegration_observed.cpp',ROOT/'outputs/v245_validation/stage01/dependencies.json']
    locks += [p for p in OUT.iterdir() if p.is_file()]
    return run.setup_guard().supervise('method_review2',[sys.executable,str(Path(__file__).resolve()),'task'],30,1536,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v247 method decision closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':review()
    else:raise SystemExit(0 if main() else 1)
