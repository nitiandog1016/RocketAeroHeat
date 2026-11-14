"""模块职责：复用本轮事件与准确成熟通量源码，形成具名下一候选；不再运行原生或CFD。"""
from pathlib import Path
import json,sys
import runPhaseBranchV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC
TREE=ROOT/'outputs/v167_validation/build01/source_tree.json'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
FILES=[SOURCE/'SU2_CFD/src/drivers/CDriver.cpp',SOURCE/'SU2_CFD/src/numerics/flow/convection/hllc.cpp',
       ROOT/'outputs/v220_validation/source/Common/src/CConfig.cpp']


def main():
    out=BASE/'outcome01';out.mkdir(exist_ok=False);tree=read(TREE);known=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    bindings=[]
    windows=[[(1739,1757)],[(618,640),(672,714),(1077,1102)],[(3808,3818),(4612,4619)]]
    for p,spans in zip(FILES,windows):
        expected=tree.get(str(p),known.get(str(p)));assert expected and sha(p)==expected,str(p)
        lines=p.read_text(encoding='utf8').splitlines();parts=[]
        for a,b in spans:parts.append('\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,min(b,len(lines))+1)))
        (out/(p.stem+'_excerpts.txt')).write_text('\n\n'.join(parts)+'\n',encoding='utf8');bindings.append(dict(path=str(p),sha256=expected,lines=spans))
    actual=read(BASE/'production_audit02/report.json');a=read(BASE/'real_audit01/report.json')
    events=read(BASE/'production_audit02/limited_state_conditions.json')['events']
    persistent=[r for r in events if r['globalID']==335030]
    assert persistent[-1]['update']==92 and (BASE/'rocket01/failure.txt').read_text().strip().endswith('what=PHASE_STEP_NO_ADMISSIBLE_PROGRESS')
    members=read(BASE/'build01/members.json');objects=[r for rows in members.values() for r in rows if 'convection_hllc.cpp.obj' in r['path'] or 'drivers_CDriver.cpp.obj' in r['path']]
    for r in objects:assert sha(r['path'])==r['sha256']
    save(out/'report.json',dict(passed=True,sourceBindings=bindings,matureObjects=objects,
        failure='PHASE_STEP_NO_ADMISSIBLE_PROGRESS',completedUpdates=92,attemptedUpdate=93,failedPointIDRecorded=False,
        lastPersistentlyLimitedPoint=persistent[-1],candidateIsNotProofOfFailedPoint=True,
        replayPassedButRealDevelopmentFailed=True,rawLimitedEvents=actual['limitedCandidates'],oldStateRecoveries=0,
        independentActualSavedLimitedComponents=0,allLimitedCandidatesCheckedButNotAllActualUpdatedStatesSaved=True,
        nativeFinalForceWritingReached=False,fullFinalFramePresent=False,latestFullSavedUpdate=1,latestSurfaceSavedUpdate=50,
        lastCompletePhysicalLogRMS=a['lastCompletePhysicalLogRMS'],sourceConfirmedGeneralHLLCFactory=True,
        existingPhaseConfigRejectsHLLC=True,generalHLLCUsesRoeAveragingAndHasNoProvenStrictPhaseDomainGuarantee=True,
        selectedNextScope='bounded general-EOS HLLC compatibility decision using existing objects and a single whole-interface control, then one separately admitted branch if justified',
        genericFluxOrJacobianReimplementationAllowed=False,currentCFDReadmission=False,physicalQualification=False))
    print(json.dumps(dict(passed=True,failedUpdate=93,nextCandidate='native general-EOS HLLC; not admitted',phaseGuardAloneInsufficient=True)))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseBranchV237.py',ROOT/'scripts/runPhaseAcceptanceV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',TREE,
          ROOT/'outputs/v220_validation/build02/dependencies.json',BASE/'build01/members.json',BASE/'production_audit02/report.json',
          BASE/'production_audit02/limited_state_conditions.json',BASE/'real_audit01/report.json',BASE/'rocket01/failure.txt']+FILES
        raise SystemExit(0 if run.run.guard.supervise('outcome',[sys.executable,str(Path(__file__).resolve()),'task'],30,1024,locks,BASE) else 1)
