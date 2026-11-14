"""模块职责：接受控制通过后，绑定唯一真实相态守卫分支；不扫描参数、不自动重试。"""
from pathlib import Path
import os,shutil,sys
import runPhaseAcceptanceV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CASE=BASE/'rocket01';BUILD=run.BUILD;CLOSED=True
CHECKPOINT=ROOT/'outputs/v222_validation/rocket01/restart_flow_000199.csv'
MESH=ROOT/'outputs/v222_validation/rocket01/mesh.su2'
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def prepare():
    """模块一：固定原199初态、原网格和v236收益合同；仅替换已测接受机制。"""
    for folder in ('process_build2','process_control','process_audit','audit01'):
        assert read(BASE/folder/'report.json')['passed'],folder
    assert read(BASE/'audit01/report.json')['realEvents']==466
    old=ROOT/'outputs/v236_validation';contract=read(old/'contract.json')
    assert sha(CHECKPOINT)==contract['sourceCheckpointSHA256'] and sha(MESH)==contract['meshSHA256'] and sha(GEOMETRY)==contract['geometrySHA256']
    bound=read(BASE/'stage02/inputs.json')
    for n,h in bound.items():assert sha(n)==h,n
    built=read(BUILD/'report.json');assert sha(BUILD/'driver.exe')==built['driverSHA256']
    CASE.mkdir(exist_ok=False)
    for name in ('case.cfg','original_faces.txt'):shutil.copyfile(old/'rocket01'/name,CASE/name)
    contract.update(branch='v222-state199-native-reset-v237-phase-feasibility',sourceObjectsRecompiled=1,
        buildSeconds=180,acceptanceControlSeconds=60,acceptanceGuard='native upper + common five-state feasibility halving',
        guardControlsPassed=True,controlReportSHA256=sha(BASE/'audit01/report.json'),
        numericalConfigurationByteIdenticalToV236=True,productionCallsitePreviouslyExecuted=False,
        forcedFinalNativeRestartRequired=True,rawIncrementObserver=True,automaticRetries=0)
    save(BASE/'contract.json',contract)
    for archive,rows in read(BUILD/'members.json').items():
        for row in rows:assert sha(row['path'])==row['sha256'];bound[row['path']]=row['sha256']
        bound[str(BUILD/archive)]=sha(BUILD/archive)
    for p in (BUILD/'driver.exe',BASE/'contract.json',BASE/'audit01/report.json',BASE/'control01/report.json',
              CASE/'case.cfg',CASE/'original_faces.txt',CHECKPOINT,MESH,GEOMETRY,Path(__file__).resolve()):bound[str(p)]=sha(p)
    save(BASE/'candidate_inputs.json',bound)
    save(BASE/'candidate_admission.json',dict(passed=True,unique=True,admittedBranch=contract['branch'],
        contractSHA256=sha(BASE/'contract.json'),inputsSHA256=sha(BASE/'candidate_inputs.json'),
        originalSurfaceAndPhysicalGatesPreserved=True,exactHistoricalContinuation=False,physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v237 real branch closed; follow NEXT_DECISION')
    run.configure()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseAcceptanceV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage=='prepare':
        locks += [BASE/'audit01/report.json',BASE/'stage02/inputs.json',BUILD/'report.json',BUILD/'members.json',CHECKPOINT,MESH,GEOMETRY]
        return run.guard.supervise('candidate_prepare',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE)
    assert stage=='run' and read(BASE/'candidate_admission.json')['passed']
    assert not (BASE/'process_rocket').exists(),'Unique real candidate already consumed'
    bound=read(BASE/'candidate_inputs.json')
    for n,h in bound.items():assert sha(n)==h,n
    locks+=list(map(Path,bound))+[BASE/'candidate_inputs.json',BASE/'candidate_admission.json']
    return run.guard.supervise('rocket',[str(BUILD/'driver.exe'),'case.cfg','.',os.path.relpath(str(CHECKPOINT),str(CASE)),'original_faces.txt'],
        2400,1024,list(dict.fromkeys(locks)),CASE,'cfd',True)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v237 real branch closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
