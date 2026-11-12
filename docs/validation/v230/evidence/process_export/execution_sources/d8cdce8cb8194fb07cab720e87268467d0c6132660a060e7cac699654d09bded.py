"""模块职责：准确绑定成熟原生循环，准备并监督唯一真实全局外流候选。"""
from pathlib import Path
import shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_SOURCE,OLD_BUILD
BASE=ROOT/'outputs/v230_validation';DOC=ROOT/'docs/validation/v230';CASE=BASE/'rocket01';BUILD=BASE/'build01'
PREV=ROOT/'outputs/v229_validation';HELPER=ROOT/'native/su2PhaseAwareAir/rocketEulerRun.cpp'
CLOSED=False

def prepare(stage='prepare01'):
    out=BASE/stage;out.mkdir(exist_ok=False);CASE.mkdir(exist_ok=True);BUILD.mkdir(exist_ok=True)
    assert read(PREV/'restart_audit01/report.json')['wholeNativeInterfacePassed']
    assert read(PREV/'build_audit02/report.json')['nativeControlAdmission']
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    olddeps=read(PREV/'build_audit02/actual_dependencies.json')
    for n,h in olddeps.items():assert sha(n)==h,n
    sources={}
    names=['SU2_CFD/src/drivers/CDriver.cpp','SU2_CFD/src/drivers/CSinglezoneDriver.cpp',
      'SU2_CFD/src/iteration/CFluidIteration.cpp','SU2_CFD/include/iteration/CIteration.hpp',
      'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl']
    for rel in names:
        expected=olddeps.get(str(PREV/'source'/rel),tree[str(OLD_SOURCE/rel)])
        candidates=[PREV/'source'/rel,ROOT/'outputs/v222_validation/source01'/rel,
          ROOT/'outputs/v221_validation/source01'/rel,ROOT/'outputs/v167_validation/build01/executed_sources'/rel,OLD_SOURCE/rel]
        p=next((p for p in candidates if p.is_file() and sha(p)==expected),None);assert p,rel
        target=BASE/'source01'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        sources[rel]=dict(path=str(p),sha256=expected)
    save(out/'sources.json',sources)
    mesh=ROOT/'outputs/v194_validation/rocket01/mesh.su2';assert sha(mesh)=='faff7bf901d330eaf53a4616014620281f54e58e8137d1a0e4943ab14804f0e7'
    if (CASE/'mesh.su2').exists():assert sha(mesh)==sha(CASE/'mesh.su2')
    else:shutil.copyfile(mesh,CASE/'mesh.su2')
    cfg={}
    for line in (ROOT/'outputs/v228_validation/prepare03/case.proposed.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(MESH_FILENAME='mesh.su2',KIND_TURB_MODEL='NONE',ITER='800',CONV_STARTITER='10',
      OUTPUT_FILES='(RESTART_ASCII)',SCREEN_WRT_FREQ_INNER='10',OUTPUT_WRT_FREQ='50')
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    cmd=read(PREV/'build02/helper_command.json')
    cmd=[str(HELPER) if x==str(ROOT/'native/su2PhaseAwareAir/eulerCompatibility.cpp') else
      str(BUILD/'rocket.d') if x==str(PREV/'build02/compatibility.d') else
      str(BUILD/'rocket.exe') if x==str(PREV/'build02/compatibility.exe') else x for x in cmd]
    filtered=[];i=0
    while i<len(cmd):
        if cmd[i] in ('-MF','-o'):i+=2;continue
        if cmd[i] in ('-MD','-Wl,--start-group','-Wl,--end-group') or cmd[i].endswith('.a'):i+=1;continue
        filtered.append(cmd[i]);i+=1
    dep=out/'helper.d';dep_cmd=filtered+['-M','-MT','dep','-MF',str(dep)]
    subprocess.check_call(dep_cmd,cwd=str(OLD_BUILD),timeout=20)
    bindings={}
    for n in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
        if p!=HELPER:
            expected=olddeps.get(str(p))
            if expected is None:
                rel=p.relative_to(PREV/'source');expected=tree[str(OLD_SOURCE/rel)]
            assert expected==h,str(p)
        bindings[str(p)]=h
    save(BUILD/'dependencies.json',bindings);save(BUILD/'command.json',cmd);save(out/'dependency_command.json',dep_cmd)
    shutil.copyfile(HELPER,out/'executed_helper.cpp')
    save(CASE/'contract.json',dict(model='rah-su2-phase-air-euler-eplus300k-v1',caseRole='finite-unconverged-research',
      meshSHA256=sha(CASE/'mesh.su2'),configSHA256=sha(CASE/'case.cfg'),updatesLimit=800,wallSecondsLimit=600,
      physicalResidual='-nativeImplicitRHS; physical energy component minus 300000 mass; before RMS',
      residualStage='evaluated for update; differs from saved post-update conservative state',
      geometryIdentity='v194-original-19752-wall-faces',firstTargetFaces=3128,heatFluxSolved=False,
      physicalQualification=False,expectedStateUpdates=[0,1]+list(range(50,801,50)),
      absoluteLogRMSLimit=-8,areaPressureLastPairLimit=.01,massEnergyRelativeBalanceLimit=1e-8,
      spatialFarfieldTested=False,automaticRetries=0))
    save(out/'report.json',dict(passed=True,dependencies=len(bindings),actualSources=len(sources),CFDStarted=False,newMeshGenerated=False))

def inputs():
    result=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',HELPER,
      PREV/'build02/members.json',PREV/'build_audit02/actual_dependencies.json',PREV/'restart_audit01/report.json',
      ROOT/'outputs/v228_validation/prepare03/case.proposed.cfg',ROOT/'outputs/v194_validation/rocket01/mesh.su2']
    for n in read(PREV/'build_audit02/actual_dependencies.json'):result.append(Path(n))
    for rows in read(PREV/'build02/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256'];result.append(Path(row['path']))
    if (BASE/'prepare01/sources.json').exists():
        for row in read(BASE/'prepare01/sources.json').values():result.append(Path(row['path']))
    return list(dict.fromkeys(result))

def main(stage):
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    locks=inputs();kind='offline';seconds=120;cwd=BASE
    if stage in ('prepare','prepare2'):cmd=[sys.executable,str(Path(__file__).resolve()),'task_'+stage]
    elif stage=='build':
        assert read(BASE/'prepare02/report.json')['passed'];cmd=read(BUILD/'command.json');kind='build'
        locks += [Path(n) for n in read(BUILD/'dependencies.json')]+[BUILD/'command.json']
    elif stage=='run':
        assert read(BASE/'process_build/report.json')['passed'];assert not (BASE/'process_run').exists()
        contract=read(CASE/'contract.json');assert sha(CASE/'mesh.su2')==contract['meshSHA256'];assert sha(CASE/'case.cfg')==contract['configSHA256']
        kind='cfd';seconds=600;cwd=CASE;cmd=[str(BUILD/'rocket.exe'),'case.cfg','.']
        locks += [BUILD/'rocket.exe',CASE/'contract.json',CASE/'case.cfg',CASE/'mesh.su2']
    else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise(stage,cmd,seconds,512,list(dict.fromkeys(locks)),cwd,kind,stage=='run')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v230 finite candidate closed; no automatic continuation')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_prepare2':prepare('prepare02')
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
