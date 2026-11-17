"""模块职责：固定方向诊断的准确输入/薄编/完整域资源监督；不重开旧KINSOL。"""
import sys,shutil,subprocess
from pathlib import Path
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v258_validation';DOC=ROOT/'docs/validation/v258';PREV=ROOT/'outputs/v257_validation'
CASE=BASE/'rocket01';BUILD=BASE/'build01';STAGE=BASE/'stage01';CLOSED=False
DRIVER=ROOT/'native/su2PhaseAwareAir/directionConsistencyV258.cpp';HEADER=ROOT/'native/su2PhaseAwareAir/FullResidualProviderV257.hpp'
WD=ROOT/'tools/runtime/su2-thermal-build/build-v117';SRC=ROOT/'tools/runtime/su2-thermal-build/source'


def prepare():
    for p in [CASE,BUILD,STAGE]:p.mkdir(exist_ok=False)
    assert read(PREV/'process_control/report.json')['passed']
    n=read(PREV/'rocket01/native_completion.json');assert n['nativeFlag']==-13 and n['linearGateRefused'] and n['observedCurrentStateChanges']==0
    for key in ['initial_unknowns.bin','linear_refused_trial.bin','initial_state.bin','initial_full_residual.bin','scaling.csv']:
        index=read(ROOT/'docs/validation/v257/evidence/local_asset_index.json');assert sha(PREV/'rocket01'/key)==index['rocket01/'+key]['sha256']
    shutil.copyfile(PREV/'rocket01/case.cfg',CASE/'case.cfg')
    prior=read(PREV/'stage02/bridge_command.json');cmd=[]
    for part in prior:
        if 'libsundials_kinsol.a' in part or 'libsundials_sunlinsolspgmr.a' in part:continue
        if part==str(ROOT/'native/su2PhaseAwareAir/kinsolCoupledV257.cpp'):part=str(DRIVER)
        if part in [str(PREV/'build01/preflight.exe'),str(PREV/'build01/preflight.d')]:part=part.replace(str(PREV/'build01'),str(BUILD))
        cmd.append(part)
    save(STAGE/'command.json',cmd)
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    sources={}
    for name in ['Common/src/linear_algebra/CSysMatrix.cpp','SU2_CFD/src/solvers/CTurbSASolver.cpp']:
        p=SRC/name;assert sha(p)==tree[str(p)];sources[str(p)]=sha(p)
    ns=ROOT/'outputs/v220_validation/source/SU2_CFD/src/solvers/CNSSolver.cpp';assert sha(ns)==read(ROOT/'outputs/v220_validation/build02/dependencies.json')[str(ns)];sources[str(ns)]=sha(ns)
    sundials=ROOT/'tools/runtime/sundials-7.4.0-8e17876/source';manifest=read(PREV/'stage02/library_source.json')
    for name in ['src/kinsol/kinsol_ls.c','src/sunlinsol/spgmr/sunlinsol_spgmr.c']:
        p=sundials/name;assert sha(p)==manifest['files'][name];sources[str(p)]=sha(p)
    components=dict(nativeReference=1549975552,ILU=640513952,adapter10Vectors=157987520,native6BlockVectors=48344544,NVector5Vectors=78993760,reserve=256*1024**2)
    forecast=sum(components.values());assert forecast<4*1024**3
    save(STAGE/'report.json',dict(passed=True,baseCommit='29d7f95a2a01ba61cdf3db46799e88fab96f94ca',accurateSourceSHA256=sources,
        sourceFindings=dict(strongWallRowsIdentity=True,physicalEnergyTransformation='C: E_native=E_physical+300000*rho; C^-1 inverse on residual and direction',
            SPGMRRightScalingRemovedBeforeUserPSolve=True,SPGMRLeftScalingAfterATimes=True,blockPreconditionerOmitsFlowSACrossBlocks=True,
            ILUCopySeparateFromSpatialMatrix=True,noSpecificInterfaceBugProved=True),
        forecastComponentsBytes=components,forecastBytes=forecast,fullResidualCallLimit=4,KINSOLExecuted=False,physicalQualification=False))
    # 定向源码摘录供下轮恢复，不复制成熟源码整树。
    ranges={'Common/src/linear_algebra/CSysMatrix.cpp':[(922,945),(1080,1101),(1229,1261)],
        'SU2_CFD/src/solvers/CTurbSASolver.cpp':[(540,572)]}
    chunks=[]
    for name,spans in ranges.items():
        lines=(SRC/name).read_text(encoding='utf8').splitlines()
        for first,last in spans:chunks.append(name+f':{first}-{last}\n'+'\n'.join(f'{i+1}: {lines[i]}' for i in range(first-1,last)))
    (STAGE/'source_excerpt.txt').write_text('\n\n'.join(chunks)+'\n',encoding='utf8')


def build():
    subprocess.check_call(read(STAGE/'command.json'),cwd=WD)
    deps={}
    for name in (BUILD/'preflight.d').read_text(encoding='utf8').replace('\\\n',' ').split(': ',1)[1].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();deps[str(p)]=sha(p)
    save(BUILD/'dependencies.json',deps);save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'preflight.exe'),physicsObjectsRecompiled=0,matureLibrariesRecompiled=0))


def main(stage):
    if CLOSED:raise SystemExit('v258 fixed direction diagnostic closed; follow NEXT_DECISION')
    assert stage in ['prepare','build','control'];BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',DRIVER,HEADER,ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        PREV/'process_control/report.json',PREV/'rocket01/native_completion.json',PREV/'rocket01/case.cfg',PREV/'stage02/bridge_command.json',
        PREV/'stage02/library_source.json',PREV/'stage02/library_build.json',PREV/'build01/members.json',ROOT/'docs/validation/v257/evidence/local_asset_index.json',
        ROOT/'outputs/v167_validation/build01/source_tree.json',ROOT/'outputs/v220_validation/build02/dependencies.json']
    locks += [PREV/'rocket01'/p for p in ['initial_unknowns.bin','linear_refused_trial.bin','initial_state.bin','initial_full_residual.bin','scaling.csv']]
    dependencies=read(PREV/'stage02/bridge_dependencies.json')
    for p,h in dependencies.items():assert sha(p)==h,p;locks.append(Path(p))
    for rows in read(PREV/'build01/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    for row in read(PREV/'stage02/library_build.json')['staticLibraries'].values():assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    if stage=='prepare':argv=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];seconds=60;kind='offline';cwd=ROOT
    else:
        assert read(STAGE/'report.json')['passed'];locks += [STAGE/'report.json',STAGE/'command.json',CASE/'case.cfg']
        for p,h in read(STAGE/'report.json')['accurateSourceSHA256'].items():assert sha(p)==h,p;locks.append(Path(p))
        if stage=='build':argv=[sys.executable,str(Path(__file__).resolve()),'task_build'];seconds=120;kind='build';cwd=ROOT
        else:
            assert read(BUILD/'report.json')['passed'] and sha(BUILD/'preflight.exe')==read(BUILD/'report.json')['engineSHA256']
            for p,h in read(BUILD/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
            locks += [BUILD/'report.json',BUILD/'preflight.exe',BUILD/'dependencies.json',ROOT/'outputs/v250_validation/rocket01/mesh.su2',ROOT/'outputs/v250_validation/rocket01/input_restart.csv']
            argv=[str(BUILD/'preflight.exe'),'case.cfg','.'];seconds=180;kind='control';cwd=CASE
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p==Path(__file__).resolve() or p in [DOC/'PLAN.md',DRIVER,HEADER]:shutil.copyfile(p,snap/(sha(p)+p.suffix))
    return supervise(guard,BASE,stage,argv,seconds,512,list(dict.fromkeys(locks)),cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v258 fixed direction diagnostic closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
