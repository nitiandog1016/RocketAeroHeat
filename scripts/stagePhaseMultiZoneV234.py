"""模块职责：独立两黏性区身份与最小注册分派，不修改旧单区合同、气体或数值公式。"""
from pathlib import Path
import difflib,shutil,subprocess,sys,time
import runMultiZoneCompatibilityV234 as base
from observeMultiZoneFailureV234 import safe_wait_procs
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=base.BASE;DOC=base.DOC;SOURCE=BASE/'mz_source';BUILD=BASE/'mz_build01';PREP=BASE/'mz_prepare01'
CLOSED=True
LINK=BASE/'mz_link02'
GUARD='''  // 独立双黏性域身份；旧单区RAH_PHASE_AIR原合同保持。
  if (Kind_FluidModel == RAH_PHASE_AIR_MZ) {
    if (SystemMeasurements != SI || Ref_NonDim != DIMENSIONAL || Time_Domain ||
        !Multizone_Problem || nZone != 2 || DiscreteAdjoint || ContinuousAdjoint ||
        Kind_Solver != MAIN_SOLVER::RANS || Kind_Turb_Model != TURB_MODEL::SA ||
        Kind_Upwind_Flow != UPWIND::ROE || nMarker_WallFunctions ||
        nMarker_Isothermal == 0 || nMarker_Riemann == 0 || nMarker_Fluid_InterfaceBound != 1 ||
        nMarker_Euler || nMarker_HeatFlux || nMarker_HeatTransfer ||
        nMarker_FarField || nMarker_Supersonic_Inlet || nMarker_Supersonic_Outlet ||
        nMarker_Inlet || nMarker_Outlet || Low_Mach_Corr || nMarker_Giles || GetDynamic_Grid())
      SU2_MPI::Error("RAH_PHASE_AIR_MZ requires direct SI steady two-zone RANS/SA, ROE, resolved thermal wall and one fluid interface", CURRENT_FUNCTION);
    for (unsigned short i=0; i<nMarker_Riemann; ++i)
      if (Kind_Data_Riemann[i] != STATIC_SUPERSONIC_INFLOW_PT)
        SU2_MPI::Error("RAH_PHASE_AIR_MZ requires ambient PT characteristic boundary", CURRENT_FUNCTION);
  }

'''


def transform(rel,text):
    if rel=='Common/include/option_structure.hpp':
        a='  RAH_PHASE_AIR = 15, // 相态标准Air，内部能量+300000J/kg';assert text.count(a)==1
        text=text.replace(a,'  RAH_PHASE_AIR_MZ = 17, // 独立双黏性域；16保留给已存在Euler身份\n'+a)
        a='  MakePair("RAH_PHASE_AIR", RAH_PHASE_AIR)';assert text.count(a)==1
        return text.replace(a,'  MakePair("RAH_PHASE_AIR_MZ", RAH_PHASE_AIR_MZ)\n'+a)
    if rel=='Common/src/CConfig.cpp':
        anchor='  /*--- Check if CoolProp is used with non-dimensionalization. ---*/';assert text.count(anchor)==1
        text=text.replace(anchor,GUARD+anchor)
        anchor='if (Kind_FluidModel == RAH_PHASE_AIR &&';assert text.count(anchor)==1
        return text.replace(anchor,'if ((Kind_FluidModel == RAH_PHASE_AIR || Kind_FluidModel == RAH_PHASE_AIR_MZ) &&')
    if rel=='SU2_CFD/src/solvers/CEulerSolver.cpp':
        assert text.count('case RAH_PHASE_AIR:')==3
        text=text.replace('case RAH_PHASE_AIR:','case RAH_PHASE_AIR_MZ:\n    case RAH_PHASE_AIR:')
        a='config->GetKind_FluidModel() == RAH_PHASE_AIR';assert text.count(a)==3
        return text.replace(a,'('+a+' || config->GetKind_FluidModel() == RAH_PHASE_AIR_MZ)')
    if rel=='SU2_CFD/src/numerics/flow/convection/roe.cpp':
        a='config->GetKind_FluidModel() == RAH_PHASE_AIR';assert text.count(a)==1
        return text.replace(a,'('+a+' || config->GetKind_FluidModel() == RAH_PHASE_AIR_MZ)')
    raise ValueError(rel)


def prepare():
    PREP.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    known=read(base.GAS/'build02/dependencies.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    changed=['Common/include/option_structure.hpp','Common/src/CConfig.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/numerics/flow/convection/roe.cpp']
    copied={};changes={};patch=[]
    for part in ('Common/include','SU2_CFD/include'):
        for p in (base.GAS/'source'/part).rglob('*'):
            if not p.is_file():continue
            rel=p.relative_to(base.GAS/'source').as_posix();h=sha(p);assert h==known.get(str(p),tree.get(str(base.OLD_SOURCE/rel))),rel
            target=SOURCE/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target);copied[str(target)]=sha(target)
    for rel in changed:
        p=base.GAS/'source'/rel;h=sha(p);assert h==known[str(p)],rel
        before=p.read_text(encoding='utf8');after=transform(rel,before);target=SOURCE/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(after,encoding='utf8')
        patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='old/'+rel,tofile='phase-mz/'+rel))
        copied[str(target)]=sha(target);changes[rel]=dict(source=str(p),originalSHA256=h,patchedSHA256=sha(target))
    subprocess.check_call(['cmd','/c','mklink','/J',str(SOURCE/'externals'),str(base.OLD_SOURCE/'externals')],stdout=subprocess.DEVNULL)
    (PREP/'phase_multi_zone.patch').write_text(''.join(patch),encoding='utf8');save(PREP/'changes.json',changes)
    commands={}
    for rel,cmd in read(base.GAS/'prepare03/commands.json').items():
        if rel not in changed:continue
        cmd=[x.replace(str(base.GAS/'source'),str(SOURCE)).replace(str(base.GAS/'build'),str(BUILD)) for x in cmd]
        obj=cmd[cmd.index('-o')+1];commands[rel]=cmd[:-2]+['-MD','-MF',obj+'.d']+cmd[-2:]
    helper=read(BASE/'build03/command.json')
    helper=[x.replace(str(base.GAS/'source'),str(SOURCE)).replace(str(base.GAS/'build02'),str(BUILD)).replace(str(BASE/'build03'),str(BUILD)) for x in helper]
    commands['helper']=helper;save(PREP/'commands.json',commands)
    bindings={}
    for j,(name,cmd) in enumerate(commands.items()):
        filtered=[];i=0
        while i<len(cmd):
            x=cmd[i]
            if x in ('-MF','-o'):i+=2;continue
            if x in ('-MD','-c','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            filtered.append(x);i+=1
        path=PREP/('unit%d.d'%j);subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(path)],cwd=str(base.OLD_BUILD),timeout=20)
        for n in path.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(n);p=(p if p.is_absolute() else base.OLD_BUILD/p).resolve();h=sha(p)
            if str(p) in copied:assert h==copied[str(p)]
            elif p!=base.CPP:assert h==known.get(str(p),tree.get(str(p))),str(p)
            bindings[str(p)]=h
    save(PREP/'dependencies.json',bindings);save(PREP/'copied.json',copied)
    case=BASE/'multi03';case.mkdir()
    for p in (BASE/'multi02').iterdir():
        if p.suffix in ('.cfg','.su2') or p.name=='mesh_mapping.json':shutil.copyfile(p,case/p.name)
    for z in (0,1):
        p=case/('zone%d.cfg'%z);t=p.read_text(encoding='utf8');assert t.count('FLUID_MODEL= RAH_PHASE_AIR\n')==1
        p.write_text(t.replace('FLUID_MODEL= RAH_PHASE_AIR\n','FLUID_MODEL= RAH_PHASE_AIR_MZ\n'),encoding='utf8')
    save(PREP/'report.json',dict(passed=True,identity='rah-su2-phase-air-mz-saneg-eplus300k-v1',enumValue=17,modifiedSourceFiles=4,
      compiledObjectsPlanned=3,numericalFormulaChanges=0,oldSingleZoneGuardPreserved=True,dependencies=len(bindings),flowUpdates=0))


def build():
    commands=read(PREP/'commands.json');receipts=[]
    for name,cmd in commands.items():
        if name=='helper':continue
        start=time.monotonic();subprocess.check_call(cmd,cwd=str(base.OLD_BUILD));receipts.append(dict(unit=name,seconds=time.monotonic()-start));save(BUILD/'compile_receipts.json',receipts)
    ar=str(Path(commands['helper'][0]).parent/'ar.exe');members=read(base.GAS/'build02/members.json');new={}
    for name,rows in members.items():
        result=[]
        for r in rows:
            p=Path(r['path']);assert sha(p)==r['sha256']
            if p.name in ('CConfig.obj','CEulerSolver.obj','roe.obj'):p=BUILD/p.name
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result],cwd=str(base.OLD_BUILD));new[name]=result
    save(BUILD/'members.json',new);subprocess.check_call(commands['helper'],cwd=str(base.OLD_BUILD))
    save(BUILD/'report.json',dict(passed=True,helperSHA256=sha(BUILD/'control.exe'),physicsObjectsUnchangedExceptDispatch=3,
      fullEngineBuilt=False,flowUpdates=0))


def prepare4():
    """只修正多区配置语法；保留首初始化拒绝及原小壳。"""
    case=BASE/'multi04';case.mkdir(exist_ok=False)
    for p in (BASE/'multi03').iterdir():
        if p.suffix in ('.cfg','.su2') or p.name=='mesh_mapping.json':shutil.copyfile(p,case/p.name)
    for z in (0,1):
        p=case/('zone%d.cfg'%z);s=p.read_text(encoding='utf8');assert s.count('\nITER= 1\n')==1
        p.write_text(s.replace('\nITER= 1\n','\nINNER_ITER= 1\n'),encoding='utf8')
    save(case/'configuration_recovery.json',dict(onlyChange='zone ITER=1 to INNER_ITER=1',flowUpdates=0))


def main(stage):
    if CLOSED:raise SystemExit('v234 multi-zone identity closed; use NEXT_DECISION')
    g=base.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda _:size(BASE)
    g.psutil.wait_procs=safe_wait_procs(g.psutil.wait_procs)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/runMultiZoneCompatibilityV234.py',
      ROOT/'scripts/observeMultiZoneFailureV234.py',DOC/'PLAN.md',base.CPP,base.GAS/'build02/dependencies.json',base.GAS/'build02/members.json']
    for group in read(base.GAS/'build02/members.json').values():
        for row in group:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    if stage=='prepare_mz':
        locks += [Path(n) for n in read(base.GAS/'build02/dependencies.json')]
        return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,128,list(dict.fromkeys(locks)),BASE)
    deps=read(PREP/'dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    locks += [Path(n) for n in deps]+[PREP/'commands.json']
    if stage in ('link_mz','multi3','prepare4','multi4','capture4'):
        # 三对象已有成功编译记录；只恢复薄驱动链接，首次容量失败原样保留。
        failed=read(BASE/'process_build_mz/report.json')
        assert failed['stopReason']=='output' and failed['inputsStable'] and not failed['liveGroupAfterCleanup']
        assert len(read(BUILD/'compile_receipts.json'))==3
        for rows in read(BUILD/'members.json').values():
            for row in rows:
                assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
        locks += [BUILD/'members.json',BUILD/'libSU2core.a',BUILD/'libSU2Common.a',BUILD/'compile_receipts.json']
    if stage=='link_mz':
        LINK.mkdir(exist_ok=False)
        cmd=read(PREP/'commands.json')['helper']
        cmd=[x.replace(str(BUILD/'control.exe'),str(LINK/'control.exe')).replace(str(BUILD/'control.d'),str(LINK/'control.d')) for x in cmd]
        save(LINK/'command.json',cmd);locks.append(LINK/'command.json')
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='identity-build')
        return g.supervise(stage,cmd,180-used,256,list(dict.fromkeys(locks)),base.OLD_BUILD,'identity-build')
    if stage=='build_mz':
        return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build'],180,128,list(dict.fromkeys(locks)),BASE,'identity-build')
    if stage=='prepare4':
        return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare4'],120,256,list(dict.fromkeys(locks)),BASE)
    if stage=='capture4':
        # 唯一同输入崩溃定位；调试器停止后直接退出，不继续或推进流动。
        case=BASE/'capture04';case.mkdir(exist_ok=False)
        for p in (BASE/'multi04').iterdir():
            if p.suffix in ('.cfg','.su2'):shutil.copyfile(p,case/p.name)
        gdb=base.OLD_SOURCE.parent/'toolchain/w64devkit/bin/gdb.exe'
        locks += [gdb,LINK/'control.exe']+list(case.iterdir())
        cmd=[str(gdb),'-batch','-ex','set pagination off','-ex','run','-ex','bt 12','-ex','info args',
          '-ex','p failure_observer::phase','-ex','p failure_observer::scenario','-ex','p failure_observer::zone',
          '--args',str(LINK/'control.exe'),'case.cfg','.','multi']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='recovery-query')
        return g.supervise(stage,cmd,60-used,256,list(dict.fromkeys(locks)),case,'recovery-query')
    if stage in ('multi3','multi4'):
        assert read(BASE/'process_link_mz/report.json')['passed']
        case=BASE/('multi03' if stage=='multi3' else 'multi04');locks += [LINK/'control.exe']+[p for p in case.iterdir() if p.is_file()]
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='recovery-query')
        return g.supervise(stage,[str(LINK/'control.exe'),'case.cfg','.','multi'],60-used,256,list(dict.fromkeys(locks)),case,'recovery-query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v234 multi-zone identity closed; use NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    elif sys.argv[1]=='task_prepare4':prepare4()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
