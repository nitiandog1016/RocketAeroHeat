"""模块职责：固定两层FAS准备、单对象构建与完整原生控制，逐任务监督并保留失败。"""
from pathlib import Path
import difflib,json,posixpath,re,shutil,subprocess,sys,time
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v245_validation';DOC=ROOT/'docs/validation/v245';CLOSED=True
GAS=ROOT/'outputs/v220_validation';OLD=ROOT/'tools/runtime/su2-thermal-build/source'
WD=OLD.parent/'build-v117';PREV=ROOT/'outputs/v235_validation';NATIVE=ROOT/'native/su2PhaseAwareAir'
STAGE=BASE/'stage01';BUILD=BASE/'build01'


def prepare():
    """模块一：绑定已有执行来源；只插只读回调，原算法字句不改。"""
    STAGE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    known=read(GAS/'build02/dependencies.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    specs=['Common/src/CConfig.cpp','Common/src/geometry/CMultiGridGeometry.cpp',
        'SU2_CFD/src/integration/CMultiGridIntegration.cpp','SU2_CFD/src/integration/CSingleGridIntegration.cpp',
        'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/iteration/CFluidIteration.cpp','SU2_CFD/include/solvers/CSolver.hpp']
    bindings=[]
    for rel in specs:
        expected=known.get(str(GAS/'source'/rel),tree.get(str(OLD/rel)))
        choices=[PREV/'source01'/rel,GAS/'source'/rel,ROOT/'outputs/v167_validation/build01/executed_sources'/rel,OLD/rel]
        p=next((p for p in choices if p.is_file() and sha(p)==expected),None);assert p,rel
        dest=STAGE/Path(rel).name;shutil.copyfile(p,dest)
        bindings.append(dict(path=str(p),sha256=expected,snapshot=str(dest),relative=rel))
    save(STAGE/'bindings.json',bindings)
    rel='SU2_CFD/src/integration/CMultiGridIntegration.cpp';p=STAGE/Path(rel).name
    original=p.read_text(encoding='utf8');s=original
    def includes(m):
        name=posixpath.normpath(posixpath.join(posixpath.dirname(rel),m.group(1)))
        q=GAS/'source'/name;assert q.is_file()
        assert sha(q)==known.get(str(q),tree.get(str(OLD/name))),name
        return '#include "'+q.as_posix()+'"'
    s=re.sub(r'#include "(\.\./[^\"]+)"',includes,s)
    s='#include "'+(NATIVE/'FASHooksV245.hpp').as_posix()+'"\n'+s
    def hook(label):
        return '\n    if(rah_fas_observer::transfer) rah_fas_observer::transfer("'+label+'",geometry_fine,solver_container_fine,geometry_coarse,solver_container_coarse,config);\n'
    anchors=[('    SetRestricted_Solution(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);','restrict_in','restrict_out'),
       ('    solver_coarse->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);',None,'coarse_preprocessed'),
       ('    SetForcing_Term(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1);','forcing_in','forcing_out'),
       ('    GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);','correction_in','correction_out'),
       ('    SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);','prolong_in','prolong_out')]
    for anchor,before,after in anchors:
        assert s.count(anchor)==1,anchor;s=s.replace(anchor,(hook(before) if before else '')+anchor+hook(after))
    patched=STAGE/'CMultiGridIntegration_observed.cpp';patched.write_text(s,encoding='utf8')
    (STAGE/'observer.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile=rel,tofile='observed/'+rel)),encoding='utf8')
    cmd=list(read(PREV/'stage01/commands.json')['SU2_CFD/src/integration/CIntegration.cpp'])
    cmd[cmd.index('-c')+1]=str(patched);cmd[cmd.index('-o')+1]=str(BUILD/'CMultiGridIntegration.obj');cmd[cmd.index('-MF')+1]=str(BUILD/'multi.d')
    helper=read(ROOT/'outputs/v236_validation/build01/command.json')
    helper=[str(NATIVE/'fasCompatibilityV245.cpp') if x.endswith('restartBranchV236.cpp') else
       x.replace(str(ROOT/'outputs/v236_validation/build01'),str(BUILD)).replace(str(PREV/'build01'),str(BUILD)) for x in helper]
    save(STAGE/'commands.json',dict(object=cmd,helper=helper))
    deps={}
    for j,command in enumerate((cmd,helper)):
        args=[];i=0
        while i<len(command):
            x=command[i]
            if x in ('-o','-MF'):i+=2;continue
            if x in ('-c','-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            args.append(x);i+=1
        d=STAGE/('deps%d.d'%j);subprocess.check_call(args+['-M','-MT','dep','-MF',str(d)],cwd=str(WD),timeout=25)
        for name in d.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            q=Path(name);q=(q if q.is_absolute() else WD/q).resolve();h=sha(q)
            if q.parent not in (STAGE,NATIVE):
                expect=known.get(str(q),tree.get(str(q)))
                if expect is None:expect=tree[str(OLD/q.relative_to(GAS/'source'))]
                assert h==expect,str(q)
            deps[str(q)]=h
    save(STAGE/'dependencies.json',deps)
    original_case=PREV/'observe02';cfg={}
    for line in (original_case/'case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(MGLEVEL='1',MGCYCLE='V_CYCLE',MG_MIN_MESHSIZE='1',MG_PRE_SMOOTH='(1, 1)',
      MG_POST_SMOOTH='(1, 0)',MG_CORRECTION_SMOOTH='(0, 0)',MG_DAMP_RESTRICTION='0.5',
      MG_DAMP_PROLONGATION='0.5',MG_SMOOTH_EARLY_EXIT='NO',MG_SMOOTH_OUTPUT='YES')
    for name in ('off01','observe01'):
        out=BASE/name;out.mkdir(exist_ok=False);shutil.copyfile(original_case/'mesh.su2',out/'mesh.su2')
        (out/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    save(STAGE/'report.json',dict(passed=True,dependencies=len(deps),numericalObjectsChanged=0,observedObjectsToCompile=1,
       updatesEach=5,expectedLayers=2,physicalQualification=False,publicBenchmarkExecuted=False))


def build():
    """模块二：仅编必要观察对象；成熟物性、边界、湍流与几何归档逐SHA复用。"""
    cmds=read(STAGE/'commands.json');subprocess.check_call(cmds['object'],cwd=str(WD))
    ar=str(Path(cmds['helper'][0]).parent/'ar.exe');members=read(PREV/'build01/members.json');changed=0;new={}
    for name,rows in members.items():
        result=[]
        for row in rows:
            p=Path(row['path']);assert sha(p)==row['sha256']
            if p.name=='integration_CMultiGridIntegration.cpp.obj':p=BUILD/'CMultiGridIntegration.obj';changed+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result],cwd=str(WD));new[name]=result
    assert changed==1;save(BUILD/'members.json',new);subprocess.check_call(cmds['helper'],cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,objectsCompiled=1,readOnlyCallbacks=True,engineSHA256=sha(BUILD/'branch.exe')))


def main(stage):
    if CLOSED:raise SystemExit('v245 fixed control closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
       NATIVE/'fasCompatibilityV245.cpp',NATIVE/'FASHooksV245.hpp',NATIVE/'SingleDomainHooksV235.hpp']
    if stage=='prepare':
        locks += [PREV/'observe02/case.cfg',PREV/'observe02/mesh.su2',PREV/'stage01/commands.json',
          ROOT/'outputs/v236_validation/build01/command.json',GAS/'build02/dependencies.json',ROOT/'outputs/v167_validation/build01/source_tree.json']
        return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,256,locks,BASE)
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for rows in read(PREV/'build01/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    locks += [STAGE/'commands.json',STAGE/'dependencies.json',PREV/'build01/members.json']
    if stage=='build':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build'],180,256,list(dict.fromkeys(locks)),BASE,'build')
    if stage in ('off','observe'):
        assert read(BASE/'process_build/report.json')['passed'];case=BASE/(stage+'01');locks += [BUILD/'branch.exe',case/'case.cfg',case/'mesh.su2']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='query')
        return guard.supervise(stage,[str(BUILD/'branch.exe'),'case.cfg','.',stage],60-used,256,list(dict.fromkeys(locks)),case,'query')
    raise ValueError('TASK_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v245 fixed control closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
