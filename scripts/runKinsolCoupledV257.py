"""模块职责：固定成熟库/薄适配构建、完整耦合控制及准确预算；不启动旧控制。"""
import sys,os,json,hashlib,tarfile,shutil,subprocess,re,urllib.request
from pathlib import Path,PurePosixPath
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v257_validation';DOC=ROOT/'docs/validation/v257';RT=ROOT/'tools/runtime/sundials-7.4.0-8e17876'
SRC=RT/'source';LIBBUILD=RT/'build';TC=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'
PREV=ROOT/'outputs/v256_validation';WD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
BUILD=BASE/'build01';CASE=BASE/'rocket01';STAGE=BASE/'stage02';COMMIT='8e17876d3b4d682b4098684b07a85b005a122f81';CLOSED=True
HEADER=ROOT/'native/su2PhaseAwareAir/FullResidualProviderV257.hpp';DRIVER=ROOT/'native/su2PhaseAwareAir/kinsolCoupledV257.cpp'


def prepare():
    assert not read(BASE/'process_prepare/report.json')['passed'];STAGE.mkdir(exist_ok=False)
    url='https://codeload.github.com/LLNL/sundials/tar.gz/'+COMMIT;archive=RT/'fixed-source.tar.gz';partial=RT/'fixed-source.partial'
    save(STAGE/'download_request.json',dict(url=url,commit=COMMIT,destination=str(archive),limitBytes=128*1024**2))
    assert archive.is_file() and SRC.is_dir() and not partial.exists();total=archive.stat().st_size
    files={};selected={'CMakeLists.txt','cmake','src','include','examples','LICENSE','NOTICE','README.md','CITATIONS.md'}
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar.getmembers():
            parts=PurePosixPath(member.name).parts
            assert parts and parts[0]=='sundials-'+COMMIT and '..' not in parts
            if len(parts)<2 or parts[1] not in selected or not member.isfile():continue
            if parts[1]=='examples' and not (member.name.endswith('.txt') or member.name.endswith('.cmake')):continue
            relative=Path(*parts[1:]);dest=(SRC/relative).resolve();assert SRC.resolve() in dest.parents
            with tar.extractfile(member) as inp:expected=hashlib.sha256(inp.read()).hexdigest()
            assert dest.is_file() and sha(dest)==expected,relative
            files[relative.as_posix()]=expected
    references={'src/kinsol/kinsol.c':PREV/'research/sundials/src/kinsol/kinsol.c',
        'src/kinsol/kinsol_ls.c':PREV/'research/sundials/src/kinsol/kinsol_ls.c',
        'LICENSE':PREV/'research/sundials/LICENSE','NOTICE':BASE/'research/NOTICE_exact'}
    for name in ['src/sunlinsol/spgmr/sunlinsol_spgmr.c','include/kinsol/kinsol.h','include/kinsol/kinsol_ls.h']:references[name]=BASE/'research'/name
    for name,p in references.items():assert sha(SRC/name)==sha(p),name
    assert (PREV/'research/sundials/NOTICE').read_bytes()==(SRC/'NOTICE').read_bytes()+b'\n'
    N=1974844;V=N*8;graph=read(ROOT/'outputs/v251_validation/stage01/report.json');assert graph['matrixBlocks']==3079394
    components=dict(measuredNativeReference=read(PREV/'process_control2/report.json')['peakObservedGroupRSSBytes'],
        ILUFlowSA=graph['matrixBlocks']*26*8,SPGMR23Vectors=23*V,KIN8VectorAllowance=8*V,
        caller4Vectors=4*V,adapter8VectorAllowance=8*V,native4BlockVectors=335726*12*8,reserve=256*1024**2)
    forecast=sum(components.values());assert forecast<4*1024**3
    configure=[str(TC/'cmake.exe'),'-S',str(SRC),'-B',str(LIBBUILD),'-G','Ninja','-DCMAKE_C_COMPILER='+str(TC/'gcc.exe'),
        '-DCMAKE_MAKE_PROGRAM='+str(TC/'ninja.exe'),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_SHARED_LIBS=OFF','-DBUILD_STATIC_LIBS=ON',
        '-DBUILD_ARKODE=OFF','-DBUILD_CVODE=OFF','-DBUILD_CVODES=OFF','-DBUILD_IDA=OFF','-DBUILD_IDAS=OFF','-DBUILD_KINSOL=ON',
        '-DBUILD_TESTING=OFF','-DBUILD_BENCHMARKS=OFF','-DEXAMPLES_ENABLE_C=OFF','-DEXAMPLES_ENABLE_CXX=OFF','-DEXAMPLES_INSTALL=OFF',
        '-DBUILD_FORTRAN_MODULE_INTERFACE=OFF','-DENABLE_OPENMP=OFF','-DENABLE_MPI=OFF','-DSUNDIALS_INDEX_SIZE=64',
        '-DSUNDIALS_PRECISION=double','-DSUNDIALS_LOGGING_LEVEL=3']
    save(STAGE/'library_commands.json',dict(configure=configure,build=[str(TC/'cmake.exe'),'--build',str(LIBBUILD),'--parallel','1']))
    save(STAGE/'library_source.json',dict(commit=COMMIT,archiveURL=url,archiveSHA256=sha(archive),archiveBytes=total,files=files,githubConnectorReferenceFiles={n:sha(p) for n,p in references.items()},license='BSD-3-Clause'))
    cfg=(PREV/'rocket02/case.cfg').read_text(encoding='utf8')
    # 两版本case目录具有相同深度，原大输入相对路径不变；预条件仅ILU0，不改变空间方法。
    cfg,n=re.subn(r'^LINEAR_SOLVER_PREC=.*$','LINEAR_SOLVER_PREC= ILU',cfg,count=1,flags=re.M);assert n==1
    cfg+='\nLINEAR_SOLVER_ILU_FILL_IN= 0\nLINEAR_SOLVER_ILU_LEVEL_SCHEDULING= NO\n'
    (CASE/'case.cfg').write_text(cfg,encoding='utf8')
    save(STAGE/'report.json',dict(passed=True,unknowns=N,vectorBytes=V,resourceComponentsBytes=components,forecastBytes=forecast,
        forecastNotMeasuredPeak=True,memoryLimitBytes=4*1024**3,KrylovDimension=20,maxRestarts=0,eta=1e-5,
        nonlinearAttemptLimit=3,physicalQualification=False,newQualifiedAreaM2=0))


def library():
    commands=read(STAGE/'library_commands.json');env=dict(os.environ);env['PATH']=str(TC)+os.pathsep+env.get('PATH','')
    subprocess.check_call(commands['configure'],env=env);subprocess.check_call(commands['build'],env=env)
    wanted=['kinsol','nvecserial','sunlinsolspgmr','core'];libs={}
    for word in wanted:
        choices=list(LIBBUILD.rglob('libsundials_'+word+'.a'));assert len(choices)==1,(word,choices)
        libs[word]=dict(path=str(choices[0]),sha256=sha(choices[0]))
    size=sum(p.stat().st_size for p in RT.rglob('*') if p.is_file());assert size<256*1024**2
    save(STAGE/'library_build.json',dict(passed=True,commit=COMMIT,staticLibraries=libs,runtimeBytes=size,systemInstalled=False))


def bridge():
    prior=read(PREV/'stage02/command.json');command=[str(DRIVER) if p==str(ROOT/'native/su2PhaseAwareAir/fullDomainResidualV256.cpp') else p.replace(str(PREV/'build01'),str(BUILD)) for p in prior]
    command[1:1]=['-I'+str(SRC/'include'),'-I'+str(LIBBUILD/'include')]
    libs=read(STAGE/'library_build.json')['staticLibraries']
    at=command.index('-Wl,--end-group');command[at:at]=[libs[k]['path'] for k in ['kinsol','sunlinsolspgmr','nvecserial','core']]
    save(STAGE/'bridge_command.json',command)
    members=read(PREV/'build01/members.json');ar=str(TC/'ar.exe')
    for name,rows in members.items():
        for row in rows:assert sha(row['path'])==row['sha256']
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in rows],cwd=WD)
    save(BUILD/'members.json',members);subprocess.check_call(command,cwd=WD)
    generated={str(p):sha(p) for p in (LIBBUILD/'include').rglob('*') if p.is_file()}
    save(BUILD/'generated_headers.json',generated)
    save(BUILD/'report.json',dict(passed=True,physicsObjectsRecompiled=0,engineSHA256=sha(BUILD/'preflight.exe'),libraryCommit=COMMIT,
        generatedHeaderSHA256=generated,dependencyFileSHA256=sha(BUILD/'preflight.d')))


def main(stage):
    if CLOSED:raise SystemExit('v257 fixed KINSOL coupled tasks closed; follow NEXT_DECISION')
    assert stage in ['prepare2','library','bridge','control'],'首准备已封存，仅具名恢复。'
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        PREV/'process_control2/report.json',PREV/'rocket02/case.cfg',ROOT/'outputs/v251_validation/stage01/report.json']
    if stage=='prepare2':
        locks += [RT/'fixed-source.tar.gz',BASE/'research/NOTICE_exact']
        argv=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];limit=120-read(BASE/'process_prepare/report.json')['elapsedSeconds'];kind='offline';cwd=ROOT
    else:
        locks += [STAGE/'report.json',STAGE/'library_source.json',STAGE/'library_commands.json',CASE/'case.cfg']
        locks += [SRC/n for n in read(STAGE/'library_source.json')['files']]
        if stage=='library':argv=[sys.executable,str(Path(__file__).resolve()),'task_library'];limit=180;kind='build';cwd=ROOT
        else:
            locks += [HEADER,DRIVER,STAGE/'library_build.json']
            locks += [p for p in (LIBBUILD/'include').rglob('*') if p.is_file()]
            locks += [Path(r['path']) for r in read(STAGE/'library_build.json')['staticLibraries'].values()]
            locks += [ROOT/'native/su2PhaseAwareAir/singleDomainRunV235.cpp',PREV/'build01/members.json']
            for rows in read(PREV/'build01/members.json').values():locks += [Path(r['path']) for r in rows]
            # 复用实际成熟对象的原头依赖；新本轮头和SUNDIALS头另绑定。
            locks += [Path(p) for p in read(PREV/'stage02/dependencies.json')]
            if stage=='bridge':argv=[sys.executable,str(Path(__file__).resolve()),'task_bridge'];limit=120;kind='build';cwd=ROOT
            elif stage=='control':
                assert read(BUILD/'report.json')['passed'];locks += [BUILD/'report.json',BUILD/'preflight.d',BUILD/'generated_headers.json',BUILD/'preflight.exe',ROOT/'outputs/v250_validation/rocket01/mesh.su2',ROOT/'outputs/v250_validation/rocket01/input_restart.csv',PREV/'rocket02/initial_state.bin',PREV/'rocket02/base_full_residual.bin']
                # 构建后按真实依赖清单补齐新头；同时核对构建时绑定，不能仅锁住已变化的文件。
                previous=read(BASE/'process_bridge/report.json');assert previous['passed']
                for p,h in previous['inputSHA256'].items():
                    if Path(p)!=Path(__file__).resolve():assert sha(p)==h,p
                dependencies={}
                for name in (BUILD/'preflight.d').read_text(encoding='utf8').replace('\\\n',' ').split(': ',1)[1].split():
                    p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();dependencies[str(p)]=sha(p);locks.append(p)
                save(STAGE/'bridge_dependencies.json',dependencies);locks.append(STAGE/'bridge_dependencies.json')
                for name,h in read(STAGE/'library_source.json')['files'].items():assert sha(SRC/name)==h,name
                for row in read(STAGE/'library_build.json')['staticLibraries'].values():assert sha(row['path'])==row['sha256']
                assert sha(BUILD/'preflight.exe')==read(BUILD/'report.json')['engineSHA256']
                argv=[str(BUILD/'preflight.exe'),'case.cfg','.'];limit=900;kind='control';cwd=CASE
            else:raise ValueError('STAGE_NOT_ADMITTED')
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ['.py','.cpp','.hpp','.md','.json','.cfg','.h','.c','.inl'] and p.stat().st_size<2*1024**2:
            # 工程/计划精确快照；固定成熟源已独立manifest绑定且不会被修改。
            if p==Path(__file__).resolve() or p in [DOC/'PLAN.md',HEADER,DRIVER] or p.parent==STAGE:
                target=snap/(sha(p)+p.suffix)
                if not target.exists():shutil.copyfile(p,target)
    return supervise(guard,BASE,stage,argv,limit,1024,list(dict.fromkeys(locks)),cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v257 fixed KINSOL coupled tasks closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_library':library()
    elif sys.argv[1]=='task_bridge':bridge()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
