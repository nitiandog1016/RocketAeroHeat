"""模块职责：一次绑定成熟ILU存储与精度语义，复用v250对象运行唯一固定真实分支。"""
from pathlib import Path
import os,re,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v251_validation';DOC=ROOT/'docs/validation/v251';CASE=BASE/'rocket01'
STAGE=BASE/'stage01';BUILD=BASE/'build01';CLOSED=True
PREV=ROOT/'outputs/v250_validation';GAS=ROOT/'outputs/v220_validation'
WD=ROOT/'tools/runtime/su2-thermal-build/build-v117';NATIVE=ROOT/'native/su2PhaseAwareAir'
DRIVER=NATIVE/'nativeILUV251.cpp';HOOK=NATIVE/'NewtonHooksV250.hpp'
REFERENCE=ROOT/'outputs/v248_validation/rocket01/state_000200.bin'
MEMBERS=PREV/'build01/members.json';FACES=PREV/'rocket01/original_faces.txt'
SOURCES={
 'factory':GAS/'source/Common/include/linear_algebra/CPreconditioner.hpp',
 'matrix':ROOT/'tools/runtime/su2-thermal-build/source/Common/src/linear_algebra/CSysMatrix.cpp',
 'krylov':ROOT/'tools/runtime/su2-thermal-build/source/Common/src/linear_algebra/CSysSolve.cpp',
 'config':GAS/'source/Common/src/CConfig.cpp',
 'newton':PREV/'stage01/CNewtonIntegration_original.cpp',
 'inner':GAS/'source/SU2_CFD/include/integration/CNewtonIntegration.hpp',
 'scalar':GAS/'source/SU2_CFD/include/solvers/CScalarSolver.inl'}

def setup_guard():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes
    guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    return guard

def prepare():
    """模块一：真实原单元边图、成熟源对象及额外存储绑定；不运行流动或小控制。"""
    for p in (CASE,STAGE,BUILD):
        if p.exists():assert not any(p.iterdir()),str(p)
        else:p.mkdir()
    source_tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    known=read(PREV/'stage01/dependencies.json')
    source_rows=[]
    ranges={'factory':[(207,251),(429,440)],'matrix':[(184,259),(336,355),(404,411),(1064,1196)],
      'krylov':[(425,455),(467,501),(526,624),(1460,1472),(1501,1540)],
      'config':[(2000,2041)],'newton':[(259,315)],'inner':[(163,189)],'scalar':[(574,599)]}
    for key,p in SOURCES.items():
        expected=known.get(str(p),source_tree.get(str(p)))
        if key=='scalar':expected=source_tree.get(str(ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/include/solvers/CScalarSolver.inl'))
        if key=='config':expected='241eb647b07a136eb87f5528960db0737ae933ef5e45579be97218d7d692a48e'
        if key=='newton':expected=read(PREV/'stage01/report.json')['originalNewtonSHA256']
        assert expected and sha(p)==expected,(key,sha(p),expected)
        lines=p.read_text(encoding='utf8').splitlines()
        source_rows.append(dict(role=key,path=str(p),sha256=expected,excerpts=[dict(first=a,last=b,text='\n'.join(lines[a-1:b])) for a,b in ranges[key]]))
    # 四面体六条边、棱柱九条边；与原生FVM邻接一致，不按单元全顶点两两连接。
    mesh=PREV/'rocket01/mesh.su2';parts={10:[],13:[]}
    with mesh.open(encoding='utf8') as f:
        for line in f:
            if line.startswith('NELEM='):
                count=int(line.split('=')[1]);break
        for _ in range(count):
            a=[int(x) for x in next(f).split()];kind=a[0];assert kind in parts
            parts[kind].append(a[1:5] if kind==10 else a[1:7])
    assert count==705490 and len(parts[10])==70574 and len(parts[13])==634916
    pairs={10:((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)),13:((0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5))}
    edges=[];N=335726
    for kind,rows in parts.items():
        a=np.array(rows,dtype=np.int64)
        for i,j in pairs[kind]:edges.append(np.minimum(a[:,i],a[:,j])*N+np.maximum(a[:,i],a[:,j]))
    keys=np.unique(np.concatenate(edges));edge_count=len(keys);blocks=N+2*edge_count
    # 所有本机成熟标量按8字节保守计算；实际驱动再次核对sizeof和真实图计数。
    extra_ilu=blocks*(25+1)*8;extra_krylov=2*(20-10)*N*5*8
    previous_peak=read(PREV/'process_rocket/report.json')['peakObservedGroupRSSBytes']
    forecast=previous_peak+extra_ilu+extra_krylov+256*1024**2
    assert forecast<4*1024**3,(forecast,edge_count)
    for name in ('mesh.su2','input_restart.csv'):os.link(str(PREV/'rocket01'/name),str(CASE/name))
    shutil.copyfile(FACES,CASE/'original_faces.txt')
    cfg=(PREV/'rocket01/case.cfg').read_text(encoding='utf8')
    for key,value in [('LINEAR_SOLVER_PREC','ILU'),('LINEAR_SOLVER_ITER','20')]:
        cfg,n=re.subn('^'+key+r'=.*$',key+'= '+value,cfg,count=1,flags=re.M);assert n==1
    cfg+='\nLINEAR_SOLVER_ILU_FILL_IN= 0\nLINEAR_SOLVER_ILU_LEVEL_SCHEDULING= NO\nLINEAR_SOLVER_INNER= NONE\n'
    (CASE/'case.cfg').write_text(cfg,encoding='utf8')
    schema=read(PREV/'rocket01/schema.json');schema.update(schema='RAH.NativeILU.v251',flowUpdateGate='native scaled report <= 1e-5; before CompleteImplicitIteration',
      SA='native BCGSTAB/ILU(0), 20/1e-5; observed after update, no SA pre-update gate')
    save(CASE/'schema.json',schema)
    command=[str(DRIVER) if x==str(NATIVE/'newtonDevelopmentV250.cpp') else x.replace(str(PREV/'build01'),str(BUILD)) for x in read(PREV/'stage01/commands.json')['helper']]
    save(STAGE/'command.json',command)
    args=[];i=0
    while i<len(command):
        x=command[i]
        if x in ('-o','-MF'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    dep=STAGE/'deps.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=str(WD),timeout=25)
    deps={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
        if p!=DRIVER:assert h==known[str(p)],str(p)
        deps[str(p)]=h
    save(STAGE/'dependencies.json',deps);shutil.copyfile(DRIVER,STAGE/DRIVER.name)
    for name in ('CNewtonIntegration_original.cpp','CNewtonIntegration_observed.cpp'):shutil.copyfile(PREV/'stage01'/name,STAGE/name)
    save(STAGE/'report.json',dict(passed=True,adoptedForOneBoundedBranch=True,productionQualified=False,sources=source_rows,
      meshSHA256=sha(mesh),points=N,edges=edge_count,matrixBlocks=blocks,extraILUFlowSABytes=extra_ilu,
      extraOuterKrylovBytes=extra_krylov,previousMeasuredPeakBytes=previous_peak,forecastWith256MiBReserve=forecast,
      memoryLimitBytes=4*1024**3,forecastIsNotPeakProof=True,outerLimit=20,reportedTolerance=1e-5,
      innerLimit=3,innerTolerance=.1,ILUFill=0,independentJacobianAccuracy=False,matureObjectsToRecompile=0))
    print(dict(passed=True,edges=edge_count,forecastGiB=forecast/1024**3))

def build():
    """模块二：薄归档逐位复用全部成熟对象（含v250观察对象），仅编译本次驱动。"""
    command=read(STAGE/'command.json');ar=str(Path(command[0]).parent/'ar.exe')
    members=read(MEMBERS)
    for name,items in members.items():
        for r in items:assert sha(r['path'])==r['sha256']
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in items],cwd=str(WD))
    save(BUILD/'members.json',members);subprocess.check_call(command,cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,matureObjectsRecompiled=0,observerObjectsRecompiled=0,engineSHA256=sha(BUILD/'preflight.exe')))

def main(mode):
    if CLOSED:raise SystemExit('v251 native ILU branch closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);g=setup_guard()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',DOC/'PLAN.md',DRIVER,HOOK,REFERENCE,MEMBERS,FACES,
      PREV/'stage01/dependencies.json',PREV/'stage01/commands.json',PREV/'stage01/report.json',PREV/'rocket01/case.cfg',PREV/'rocket01/schema.json',
      PREV/'rocket01/mesh.su2',PREV/'rocket01/input_restart.csv',PREV/'process_rocket/report.json',ROOT/'outputs/v167_validation/build01/source_tree.json']+list(SOURCES.values())
    if mode=='prepare':return supervise(g,BASE,'prepare2',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,1536,locks,BASE,'offline')
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for items in read(MEMBERS).values():
        for r in items:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    locks += [STAGE/'command.json',STAGE/'dependencies.json',STAGE/'report.json',CASE/'case.cfg',CASE/'mesh.su2',CASE/'input_restart.csv',CASE/'original_faces.txt',CASE/'schema.json']
    if mode=='build':return supervise(g,BASE,mode,[sys.executable,str(Path(__file__).resolve()),'task_build'],120,1536,list(dict.fromkeys(locks)),BASE,'build')
    if mode=='run':
        assert read(BASE/'process_build/report.json')['passed'] and read(BASE/'process_build/terminal_capacity.json')['passed']
        capacity=read(PREV/'stage01/report.json')['capacityComponentsBytes']
        capacity['refusedDirectionBinaryAllowance']=24+335726*17*8
        assert sum(capacity.values())<1536*1024**2
        save(BASE/'run_admission.json',dict(capacityComponentsBytes=capacity,projectedLogicalBytes=sum(capacity.values()),
          limitBytes=1536*1024**2,maximumUpdates=200,timeSeconds=1800,normalStopReserveSeconds=120,
          formalBenefitWindowNotGuaranteed=True,oneFixedCombination=True,automaticRetry=False))
        locks.append(BASE/'run_admission.json')
        locks += [BUILD/'preflight.exe',BUILD/'members.json']
        return supervise(g,BASE,'rocket',[str(BUILD/'preflight.exe'),'case.cfg','.',str(REFERENCE),'original_faces.txt'],1800,1536,list(dict.fromkeys(locks)),CASE,'cfd')
    raise ValueError('MODE_NOT_ADMITTED')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v251 native ILU branch closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
