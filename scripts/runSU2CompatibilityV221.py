"""模块职责：一次完整边界控制的准备/构建/查询监督，禁止流场推进。"""
from pathlib import Path
import json,shutil,subprocess,sys
import runSU2PhaseAirV220 as supervision
from prepareSU2PhaseAirV220 import ROOT,SOURCE,OLD_SOURCE,OLD_BUILD,sha,read,save
BASE=ROOT/'outputs/v221_validation';BUILD=BASE/'build'
supervision.BASE=BASE
EXECUTION_CLOSED=True

def compiler_command():
    command=read(ROOT/'outputs/v220_validation/process_recover_helper/report.json')['argv']
    old=ROOT/'outputs/v220_validation/build02'
    return [str(BUILD/'compatibility.d') if x==str(old/'probe02.d') else
            str(BUILD/'compatibility.exe') if x==str(old/'phase_probe02.exe') else
            str(ROOT/'native/su2PhaseAwareAir/compatibility.cpp') if x==str(ROOT/'native/su2PhaseAwareAir/probe.cpp') else x for x in command]

def prepare():
    BUILD.mkdir(exist_ok=False);out=BASE/'control01';out.mkdir(exist_ok=False)
    mesh=ROOT/'outputs/v166_validation/small01/mesh.su2';shutil.copyfile(mesh,out/'mesh.su2')
    old=(ROOT/'outputs/v166_validation/small01/case.cfg').read_text(encoding='utf8')
    cfg={}
    for line in old.splitlines():
        if '=' in line:
            k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    for k in ('SST_OPTIONS','FREESTREAM_TURBULENCEINTENSITY','FREESTREAM_TURB2LAMVISCRATIO',
              'MU_REF','MU_T_REF','SUTHERLAND_CONSTANT'):cfg.pop(k,None)
    cfg.update(SA_OPTIONS='(NEGATIVE, WITHFT2)',KIND_TURB_MODEL='SA',RESTART_SOL='NO',FLUID_MODEL='RAH_PHASE_AIR',
        VISCOSITY_MODEL='RAH_PHASE_LINEAR',CONDUCTIVITY_MODEL='RAH_PHASE_LINEAR',FREESTREAM_NU_FACTOR='3',
        FREESTREAM_PRESSURE='1197.0279007383137',FREESTREAM_TEMPERATURE='226.50908361133003',
        MARKER_RIEMANN='(farfield, STATIC_SUPERSONIC_INFLOW_PT, 1197, 350, 2.001, 0.401, 0.299)',
        ITER='1',CONV_FIELD='(RMS_DENSITY, RMS_ENERGY, RMS_NU_TILDE)',OUTPUT_FILES='(RESTART_ASCII)',
        OUTPUT_WRT_FREQ='1',VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)',
        WRT_RESTART_OVERWRITE='YES',SOLUTION_FILENAME='restart_flow',RESTART_FILENAME='restart_flow')
    (out/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    cmd=compiler_command();dep=BUILD/'inputs.d';filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(OLD_BUILD),timeout=20)
    names=dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split();bindings={}
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');known=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n in names:
        p=Path(n);p=p if p.is_absolute() else OLD_BUILD/p;p=p.resolve();h=sha(p);bindings[str(p)]=h
        if str(p) in known:assert h==known[str(p)],n
        else:
            try:relative=p.relative_to(SOURCE)
            except ValueError:continue
            assert h==tree[str(OLD_SOURCE/relative)],n
    save(BUILD/'inputs.json',bindings);save(BUILD/'command.json',cmd)
    save(BUILD/'prepare_report.json',dict(passed=True,dependencies=len(bindings),meshSHA256=sha(mesh),configSHA256=sha(out/'case.cfg')))

def main(stage):
    BASE.mkdir(exist_ok=True);locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2PhaseAirV220.py',
        ROOT/'docs/validation/v221/PLAN.md',ROOT/'native/su2PhaseAwareAir/compatibility.cpp']
    if stage=='prepare':
        return supervision.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],'offline',30,64,locks,BASE)
    inputs=read(BUILD/'inputs02.json');assert all(sha(Path(n))==h for n,h in inputs.items())
    locks += [Path(n) for n in inputs]+[BUILD/'inputs02.json']
    members=read(ROOT/'outputs/v220_validation/build02/members.json')
    for group in members.values():
        for row in group:assert sha(Path(row['path']))==row['sha256'];locks.append(Path(row['path']))
    if stage=='build2':
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        return supervision.supervise(stage,compiler_command(),'build',min(120,240-used),256,locks,BASE)
    if stage=='control':
        assert read(BASE/'process_build2/report.json')['passed']
        out=BASE/'control01';locks += [BUILD/'compatibility.exe',out/'case.cfg',out/'mesh.su2']
        return supervision.supervise(stage,[str(BUILD/'compatibility.exe'),str(out/'case.cfg'),str(out),'control'],'query',30,64,locks,out)
    if stage=='restart':
        from su2PhaseAirContract import bind,read_restart
        assert read(BASE/'audit01/report.json')['passed']
        origin=BASE/'control01';out=BASE/'restart01';out.mkdir(exist_ok=False)
        metadata=read(origin/'restart_flow.rah.json')
        expected=bind(origin/'restart_flow.csv',origin/'mesh.su2',origin/'case.cfg',
            ROOT/'outputs/v220_validation/build02/SU2_CFD.exe',BUILD/'compatibility.exe','manufactured-control')
        read_restart(origin/'restart_flow.csv',metadata,expected,origin/'mesh.su2')
        for n in ('mesh.su2','restart_flow.csv','restart_flow.rah.json'):shutil.copyfile(origin/n,out/n)
        cfg=(origin/'case.cfg').read_text(encoding='utf8');assert cfg.count('RESTART_SOL= NO')==1
        (out/'case.cfg').write_text(cfg.replace('RESTART_SOL= NO','RESTART_SOL= YES'),encoding='utf8')
        locks += [BUILD/'compatibility.exe',ROOT/'scripts/su2PhaseAirContract.py',out/'case.cfg',out/'mesh.su2',
                  out/'restart_flow.csv',out/'restart_flow.rah.json',origin/'case.cfg']
        return supervision.supervise(stage,[str(BUILD/'compatibility.exe'),str(out/'case.cfg'),str(out),'restart'],'query',30,64,locks,out)
    raise ValueError('Stage not admitted')
if __name__=='__main__':
    if EXECUTION_CLOSED:raise SystemExit('v221 control closed; next task is the separately budgeted real rocket candidate')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
