"""模块职责：新合格前体网格的一次成熟Euler发展，复用数值对象并绑定实际源。"""
from pathlib import Path
import difflib,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_BUILD
BASE=ROOT/'outputs/v242_validation';DOC=ROOT/'docs/validation/v242'
PREP=BASE/'flow_prepare01';BUILD=BASE/'flow_build01';CASE=BASE/'run01'
OLD=ROOT/'outputs/v240_validation';NATIVE=ROOT/'native/su2PhaseAwareAir'
CLOSED=True


def prepare():
    """模块一：据实际独立网格审核生成唯一薄驱动和原面映射，数值组合原样复用。"""
    audit=read(BASE/'mesh_audit01/report.json');assert audit['passed']
    assert sha(BASE/'tetgen01/mesh.su2')==audit['meshSHA256']
    for folder in (PREP,BUILD,CASE):folder.mkdir(exist_ok=False)
    header=NATIVE/'forebodyEulerV240.hpp';source=OLD/'prepare02/forebodyEuler.cpp'
    previous=read(OLD/'process_map_build/report.json')['inputSHA256']
    assert sha(header)==previous[str(header)] and sha(source)==previous[str(source)]
    old_header=header.read_text(encoding='utf8');new_header=old_header.replace('namespace region240','namespace region242').replace('checked!=807','checked!=%d'%audit['outletNodes'])
    assert new_header!=old_header and 'checked!=783' in new_header
    (PREP/'forebodyEulerV242.hpp').write_text(new_header,encoding='utf8')
    original=source.read_text(encoding='utf8');text=original.replace('forebodyEulerV240.hpp','forebodyEulerV242.hpp').replace('region240::','region242::').replace('g->GetnPoint()!=6111','g->GetnPoint()!=%d'%audit['nodes'])
    assert 'g->GetnPoint()!=5354' in text and 'if(argc!=3)' in text
    (PREP/'forebodyEuler.cpp').write_text(text,encoding='utf8')
    for name,before,after in [('driver',original,text),('header',old_header,new_header)]:
        (PREP/(name+'.patch')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='v240/'+name,tofile='v242/'+name)),encoding='utf8')
    command=read(OLD/'build01/link_command.json')
    replace={str(source):str(PREP/'forebodyEuler.cpp'),str(OLD/'build01/driver.d'):str(BUILD/'driver.d'),str(OLD/'build01/forebody.exe'):str(BUILD/'forebody.exe')}
    command=[replace.get(s,s) for s in command];save(BUILD/'link_command.json',command)
    # -M只提取真实包含依赖；已成熟对象/归档依照原成功回执逐SHA绑定。
    filtered=[];i=0
    while i<len(command):
        if command[i] in ('-MF','-o'):i+=2;continue
        if command[i] in ('-MD','-Wl,--start-group','-Wl,--end-group') or command[i].endswith(('.a','.obj')):i+=1;continue
        filtered.append(command[i]);i+=1
    dep=PREP/'driver.d';dcmd=filtered+['-M','-MT','dep','-MF',str(dep)]
    save(PREP/'dependency_command.json',dcmd);subprocess.check_call(dcmd,cwd=str(OLD_BUILD),timeout=20)
    known=read(OLD/'fix01/dependencies.json');new=[PREP/'forebodyEuler.cpp',PREP/'forebodyEulerV242.hpp'];bound={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
        if p not in new:assert known[str(p)]==h,str(p)
        bound[str(p)]=h
    for s in command:
        if s.endswith(('.a','.obj')):assert sha(s)==known[s];bound[s]=sha(s)
    shutil.copyfile(OLD/'run01/case.cfg',CASE/'case.cfg');shutil.copyfile(BASE/'tetgen01/mesh.su2',CASE/'mesh.su2')
    z=np.load(str(BASE/'tetgen01/mesh.npz'))
    np.savetxt(str(CASE/'original_faces.csv'),np.column_stack([z['originalFaceIDs'],z['wall'],z['originalArea']]),
      delimiter=',',fmt=['%d','%d','%d','%d','%.17g'],comments='',header='originalFaceID,node0,node1,node2,originalAreaM2')
    save(PREP/'dependencies.json',bound)
    save(PREP/'report.json',dict(passed=True,meshAuditSHA256=sha(BASE/'mesh_audit01/report.json'),
      nodes=audit['nodes'],outletNodes=audit['outletNodes'],boundDependencies=len(bound),matureObjectsCompiled=0,
      configBitEqualToV240=True,correctedGlobalToLocalMappingReused=True,newUniformInitialState=True,
      numericalFormulaChanges=False,physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v242 flow fixed execution closed; use NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=CASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p)==CASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/prepareSU2PhaseAirV220.py',DOC/'PLAN.md',DOC/'CFD_PLAN.md',
      BASE/'mesh_audit01/report.json',BASE/'tetgen01/mesh.su2',BASE/'tetgen01/mesh.npz']
    if stage=='prepare':
        locks += [NATIVE/'forebodyEulerV240.hpp',OLD/'prepare02/forebodyEuler.cpp',OLD/'process_map_build/report.json',OLD/'fix01/dependencies.json',OLD/'build01/link_command.json',OLD/'run01/case.cfg']
        command=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];seconds=45;kind='offline';cwd=BASE
    else:
        for p,h in read(PREP/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
        locks += [PREP/'dependencies.json',PREP/'report.json',BUILD/'link_command.json']
        if stage=='build':command=read(BUILD/'link_command.json');seconds=60;kind='build-flow';cwd=OLD_BUILD
        elif stage=='run':
            assert read(BASE/'process_flow_build/report.json')['passed']
            locks += [BUILD/'forebody.exe']+[p for p in CASE.iterdir() if p.is_file()]
            save(BASE/'flow_admission.json',dict(unique=True,updatesLimit=400,secondsLimit=600,meshSHA256=sha(CASE/'mesh.su2'),
              executableSHA256=sha(BUILD/'forebody.exe'),sourceReportSHA256=sha(PREP/'report.json'),physicalQualification=False))
            locks.append(BASE/'flow_admission.json');command=[str(BUILD/'forebody.exe'),'case.cfg','.'];seconds=600;kind='cfd';cwd=CASE
        else:raise ValueError('NO_ADMISSION')
    snap=BASE/('inputs_flow_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ('.py','.md','.cpp','.hpp'):(snap/(sha(p)+p.suffix)).write_bytes(p.read_bytes())
    return guard.supervise('flow_'+stage,command,seconds,256,locks,cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v242 flow fixed execution closed; use NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
