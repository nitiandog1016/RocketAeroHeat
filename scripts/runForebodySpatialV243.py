"""模块职责：唯一细网格真实外流对照，按原父面聚合观察，复用成熟数值对象。"""
from pathlib import Path
import difflib,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_BUILD
BASE=ROOT/'outputs/v243_validation';DOC=ROOT/'docs/validation/v243';OLD=ROOT/'outputs/v242_validation'
PREP=BASE/'flow_prepare01';BUILD=BASE/'flow_build01';CASE=BASE/'run01';CLOSED=True


def once(s,a,b):
    assert s.count(a)==1,a[:80];return s.replace(a,b)


def prepare():
    """模块一：父子观察适配不影响真实推进算子，源和所有成熟归档逐SHA绑定。"""
    audit=read(BASE/'reconnected_audit01/report.json');assert audit['passed']
    for p in (PREP,BUILD,CASE):p.mkdir(exist_ok=False)
    header=OLD/'flow_prepare01/forebodyEulerV242.hpp';source=OLD/'flow_prepare01/forebodyEuler.cpp'
    known=read(OLD/'flow_prepare01/dependencies.json');assert sha(header)==known[str(header)] and sha(source)==known[str(source)]
    original=header.read_text(encoding='utf8');h=original.replace('namespace region242','namespace region243').replace('checked!=783','checked!=3022')
    h=h.replace('#include <cstring>','#include <cstring>\n#include <map>')
    h=once(h,'double area;};','double area;size_t parent;};')
    h=once(h,'std::vector<Face> faces;','std::vector<Face> faces;\nstd::vector<double> parent_area;\nstd::map<unsigned long,size_t> parent_index;')
    h=once(h,'q.a=local[q.a];q.b=local[q.b];q.c=local[q.c];faces.push_back(q);',
      '''q.a=local[q.a];q.b=local[q.b];q.c=local[q.c];
    // 原面压力必须先按子面积聚合，不能把细化后的子面当新物理原面。
    auto inserted=parent_index.emplace(q.id,parent_index.size());q.parent=inserted.first->second;
    if(inserted.second)parent_area.push_back(0.0);
    parent_area[q.parent]+=q.area;faces.push_back(q);''')
    h=once(h,'faces.size()!=4434||cfg->GetRestart()','faces.size()!=17736||parent_area.size()!=4434||cfg->GetRestart()')
    h=once(h,'std::vector<double> result;for(auto q:faces)result.push_back((p[q.a]+p[q.b]+p[q.c])/3);',
      'std::vector<double> result(parent_area.size(),0.0);\n  for(auto q:faces)result[q.parent]+=q.area*(p[q.a]+p[q.b]+p[q.c])/3;\n  for(size_t i=0;i<result.size();++i)result[i]/=parent_area[i];')
    lo=h.index('inline void readback(');hi=h.index('inline void record(',lo);h=h[:lo]+h[hi:]
    h=once(h,'for(size_t i=0;i<faces.size();++i){auto q=faces[i];double value=surface[i];','for(size_t i=0;i<parent_area.size();++i){double value=surface[i];')
    h=h.replace('q.area*std::abs(value-previous_pressure[i])','parent_area[i]*std::abs(value-previous_pressure[i])').replace('q.area*std::abs(previous_pressure[i])','parent_area[i]*std::abs(previous_pressure[i])')
    target_header=PREP/'forebodyEulerV243.hpp';target_header.write_text(h,encoding='utf8')
    cpp=source.read_text(encoding='utf8');new=cpp.replace('forebodyEulerV242.hpp','forebodyEulerV243.hpp').replace('region242::','region243::').replace('g->GetnPoint()!=5354','g->GetnPoint()!=34956')
    target=PREP/'forebodyEuler.cpp';target.write_text(new,encoding='utf8')
    for name,before,after in [('driver',cpp,new),('header',original,h)]:
        (PREP/(name+'.patch')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='v242/'+name,tofile='v243/'+name)),encoding='utf8')
    command=read(OLD/'flow_build01/link_command.json');replacements={str(source):str(target),str(OLD/'flow_build01/forebody.exe'):str(BUILD/'forebody.exe'),str(OLD/'flow_build01/driver.d'):str(BUILD/'driver.d')}
    command=[replacements.get(s,s) for s in command];save(BUILD/'link_command.json',command)
    filtered=[];i=0
    while i<len(command):
        if command[i] in ('-MF','-o'):i+=2;continue
        if command[i] in ('-MD','-Wl,--start-group','-Wl,--end-group') or command[i].endswith(('.a','.obj')):i+=1;continue
        filtered.append(command[i]);i+=1
    dep=PREP/'driver.d';dcmd=filtered+['-M','-MT','dep','-MF',str(dep)];save(PREP/'dependency_command.json',dcmd)
    subprocess.check_call(dcmd,cwd=str(OLD_BUILD),timeout=20);bound={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();v=sha(p)
        if p not in (target,target_header):assert known[str(p)]==v,str(p)
        bound[str(p)]=v
    for s in command:
        if s.endswith(('.a','.obj')):assert sha(s)==known[s];bound[s]=sha(s)
    cfg=(OLD/'run01/case.cfg').read_text(encoding='utf8');cfg=once(cfg,'ITER= 400','ITER= 800');cfg=once(cfg,'OUTPUT_WRT_FREQ= 50','OUTPUT_WRT_FREQ= 800')
    (CASE/'case.cfg').write_text(cfg,encoding='utf8');shutil.copyfile(BASE/'reconnected_audit01/mesh.su2',CASE/'mesh.su2')
    z=np.load(str(BASE/'reconnected01/mesh.npz'));order=np.argsort(z['originalFaceIDs'],kind='stable')
    np.savetxt(str(CASE/'original_faces.csv'),np.column_stack([z['originalFaceIDs'][order],z['wall'][order],z['childArea'][order]]),
      delimiter=',',fmt=['%d','%d','%d','%d','%.17g'],comments='',header='originalFaceID,node0,node1,node2,childAreaM2')
    save(PREP/'dependencies.json',bound)
    save(PREP/'report.json',dict(passed=True,matureObjectsCompiled=0,thinDriverOnly=True,sourceParentFaces=4434,computationalWallFaces=17736,
      originalFaceAreaAggregation=True,numericalFormulasChanged=False,uniformInitialState=True,meshSHA256=sha(CASE/'mesh.su2'),physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v243 fixed spatial flow closed; use NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=CASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p)==CASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/prepareSU2PhaseAirV220.py',DOC/'PLAN.md',DOC/'CFD_PLAN.md',DOC/'METHOD_DECISION.md',
      BASE/'reconnected_audit01/report.json',BASE/'reconnected_audit01/mesh.su2',BASE/'reconnected01/mesh.npz']
    if stage=='prepare':
        locks += [OLD/'flow_prepare01'/n for n in ('forebodyEuler.cpp','forebodyEulerV242.hpp','dependencies.json')]+[OLD/'flow_build01/link_command.json',OLD/'run01/case.cfg']
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];sec=45;kind='offline';cwd=BASE
    else:
        for p,h in read(PREP/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
        locks += [PREP/'dependencies.json',PREP/'report.json',BUILD/'link_command.json']
        if stage=='build':cmd=read(BUILD/'link_command.json');sec=60;kind='build';cwd=OLD_BUILD
        elif stage=='run':
            assert read(BASE/'process_flow_build/report.json')['passed'];locks += [BUILD/'forebody.exe']+[p for p in CASE.iterdir() if p.is_file()]
            save(BASE/'flow_admission.json',dict(unique=True,updatesLimit=800,secondsLimit=600,outputMiB=384,meshSHA256=sha(CASE/'mesh.su2'),executableSHA256=sha(BUILD/'forebody.exe'),physicalQualification=False))
            locks.append(BASE/'flow_admission.json');cmd=[str(BUILD/'forebody.exe'),'case.cfg','.'];sec=600;kind='cfd';cwd=CASE
        else:raise ValueError('NO_ADMISSION')
    snap=BASE/('inputs_flow_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ('.py','.md','.cpp','.hpp'):(snap/(sha(p)+p.suffix)).write_bytes(p.read_bytes())
    return guard.supervise('flow_'+stage,cmd,sec,384,locks,cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 fixed spatial flow closed; use NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
