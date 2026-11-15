"""模块职责：精确绑定成熟Euler对象，只增加真实出口阶段守卫，监督唯一前体区域发展。"""
from pathlib import Path
import difflib,json,shutil,subprocess,sys
import numpy as np
from scipy.spatial import ConvexHull
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_BUILD
from phaseEulerRegionalVolume import read_tetra_mesh
BASE=ROOT/'outputs/v240_validation';DOC=ROOT/'docs/validation/v240'
CASE=BASE/'run01';BUILD=BASE/'build01';PREP=BASE/'prepare01';NATIVE=ROOT/'native/su2PhaseAwareAir'
OLD=ROOT/'outputs/v229_validation';ORIGINAL=OLD/'source/SU2_CFD/src/solvers/CEulerSolver_observed.cpp'
CLOSED=False


def clipped_volume_audit():
    """模块一：独立逐裁切母单元凸包核对真实分割体积，不只检查总闭合。"""
    source=ROOT/'outputs/v232_validation/rocket01/mesh.su2';x,c,_=read_tetra_mesh(source)
    z=np.load(str(BASE/'region01/provenance.npz'));xx=z['xyz'];cc=z['cells'];parent=z['motherCell'];cut=float(z['cutX'])
    dx=xx[cc[:,1:]]-xx[cc[:,0,None]]
    volume=np.linalg.det(dx)/6;assert np.all(volume>0)
    actual=np.bincount(parent,weights=volume,minlength=len(c));maximum=0.
    for cid in z['partialMotherCells']:
        v=x[c[cid]];points=[p for p in v if p[0]>=cut]
        for i in range(4):
            for j in range(i):
                a,b=v[i],v[j]
                if (a[0]-cut)*(b[0]-cut)<0:points.append(((b[0]-cut)*a+(cut-a[0])*b)/(b[0]-a[0]))
        expected=ConvexHull(np.asarray(points)).volume
        error=abs(actual[cid]-expected)/expected;maximum=max(maximum,error)
    assert maximum<2e-10,maximum
    mapping=np.sum(x[z['sourceNode']]*z['sourceWeight'][:,:,None],axis=1)
    linear=float(np.max(abs(mapping-xx)));assert linear<5e-14
    longest=np.max(np.linalg.norm(xx[cc[:,[0,0,0,1,1,2]]]-xx[cc[:,[1,2,3,2,3,3]]],axis=2),axis=1)
    quality=6*np.sqrt(2)*volume/longest**3
    assert np.all(quality>1e-12)
    return dict(passed=True,independentConvexHullMothers=len(z['partialMotherCells']),maxRelativeMotherVolumeError=maximum,
      maxLinearCoordinateError=linear,minLongestEdgeVolumeQuality=float(min(quality)),
      qualityDefinition='6*sqrt(2)*V/maxEdge^3; positive screening only, not grid independence',physicalQualification=False)


def prepare():
    PREP.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False);CASE.mkdir(exist_ok=False)
    region=read(BASE/'region01/report.json');assert region['candidateOneWayOutletPassed']
    save(PREP/'geometry_audit.json',clipped_volume_audit())
    exact=read(OLD/'process_build/report.json')['inputSHA256'];assert sha(ORIGINAL)==exact[str(ORIGINAL)]
    before=ORIGINAL.read_text(encoding='utf8');text='#include "RiemannOutletObserverV240.hpp"\n'+before
    anchor='      /*--- Compute the thermodynamic state in u_b ---*/'
    # 只修改普通BC_Riemann一次；Turbo等其他边界保持原实现。
    start=text.index('void CEulerSolver::BC_Riemann(');end=text.index('void CEulerSolver::BC_TurboRiemann(',start)
    segment=text[start:end];assert segment.count(anchor)==1
    injection='''      if (Marker_Tag == "outlet" && rah_outlet_v240::hook)
        rah_outlet_v240::hook(config->GetInnerIter()+1, geometry->nodes->GetGlobalIndex(iPoint),
                              ProjVelocity_i, SoundSpeed_i, u_i, u_b, Normal);

'''
    text=text[:start]+segment.replace(anchor,injection+anchor)+text[end:]
    modified=PREP/'CEulerSolver.cpp';modified.write_text(text,encoding='utf8')
    (PREP/'boundary_observer.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=str(ORIGINAL),tofile='v240/CEulerSolver.cpp')),encoding='utf8')
    old_driver=NATIVE/'rocketEulerRun.cpp';base=old_driver.read_text(encoding='utf8');driver=base
    driver=driver.replace('#include "rocketEulerFailure.hpp"','#include "rocketEulerFailure.hpp"\n#include "forebodyEulerV240.hpp"')
    driver=driver.replace('tag!="wall"&&tag!="farfield"','tag!="wall"&&tag!="farfield"&&tag!="outlet"')
    driver=driver.replace('g->GetnPoint()!=18027','g->GetnPoint()!=6111')
    line='    observation::observe(geometry[zone][inst][0],solver[zone][inst][0][FLOW_SOL],config[zone]);'
    assert driver.count(line)==1;driver=driver.replace(line,line+'\n    region240::record(geometry[zone][inst][0],solver[zone][inst][0][FLOW_SOL],config[zone]);')
    anchor='\n};\n\nclass RocketEulerDriver'
    monitor='''
  bool Monitor(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,
    CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
    return CFluidIteration::Monitor(o,in,g,s,n,c,sm,gm,b,z,inst)||region240::benefit_stop;
  }
'''
    assert driver.count(anchor)==1;driver=driver.replace(anchor,monitor+anchor)
    driver=driver.replace('    observation::frame(g,f,0,false);','    region240::initialize(g,f,cfg);\n    observation::frame(g,f,0,false);')
    anchor='    std::ofstream result(observation::folder+"/native_completion.json");'
    driver=driver.replace(anchor,'''    const bool forced=output_container[0]->SetResultFiles(g,cfg,solver_container[0][0][0],cfg->GetInnerIter(),true);
    if(!forced)throw std::runtime_error("REGION_FINAL_NATIVE_SAVE");
    rah_outlet_v240::hook=nullptr;
'''+anchor)
    driver=driver.replace('<<observation::last<<",\\\"initialSaved\\\":true', '<<observation::last<<",\\\"benefitStop\\\":"<<(region240::benefit_stop?"true":"false")<<",\\\"forcedFinalNativeOutput\\\":true,\\\"initialSaved\\\":true')
    source=PREP/'forebodyEuler.cpp';source.write_text(driver,encoding='utf8')
    (PREP/'driver.patch').write_text(''.join(difflib.unified_diff(base.splitlines(True),driver.splitlines(True),fromfile=str(old_driver),tofile='v240/forebodyEuler.cpp')),encoding='utf8')
    cmd=read(OLD/'build02/compile_command.json');replacement={str(ORIGINAL):str(modified),str(OLD/'build02/CEulerSolver.obj'):str(BUILD/'boundary.obj'),str(OLD/'build02/CEulerSolver.obj.d'):str(BUILD/'boundary.d')}
    object_cmd=[replacement.get(s,s) for s in cmd]+['-I'+str(ORIGINAL.parent),'-I'+str(NATIVE)]
    original_link=read(ROOT/'outputs/v231_validation/recovery_build03/command.json')
    replacement={str(old_driver):str(source),str(ROOT/'outputs/v231_validation/recovery_build03/rocket.exe'):str(BUILD/'forebody.exe'),str(ROOT/'outputs/v231_validation/recovery_build03/rocket.d'):str(BUILD/'driver.d')}
    link_cmd=[replacement.get(s,s) for s in original_link];link_cmd.insert(link_cmd.index('-Wl,--start-group'),str(BUILD/'boundary.obj'))
    bindings={};dep_commands=[]
    known=read(ROOT/'docs/validation/v231/evidence/actual_dependencies.json');known.update(read(OLD/'build02/dependencies.json'));known.update(exact)
    new=[modified,source,NATIVE/'forebodyEulerV240.hpp',NATIVE/'RiemannOutletObserverV240.hpp']
    for name,command in [('boundary',object_cmd),('driver',link_cmd)]:
        filtered=[];i=0
        while i<len(command):
            if command[i] in ('-MF','-o'):i+=2;continue
            if command[i] in ('-MD','-c','-Wl,--start-group','-Wl,--end-group') or command[i].endswith(('.a','.obj')):i+=1;continue
            filtered.append(command[i]);i+=1
        dep=PREP/(name+'.d');dcmd=filtered+['-M','-MT','dep','-MF',str(dep)];subprocess.check_call(dcmd,cwd=str(OLD_BUILD),timeout=20);dep_commands.append(dcmd)
        for n in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
            if p not in new:assert known[str(p)]==h,str(p)
            bindings[str(p)]=h
    # 归档及替代primitive对象都绑定既有成功构建，不改静态库内容。
    previous=read(ROOT/'outputs/v232_validation/process_build/report.json')['inputSHA256']
    for s in original_link:
        if s.endswith(('.a','.obj')):
            p=Path(s);assert sha(p)==previous[str(p)];bindings[str(p)]=sha(p)
    for name in ('case.cfg','mesh.su2','original_faces.csv'):shutil.copyfile(BASE/'region01'/name,CASE/name)
    save(PREP/'dependencies.json',bindings);save(PREP/'dependency_commands.json',dep_commands)
    save(BUILD/'object_command.json',object_cmd);save(BUILD/'link_command.json',link_cmd)
    save(PREP/'report.json',dict(passed=True,exactOriginalSolver=str(ORIGINAL),originalSolverSHA256=sha(ORIGINAL),
      newBoundarySourceSHA256=sha(modified),newDriverSourceSHA256=sha(source),boundDependencies=len(bindings),
      numericalFormulasChanged=False,readOnlyBoundaryHookOnly=True,matureObserverObjectsToCompile=1,otherMatureObjectsRecompiled=0,
      original800StateUsedForRestart=False,CFDAdmitted=False,physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v240 regional execution closed; follow NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'REGIONAL_PLAN.md',
      NATIVE/'forebodyEulerV240.hpp',NATIVE/'RiemannOutletObserverV240.hpp',NATIVE/'rocketEulerRun.cpp',ORIGINAL]
    if stage=='prepare':
        locks += [BASE/'region01'/s for s in ('mesh.su2','case.cfg','original_faces.csv','provenance.npz','report.json')]
        locks += [OLD/'build02/compile_command.json',OLD/'build02/dependencies.json',OLD/'process_build/report.json',
          ROOT/'outputs/v231_validation/recovery_build03/command.json',ROOT/'docs/validation/v231/evidence/actual_dependencies.json',
          ROOT/'outputs/v232_validation/process_build/report.json',ROOT/'outputs/v232_validation/rocket01/mesh.su2',ROOT/'scripts/phaseEulerRegionalVolume.py']
        return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,256,locks,BASE)
    bindings=read(PREP/'dependencies.json')
    for p,h in bindings.items():assert sha(p)==h,p;locks.append(Path(p))
    locks += [PREP/'dependencies.json',PREP/'report.json',PREP/'geometry_audit.json']
    if stage in ('object','link'):
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        command=BUILD/('object_command.json' if stage=='object' else 'link_command.json');locks.append(command)
        if stage=='link':assert read(BASE/'process_object/report.json')['passed'];locks.append(BUILD/'boundary.obj')
        return guard.supervise(stage,read(command),120-used,256,locks,OLD_BUILD,'build')
    if stage=='run':
        assert read(BASE/'process_link/report.json')['passed'] and not (BASE/'process_run').exists()
        locks += [BUILD/'forebody.exe',BUILD/'boundary.obj']+[p for p in CASE.iterdir() if p.is_file()]
        save(BASE/'regional_admission.json',dict(unique=True,updatesLimit=400,secondsLimit=600,physicalQualification=False,
          geometryAuditSHA256=sha(PREP/'geometry_audit.json'),sourceReportSHA256=sha(PREP/'report.json'),
          executableSHA256=sha(BUILD/'forebody.exe'),normalOutletDynamicNativeGuard=True))
        locks.append(BASE/'regional_admission.json')
        return guard.supervise(stage,[str(BUILD/'forebody.exe'),'case.cfg','.'],600,256,locks,CASE,'cfd',False)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v240 regional execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
