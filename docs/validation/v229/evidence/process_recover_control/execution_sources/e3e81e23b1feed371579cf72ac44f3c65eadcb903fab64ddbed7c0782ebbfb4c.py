"""模块职责：只补Euler公开边界状态观测，保留首构建/首失败，不更改数值公式。"""
from pathlib import Path
import difflib,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runPhaseEulerCompatibilityV229 import ROOT,BASE,DOC,SOURCE,BUILD,PREV,OLD_BUILD,HELPER,read,save,sha,locks,size
NEW=BASE/'build02';UNIT='SU2_CFD/src/solvers/CEulerSolver.cpp'
OBSERVED=SOURCE/'SU2_CFD/src/solvers/CEulerSolver_observed.cpp'
BLOCK='''      // RocketAeroHeat：Euler公开观测缓存；只保存已计算状态，不参与任何残差或Jacobian。
      if (config->GetKind_FluidModel() == RAH_PHASE_AIR_EULER) {
        auto observed = GetCharacPrimVar(val_marker, iVertex);
        observed[prim_idx.Temperature()] = Temperature_b;
        for (auto d = 0u; d < nDim; ++d) observed[prim_idx.Velocity() + d] = Velocity_b[d];
        observed[prim_idx.Pressure()] = Pressure_b;
        observed[prim_idx.Density()] = Density_b;
        observed[prim_idx.Enthalpy()] = Enthalpy_b;
        observed[nDim+4] = GetFluidModel()->GetSoundSpeed();
      }

'''

def prepare():
    NEW.mkdir(exist_ok=False);out=BASE/'recovery01';out.mkdir(exist_ok=False)
    original=(SOURCE/UNIT).read_text(encoding='utf8')
    begin=original.index('void CEulerSolver::BC_Riemann(');pos=original.index('      /*--- Compute the residuals ---*/',begin)
    changed=original[:pos]+BLOCK+original[pos:];assert changed.replace(BLOCK,'')==original
    OBSERVED.write_text(changed,encoding='utf8')
    (out/'observation.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='a/'+UNIT,tofile='b/'+UNIT)),encoding='utf8')
    deps=read(BASE/'prepare01/dependencies.json');deps.pop(str(SOURCE/UNIT));deps[str(OBSERVED)]=sha(OBSERVED)
    assert all(sha(n)==h for n,h in deps.items());save(NEW/'dependencies.json',deps)
    cmd=read(BASE/'prepare01/commands.json')[UNIT]
    cmd=[str(OBSERVED) if x==str(SOURCE/UNIT) else x.replace(str(BUILD),str(NEW)) for x in cmd]
    save(NEW/'compile_command.json',cmd)
    save(out/'report.json',dict(passed=True,originalSHA256=sha(SOURCE/UNIT),observedSourceSHA256=sha(OBSERVED),
      newSource=str(OBSERVED),observationOnly=True,unchangedWithoutObservationBlock=True,originalFailurePreserved=True))
    control=BASE/'control02';control.mkdir()
    for n in ('mesh.su2','case.cfg'):shutil.copyfile(BASE/'control01'/n,control/n)

def build():
    subprocess.check_call(read(NEW/'compile_command.json'),cwd=str(OLD_BUILD))
    members=read(BUILD/'members.json');updated={};ar=Path(read(NEW/'compile_command.json')[0]).parent/'ar.exe'
    for group,rows in members.items():
        result=[]
        for row in rows:
            p=Path(row['path']);assert sha(p)==row['sha256']
            if p==BUILD/'CEulerSolver.obj':p=NEW/'CEulerSolver.obj'
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([str(ar),'crsT',str(NEW/group)]+[r['path'] for r in result],cwd=str(OLD_BUILD));updated[group]=result
    save(NEW/'members.json',updated)
    link=[x.replace(str(BUILD),str(NEW)) for x in read(BUILD/'link_command.json')]
    save(NEW/'link_command.json',link);subprocess.check_call(link,cwd=str(OLD_BUILD))
    helper=[x.replace(str(BUILD),str(NEW)) for x in read(BASE/'prepare01/commands.json')['helper']]
    save(NEW/'helper_command.json',helper);subprocess.check_call(helper,cwd=str(OLD_BUILD))
    save(NEW/'report.json',dict(passed=True,engineSHA256=sha(NEW/'SU2_CFD.exe'),helperSHA256=sha(NEW/'compatibility.exe'),
      identity='rah-su2-phase-air-euler-eplus300k-v1',compiled=True,observationOnlyRecovery=True,nativeCompatibility=False,CFDRuns=0))

def main(stage):
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;guard.size_bytes=lambda _:size()
    inputs=locks()+[Path(__file__).resolve()];kind='offline';seconds=120;cwd=BASE
    if stage=='prepare':cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare']
    else:
        inputs += [Path(n) for n in read(NEW/'dependencies.json')]+[NEW/'dependencies.json']
        if stage=='build':kind='build';cmd=[sys.executable,str(Path(__file__).resolve()),'task_build']
        elif stage=='control':
            assert read(NEW/'report.json')['passed'] and read(BASE/'build_audit02/report.json')['nativeControlAdmission']
            kind='query';seconds=30;cwd=BASE/'control02'
            inputs += [NEW/'compatibility.exe',NEW/'SU2_CFD.exe',cwd/'case.cfg',cwd/'mesh.su2']
            cmd=[str(NEW/'compatibility.exe'),'case.cfg','.','control']
        else:raise ValueError('No CFD entry')
    return guard.supervise('recover_'+stage,cmd,seconds,512,list(dict.fromkeys(inputs)),cwd,kind)

if __name__=='__main__':
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
