"""准确绑定既有SU2二进制与执行源码，判定备用Euler是否可直接采用。"""
from pathlib import Path
import json,shutil,time
from readEilmerNative import sha,need
from prepareHybridRocketEulerNative import save
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.monotonic();out=ROOT/'outputs/v198_validation/su2_review01';out.mkdir(parents=True,exist_ok=False)
    package=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-observed-v1'
    runtime=json.loads((package/'runtime.json').read_text());treefile=ROOT/'outputs/v163_validation/build03/source_tree.json'
    tree=json.loads(treefile.read_text());need(sha(treefile)==runtime['sourceTreeSHA256'],'准确构建树不符')
    need(runtime['sourceCommit']=='7048fd35f56d4bf52e586429d6746a36a10255f8','固定SU2版本不符')
    for item in runtime['files']:need(sha(package/item['path'])==item['sha256'],'原运行包改变')
    root=ROOT/'tools/runtime/su2-thermal-build/source';archive=ROOT/'outputs/v163_validation/build03/executed_sources'
    names=['Common/src/CConfig.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
           'SU2_CFD/src/variables/CEulerVariable.cpp','SU2_CFD/src/solvers/CSolverFactory.cpp',
           'SU2_CFD/include/fluid/CSU2IdealAir.hpp','SU2_CFD/include/fluid/IdealAir.hpp']
    binding={};texts={}
    for name in names:
        source=archive/name if (archive/name).exists() else root/name
        need(sha(source)==tree[str(root/name)],'准确源不匹配: '+name)
        target=out/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        binding[name]=dict(source=str(source),sha256=sha(source));texts[name]=source.read_text(encoding='utf-8')
    config=texts[names[0]];euler=texts[names[1]];var=texts[names[2]]
    block=config.split('if (Kind_FluidModel == IDEAL_AIR_VARIABLE_CP) {',1)[1].split('/*--- Check if CoolProp',1)[0]
    need('(Kind_Solver != MAIN_SOLVER::NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::RANS)' in block,'限制合同改变')
    prim=var.split('bool CEulerVariable::SetPrimVar',1)[1].split('void CEulerVariable::SetSecondaryVar',1)[0]
    need(prim.index('SetTDState_rhoe')<prim.index('if (check_dens') and 'catch' not in prim,'Euler拒绝路径改变')
    far=euler.split('void CEulerSolver::BC_Far_Field(',1)[1].split('void CEulerSolver::BC_Riemann(',1)[0]
    need('sqrt(Gamma*Pressure_Bound/Density_Bound)' in far and 'Gamma_Minus_One' in far,'远场实现改变')
    need('throw std::domain_error(message)' in texts[names[5]],'气体拒绝行为改变')
    r=dict(passed=True,decision='existing adapted SU2 Euler entry rejected; proceed to named OpenFOAM13 contract',
        reasons=['actual CConfig explicitly rejects EULER for IDEAL_AIR_VARIABLE_CP',
                 'CEulerVariable property exception precedes original rollback and is not caught',
                 'BC_Far_Field uses constant-Gamma invariants; not the desired variable-cp ambient contract',
                 'current Riemann whitelist is static PT supersonic inlet or static pressure outlet only'],
        sourceBinding=binding,runtimeSHA256=sha(package/'runtime.json'),sourceTreeSHA256=sha(treefile),
        assessmentScope='exact existing adapted package, not a claim that all SU2 formulations are impossible',
        CFDRuns=0,builds=0,physicalQualification=False,newSolvedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    shutil.copyfile(__file__,out/'executed_entry.py');save(out/'report.json',r)
    print(json.dumps({k:v for k,v in r.items() if k not in ['sourceBinding','reasons']}))

if __name__=='__main__':main()
