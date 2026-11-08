"""模块职责：准确源/对象绑定与隔离补丁；不写共享引擎，不启动编译或CFD。"""
from pathlib import Path
import difflib,hashlib,json,os,shutil,subprocess,time
from phaseAwareViscousAir import PhaseAwareViscousAir
from buildSU2IdealAirFirstSample import arguments

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v220_validation'
OLD_SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
OLD_BUILD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
SOURCE=BASE/'source';BUILD=BASE/'build'
BOOST=ROOT/'tools/runtime/cgal/6.0.1-native/boost'
CPP=['Common/src/CConfig.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/solvers/CNSSolver.cpp',
     'SU2_CFD/src/variables/CNSVariable.cpp','SU2_CFD/src/numerics/flow/convection/roe.cpp']
PATCHED=CPP+['Common/include/option_structure.hpp']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')
def once(text,old,new):
    assert text.count(old)==1,(old[:80],text.count(old));return text.replace(old,new)

def patch(name,text):
    if name=='Common/include/option_structure.hpp':
        text=once(text,'  IDEAL_AIR_VARIABLE_CP = 14,','  RAH_PHASE_AIR = 15, // 相态标准Air，内部能量+300000J/kg\n  IDEAL_AIR_VARIABLE_CP = 14,')
        text=once(text,'  MakePair("IDEAL_AIR_VARIABLE_CP", IDEAL_AIR_VARIABLE_CP)',
            '  MakePair("RAH_PHASE_AIR", RAH_PHASE_AIR)\n  MakePair("IDEAL_AIR_VARIABLE_CP", IDEAL_AIR_VARIABLE_CP)')
        for model in ('VISCOSITYMODEL','CONDUCTIVITYMODEL'):
            # 追加显式20，不移动任何旧枚举值。
            start=text.index('enum class '+model+' {');end=text.index('};',start)
            text=text[:end]+'  RAH_PHASE_LINEAR = 20, // 相态物性线性输运\n'+text[end:]
            needle='  MakePair("CANTERA", '+model+'::CANTERA)'
            text=once(text,needle,needle+'\n  MakePair("RAH_PHASE_LINEAR", '+model+'::RAH_PHASE_LINEAR)')
    elif name=='Common/src/CConfig.cpp':
        start=text.index('  if (Kind_FluidModel == IDEAL_AIR_VARIABLE_CP) {')
        end=text.index('\n\n  /*--- Check if CoolProp',start)
        block=text[start:end].replace('IDEAL_AIR_VARIABLE_CP','RAH_PHASE_AIR').replace('TURB_MODEL::SST','TURB_MODEL::SA').replace('NS/SST','NS/SA')
        block=block.replace('Low_Mach_Corr ||','nMarker_FarField || nMarker_Supersonic_Inlet || nMarker_Supersonic_Outlet || nMarker_Inlet || nMarker_Outlet || Low_Mach_Corr ||')
        block=block.replace('Kind_Data_Riemann[i] != STATIC_SUPERSONIC_INFLOW_PT && Kind_Data_Riemann[i] != STATIC_PRESSURE',
                            'Kind_Data_Riemann[i] != STATIC_SUPERSONIC_INFLOW_PT')
        block=block.replace('static PT inlet and static pressure outlet','ambient PT characteristic boundary')
        text=text[:end]+'\n\n  // 新气体只扩独立身份，旧合同原样保留。\n'+block+text[end:]
        needle='    saParsedOptions = ParseSAOptions(SA_Options, nSA_Options, rank);'
        text=once(text,needle,needle+'''
    if (Kind_FluidModel == RAH_PHASE_AIR &&
        (saParsedOptions.version != SA_OPTIONS::NEG || !saParsedOptions.ft2 ||
         saParsedOptions.qcr2000 || saParsedOptions.comp || saParsedOptions.rot || saParsedOptions.bc))
      SU2_MPI::Error("RAH_PHASE_AIR requires canonical SA-neg WITHFT2", CURRENT_FUNCTION);''')
    elif name=='SU2_CFD/src/solvers/CEulerSolver.cpp':
        text=once(text,'#include "../../include/fluid/CSU2IdealAir.hpp"',
            '#include "../../include/fluid/CSU2IdealAir.hpp"\n#include "../../include/fluid/CSU2PhaseAwareAir.hpp"')
        text=once(text,'    case IDEAL_AIR_VARIABLE_CP:\n      config->SetGas_Constant',
            '    case RAH_PHASE_AIR:\n      config->SetGas_Constant(rah_phase_air::R);\n      auxFluidModel = new CSU2PhaseAwareAir();\n      break;\n    case IDEAL_AIR_VARIABLE_CP:\n      config->SetGas_Constant')
        text=once(text,'      case IDEAL_AIR_VARIABLE_CP:\n        FluidModel[thread]',
            '      case RAH_PHASE_AIR:\n        FluidModel[thread] = new CSU2PhaseAwareAir();\n        break;\n      case IDEAL_AIR_VARIABLE_CP:\n        FluidModel[thread]')
        text=once(text,'    case IDEAL_AIR_VARIABLE_CP:\n      ModelTable',
            '    case RAH_PHASE_AIR:\n      ModelTable << "RAH_PHASE_AIR 60..2000K phase domain e+300000";\n      break;\n    case IDEAL_AIR_VARIABLE_CP:\n      ModelTable')
        text=text.replace('config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP',
                          '(config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP || config->GetKind_FluidModel() == RAH_PHASE_AIR)')
    elif name=='SU2_CFD/src/solvers/CNSSolver.cpp':
        text=once(text,'config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP',
                  '(config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP || config->GetKind_FluidModel() == RAH_PHASE_AIR)')
    elif name=='SU2_CFD/src/variables/CNSVariable.cpp':
        text=once(text,'#include "../../include/fluid/CSU2IdealAir.hpp"',
            '#include "../../include/fluid/CSU2IdealAir.hpp"\n#include "../../include/fluid/CSU2PhaseAwareAir.hpp"')
        text=once(text,'dynamic_cast<CSU2IdealAir*>(FluidModel) == nullptr',
                  '(dynamic_cast<CSU2IdealAir*>(FluidModel) == nullptr && dynamic_cast<CSU2PhaseAwareAir*>(FluidModel) == nullptr)')
        text=once(text,'dynamic_cast<CSU2IdealAir*>(FluidModel) != nullptr',
                  '(dynamic_cast<CSU2IdealAir*>(FluidModel) != nullptr || dynamic_cast<CSU2PhaseAwareAir*>(FluidModel) != nullptr)')
    elif name.endswith('/roe.cpp'):
        start=text.index('CUpwGeneralRoe_Flow::ComputeResidual(')
        text=text[:start]+once(text[start:],'  if (RoeSoundSpeed2 <= 0.0) {',
            '  // 新身份非法平均明确拒绝；合法Roe公式及旧模型行为保持。\n'
            '  if (config->GetKind_FluidModel() == RAH_PHASE_AIR &&\n'
            '      (!std::isfinite(RoeSoundSpeed2) || RoeSoundSpeed2 <= 0.0))\n'
            '    throw std::domain_error("RAH_PHASE_ROE_AVERAGE_REJECTED");\n'
            '  if (RoeSoundSpeed2 <= 0.0) {')
    else:raise ValueError(name)
    return text

def main():
    start=time.monotonic();out=BASE/'prepare03';out.mkdir(parents=True,exist_ok=False)
    treepath=ROOT/'outputs/v167_validation/build01/source_tree.json'
    assert sha(treepath)=='3b3683a359b1a3f4f5fc6b6f46ac92b5ed3e3c75f5e45c086de6d87e7fd5c540'
    fulltree=read(treepath)
    tree={n:h for n,h in fulltree.items() if Path(n).relative_to(OLD_SOURCE).as_posix().startswith(('Common/include/','SU2_CFD/include/')) or Path(n).relative_to(OLD_SOURCE).as_posix() in CPP}
    origins={};missing=[]
    versions=['v167_validation/build01','v164_validation/build01','v163_validation/build03',
              'v165_validation/build01','v154_validation/build01','v158_validation/build01','v162_validation/build01']
    for name,digest in tree.items():
        relative=Path(name).relative_to(OLD_SOURCE)
        candidates=[Path(name)]+[ROOT/'outputs'/v/'executed_sources'/relative for v in versions]+[ROOT/'third_party/ideal_air'/relative.name]
        exact=next((p for p in candidates if p.is_file() and sha(p)==digest),None)
        if exact is None:missing.append(str(relative))
        else:origins[relative.as_posix()]=exact
    save(out/'missing.json',missing);assert not missing,missing
    objects=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    assert len(objects)==298
    for name,digest in objects.items():assert sha(OLD_BUILD/name)==digest,name
    SOURCE.mkdir(exist_ok=True);BUILD.mkdir(exist_ok=True)
    # 原externals只读借用；不是复制/修改第三方目录。编译后按依赖文件再绑定实际读取者。
    if (SOURCE/'externals').exists():assert os.path.samefile(str(SOURCE/'externals'),str(OLD_SOURCE/'externals'))
    else:subprocess.check_call(['cmd','/c','mklink','/J',str(SOURCE/'externals'),str(OLD_SOURCE/'externals')],stdout=subprocess.DEVNULL)
    changes={};origin_rows=[]
    for name,p in origins.items():
        dest=SOURCE/name;dest.parent.mkdir(parents=True,exist_ok=True)
        if name in PATCHED:
            old=p.read_text(encoding='utf8');new=patch(name,old)
            if dest.exists():assert dest.read_text(encoding='utf8')==new,name
            else:dest.write_text(new,encoding='utf8')
            changes[name]=dict(originalSHA256=sha(p),patchedSHA256=sha(dest),
                diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name)))
        elif dest.exists():assert sha(dest)==sha(p),name
        else:os.link(str(p),str(dest))
        origin_rows.append(dict(path=name,origin=str(p),sha256=sha(p)))
    air=PhaseAwareViscousAir();rows=[]
    for i,t in enumerate(air.T):rows.append([t,air.Cp[i],air.H[i],air.S[i],air.mu[i],air.kappa[i]])
    include=SOURCE/'SU2_CFD/include/fluid'
    table='// 由固定v206数据生成；原始数据SHA: '+sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json')+'\n#pragma once\n#include <array>\nnamespace rah_phase_air {\nstruct Row { double t,cp,h,s,mu,k; };\ninline constexpr std::array<Row,'+str(len(rows))+'> table = {{\n'
    table+=''.join('{'+','.join(format(float(x),'.17g') for x in row)+'},\n' for row in rows)+'}};\n}\n'
    (include/'PhaseAirTable.hpp').write_text(table,encoding='utf8')
    for name in ('PhaseAwareAir.hpp','CSU2PhaseAwareAir.hpp'):shutil.copyfile(ROOT/'native/su2PhaseAwareAir'/name,include/name)
    entries=read(OLD_BUILD/'compile_commands.json');commands={};selected={}
    for name in CPP:
        row=next(r for r in entries if r['file'].replace('\\','/').endswith('/'+name));tokens=arguments(row['command']);flags=[];i=1
        while i<len(tokens):
            word=tokens[i]
            if word in ('-MD','-c'):i+=1;continue
            if word in ('-MQ','-MF','-o'):i+=2;continue
            if word.endswith('.cpp'):i+=1;continue
            if word.startswith('-I'):
                old=(OLD_BUILD/word[2:]).resolve()
                try:new=SOURCE/old.relative_to(OLD_SOURCE)
                except ValueError:new=old
                word='-I'+str(new)
            flags.append(word);i+=1
        obj=BUILD/(Path(name).stem+'.obj')
        commands[name]=[tokens[0]]+flags+['-I'+str(BOOST),'-o',str(obj),'-c',str(SOURCE/name)]
        selected[row['output'].replace('\\','/')]=str(obj)
    save(out/'commands.json',commands);save(out/'selected_objects.json',selected);save(out/'source_origins.json',origin_rows)
    save(out/'changes.json',changes)
    report=dict(passed=True,sourceFiles=len(tree),historicalFullTreeFiles=len(fulltree),sourceLogicalBytes=sum((SOURCE/n).stat().st_size for n in origins),
        modifiedFiles=len(changes),objectsVerified=len(objects),objectsToCompile=len(selected),tableRows=len(rows),
        boostVersion=108700,boostVersionSHA256=sha(BOOST/'boost/version.hpp'),boostRootSHA256=sha(BOOST/'boost/math/tools/toms748_solve.hpp'),
        sharedSourceModified=False,sharedObjectsModified=False,seconds=time.monotonic()-start,nativeQueries=0,CFDRuns=0)
    save(out/'report.json',report);print(json.dumps(report))

if __name__=='__main__':main()
