"""在低速 CM-011 三维扇区的一次残差装配中读取边界通量；不改变流场方程。"""
import argparse
import gzip
import json
from pathlib import Path
import re
import shutil
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import BIN, BUILD, FILES, OFFICIAL, SRC, TOOL, arguments, invoke, replace, save, sha

BASE = ROOT/"outputs/v273_validation"


def instrument(original):
    before,tail = original.split("void CIntegration::Space_Integration(",1)
    body,after = tail.split("void CIntegration::Time_Integration(",1)
    anchor = "  /*--- Compute inviscid residuals ---*/"
    body = replace(body,anchor,r'''
  bool trace = std::getenv("RAH_CM011_LOW_SPEED_LEDGER") && iMesh == MESH_0 && MainSolver == FLOW_SOL &&
      !config->GetAxisymmetric() && geometry->GetnDim()==3 && geometry->GetnPointDomain()==12341 &&
      omp_get_num_threads()==1 && config->GetInnerIter()==0 && config->GetEnergy_Equation() &&
      config->GetKind_Solver()==MAIN_SOLVER::INC_NAVIER_STOKES;
  std::FILE* tally = trace ? std::fopen("sector_stages.csv","wb") : nullptr;
  if (trace && !tally) std::abort();
  if (tally) std::fprintf(tally,"stage,mass,enthalpy\n");
  auto snapshot = [&](const char* name) {
    if (!tally) return;
    long double m=0, e=0;
    for(unsigned long i=0;i<geometry->GetnPointDomain();++i) {
      m += double(solver_container[MainSolver]->LinSysRes(i,0));
      e += double(solver_container[MainSolver]->LinSysRes(i,4));
    }
    std::fprintf(tally,"%s,%.17g,%.17g\n",name,double(m),double(e));
  };
  if (tally) {
    auto file = std::fopen("observed_input.csv","wb");
    if (!file) std::abort();
    std::fprintf(file,"id,p,u,v,w,h\n");
    for(unsigned long i=0;i<geometry->GetnPointDomain();++i) {
      auto s=solver_container[MainSolver]->GetNodes()->GetSolution(i);
      std::fprintf(file,"%lu,%.17g,%.17g,%.17g,%.17g,%.17g\n",geometry->nodes->GetGlobalIndex(i),
        double(s[0]),double(s[1]),double(s[2]),double(s[3]),double(s[4]));
    }
    std::fclose(file);
  }
  snapshot("start");
'''+anchor)
    for marker,name in (("/*--- Compute source term residuals","internal"),
                        ("/*--- Weak boundary conditions","before_far"),
                        ("/*--- Modification of the system","after_far"),
                        ("/*--- Strong boundary conditions","before_wall"),
                        ("/*--- Complete residuals for periodic","after_wall")):
        body = replace(body,marker,f'snapshot("{name}");\n  '+marker)
    far = "solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);"
    body = replace(body,far,far+r'''
        if (tally) {
          auto file = std::fopen("far_reference.csv","wb");
          if (!file) std::abort();
          std::fprintf(file,"id,href\n");
          for(unsigned long v=0;v<geometry->nVertex[iMarker];++v) {
            auto i=geometry->vertex[iMarker][v]->GetNode();
            auto prim=solver_container[MainSolver]->GetCharacPrimVar(iMarker,v);
            std::fprintf(file,"%lu,%.17g\n",geometry->nodes->GetGlobalIndex(i),double(prim[geometry->GetnDim()+3]));
          }
          std::fclose(file);
        }''')
    sym = "solver_container[MainSolver]->BC_Sym_Plane(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);"
    body = replace(body,sym,'{ '+sym+' snapshot(config->GetMarker_All_TagBound(iMarker).c_str()); }')
    body = replace(body,"//AD::ResumePreaccumulation", 'snapshot("end");\n  if(tally) std::fclose(tally);\n  //AD::ResumePreaccumulation')
    return before+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+after


def build(out):
    """只编译修改的集成对象，链接 v272 已核对的缓存，原安装和源码保持只读。"""
    out.mkdir(parents=True,exist_ok=False)
    previous = json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    locks = previous["cachedBuildInputs"]
    assert all(sha(Path(n))==h for n,h in locks.items())
    archive = TOOL/"su2.zip"
    assert sha(archive)==previous["sourceArchiveSHA256"]
    path = SRC/FILES[0]
    with zipfile.ZipFile(archive) as zip_file:
        member, = [n for n in zip_file.namelist() if n.endswith("/"+FILES[0])]
        assert zip_file.read(member)==path.read_bytes()
    local = out/path.name
    local.write_text(instrument(path.read_text("utf8")),encoding="utf8")
    entries = json.loads((BUILD/"compile_commands.json").read_bytes())
    entry, = [e for e in entries if e["file"].replace("\\","/").endswith("/"+FILES[0])]
    cmd = arguments(entry["command"])
    obj = out/(path.stem+".obj")
    for option,value in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):
        cmd[cmd.index(option)+1] = str(value)
    cmd[-1:] = ["-iquote",str(path.parent),str(local)]
    invoke(path.stem,cmd,BUILD,out,120,3)
    link = re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",
                     (BUILD/"build.ninja").read_text("utf8"),re.M)
    exe = out/"CM011_sector_observe.exe"
    invoke("link",[BIN/"g++.exe","-o",exe,link[1],obj,*arguments(link[2])],BUILD,out,60,3)
    assert all(sha(Path(n))==h for n,h in locks.items())
    save(out/"build.json",dict(originalSHA256=sha(path),instrumentedSHA256=sha(local),executableSHA256=sha(exe),
         previousBuildReceiptSHA256=sha(ROOT/"docs/validation/v272/evidence/build.json.gz")))


def sample(parent, out, build_folder):
    from readCM011TunnelSector import state_at_frame
    frame = state_at_frame(parent)
    if not frame["normalTerminal"] or not frame["passesResidualAndHeatWindow"]:
        raise ValueError("本轮单步观测要求三维正常终场通过残差与热稳定窗口")
    rec = json.loads((parent/"engine_process.json").read_bytes())
    assert rec["exitCode"]==0 and rec["status"]=="terminal" and not rec["timedOut"] and all(rec["inputsStable"].values())
    assert sha(parent/"restart.csv")==rec["outputSHA256"]["restart.csv"]
    exe = OFFICIAL if build_folder is None else build_folder/"CM011_sector_observe.exe"
    expected = rec["inputSHA256"][str(OFFICIAL)] if build_folder is None else json.loads((build_folder/"build.json").read_bytes())["executableSHA256"]
    assert sha(exe)==expected
    out.mkdir(parents=True,exist_ok=False)
    for n in ("case.cfg","mesh.su2"):
        assert sha(parent/n)==rec["inputSHA256"][str(parent/n)]
    shutil.copyfile(parent/"mesh.su2",out/"mesh.su2")
    shutil.copyfile(parent/"restart.csv",out/"initial_flow.csv")
    config,count = re.subn(r"^ITER= (16000|5000)$","ITER= 1",(parent/"case.cfg").read_text("ascii"),flags=re.M)
    assert count==1
    (out/"case.cfg").write_text(config,encoding="ascii")
    locks = {str(p):sha(p) for p in (exe,out/"mesh.su2",out/"case.cfg",out/"initial_flow.csv")}
    invoke("step",[exe,"case.cfg","--threads","1"],out,out,30,1,trace=build_folder is not None)
    assert all(sha(Path(p))==h for p,h in locks.items())
    save(out/"sample.json",dict(parentRestartSHA256=sha(parent/"restart.csv"),inputSHA256=locks,
         outputSHA256={p.name:sha(p) for p in out.iterdir() if p.is_file() and str(p) not in locks}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("build","sample"))
    parser.add_argument("new_output",type=Path)
    parser.add_argument("--build",type=Path)
    parser.add_argument("--run",choices=("sector01","sector02"),default="sector02")
    args=parser.parse_args()
    out=args.new_output.resolve()
    out.relative_to(BASE)
    if args.action=="build": build(out)
    else: sample(BASE/args.run,out,args.build.resolve() if args.build else None)
