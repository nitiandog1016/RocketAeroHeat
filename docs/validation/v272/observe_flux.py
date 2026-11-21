"""在独立目录构建两个只读观测对象，复用 SU2 编译缓存，单步回放 CM-011 低速终场。"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import zipfile

import psutil

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts"))
from buildSU2IdealAirFirstSample import arguments
from runEmbeddedSU2 import sha

TOOL = ROOT/"tools/runtime/su2-thermal-build"
SRC, BUILD = TOOL/"source", TOOL/"build-v117"
BIN = TOOL/"toolchain/w64devkit/bin"
OFFICIAL = ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
BASE = ROOT/"outputs/v272_validation"
FILES = ["SU2_CFD/src/integration/CIntegration.cpp", "SU2_CFD/src/solvers/CIncEulerSolver.cpp"]
RESTART = {"fine02": "50745a4193f291386f8d6b25351ddd2e97f57a0162706b28610ab787a448b56b",
           "outer02": "5726ef71515a0f7c08d3a290e971555a9ef64b093166cb2dbf407381dfb64863"}
GATE = ('std::getenv("RAH_CM011_LOW_SPEED_LEDGER") && config->GetAxisymmetric() && '
        'config->GetKind_Solver()==MAIN_SOLVER::INC_NAVIER_STOKES && config->GetEnergy_Equation() && '
        'geometry->GetnDim()==2 && omp_get_num_threads()==1 && config->GetInnerIter()==0')


def save(path, value):
    with path.open("x", encoding="utf8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def replace(text, old, new):
    if text.count(old) != 1:
        raise ValueError("观测位置与审核版本不一致："+old[:80])
    return text.replace(old, new)


def instrument(rel, original):
    if rel == FILES[0]:
        start, tail = original.split("void CIntegration::Space_Integration(", 1)
        body, rest = tail.split("void CIntegration::Time_Integration(", 1)
        anchor = "  /*--- Compute inviscid residuals ---*/"
        body = replace(body, anchor, '''
  const bool trace = GATE && iMesh == MESH_0 && MainSolver == FLOW_SOL;
  auto snapshot = [&](const char* name) {
    if (!trace) return;
    auto file = std::fopen((std::string("stage_")+name+".csv").c_str(),"wb");
    if (!file) std::abort();
    std::fprintf(file,"id,x,r,mass,energy\\n");
    auto flow = solver_container[MainSolver];
    for(unsigned long i=0; i<geometry->GetnPointDomain(); ++i) {
      auto coord=geometry->nodes->GetCoord(i);
      std::fprintf(file,"%lu,%.17g,%.17g,%.17g,%.17g\\n",geometry->nodes->GetGlobalIndex(i),
        double(coord[0]),double(coord[1]),double(flow->LinSysRes(i,0)),double(flow->LinSysRes(i,3)));
    }
    std::fclose(file);
  };
  snapshot("start");
  if (trace) {
    auto file=std::fopen("input_state.csv","wb");
    if (!file) std::abort();
    std::fprintf(file,"id,x,r,p,u,v,h\\n");
    auto nodes=solver_container[MainSolver]->GetNodes();
    for(unsigned long i=0; i<geometry->GetnPointDomain(); ++i){
      auto c=geometry->nodes->GetCoord(i); auto s=nodes->GetSolution(i);
      std::fprintf(file,"%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\\n",
          geometry->nodes->GetGlobalIndex(i),double(c[0]),double(c[1]),
          double(s[0]),double(s[1]),double(s[2]),double(s[3]));
    }
    std::fclose(file);
  }
'''.replace("GATE", GATE)+anchor)
        for marker, name in [("/*--- Compute viscous residuals ---*/", "convection"),
                             ("/*--- Compute source term residuals ---*/", "diffusion"),
                             ("/*--- Add viscous and convective residuals", "source"),
                             ("/*--- Weak boundary conditions ---*/", "before_far"),
                             ("/*--- Modification of the system", "after_far"),
                             ("/*--- Strong boundary conditions", "before_wall"),
                             ("/*--- Complete residuals for periodic", "after_wall"),
                             ("//AD::ResumePreaccumulation", "end")]:
            body = replace(body, marker, f'snapshot("{name}");\n  '+marker)
        return start+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+rest
    head, tail = original.split("void CIncEulerSolver::BC_Far_Field(", 1)
    body, rest = tail.split("void CIncEulerSolver::BC_Inlet(", 1)
    marker = "  /*--- Loop over all the vertices on this boundary marker ---*/"
    body = replace(body, marker, f'''  std::FILE* traceFile=nullptr;
  if ({GATE} && !pressure_based) {{
    traceFile=std::fopen("far_terms.csv","wb");
    if (!traceFile) std::abort();
    std::fprintf(traceFile,"term,id,x,r,nx,nr,mass,energy,href,theta,u,v\\n");
  }}
'''+marker)
    for key, term, value in [("LinSysRes.AddBlock(iPoint, residual);", "conv", "residual"),
                             ("LinSysRes.SubtractBlock(iPoint, residual_v);", "visc", "residual_v")]:
        # 对流在本函数出现两处；固定算例走人工可压 FDS 的非 pressure_based 分支。
        i = body.rindex(key)
        snippet = f'''
    if (traceFile) {{
      auto c=geometry->nodes->GetCoord(iPoint);
      std::fprintf(traceFile,"{term},%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\\n",
        geometry->nodes->GetGlobalIndex(iPoint),double(c[0]),double(c[1]),double(Normal[0]),double(Normal[1]),
        double({value}[0]),double({value}[3]),double(V_infty[prim_idx.Enthalpy()]),
        double(V_domain[prim_idx.Temperature()]),double(V_domain[prim_idx.Velocity()]),double(V_domain[prim_idx.Velocity()+1]));
    }}
'''
        body = body[:i]+snippet+body[i:]
    body = replace(body, "  END_SU2_OMP_FOR", "  END_SU2_OMP_FOR\n  if(traceFile) std::fclose(traceFile);")
    return head+"void CIncEulerSolver::BC_Far_Field("+body+"void CIncEulerSolver::BC_Inlet("+rest


def invoke(name, args, cwd, out, seconds, memory_gib, trace=False, bin_directory=BIN):
    env = dict(os.environ, OMP_NUM_THREADS="1", PATH=str(bin_directory)+os.pathsep+os.environ["PATH"])
    env.pop("RAH_CM011_LOW_SPEED_LEDGER", None)
    if trace:
        env["RAH_CM011_LOW_SPEED_LEDGER"]="1"
    start = time.monotonic()
    peak, reason, family = 0, None, {}
    args = list(map(str, args))
    with (out/(name+".log")).open("xb") as log:
        child = subprocess.Popen(args, cwd=cwd, stdout=log, stderr=log, env=env,
                                 creationflags=subprocess.CREATE_NO_WINDOW)
        parent = psutil.Process(child.pid)
        family[parent.pid] = parent
        try:
            while child.poll() is None:
                try:
                    for p in parent.children(recursive=True): family[p.pid]=p
                except psutil.NoSuchProcess:
                    pass
                used = 0
                for p in family.values():
                    try: used += p.memory_info().rss
                    except psutil.NoSuchProcess: pass
                peak = max(peak, used)
                if time.monotonic()-start > seconds or peak > memory_gib*1024**3:
                    reason="resource budget"
                    break
                time.sleep(.05)
        finally:
            if child.poll() is None:
                for p in reversed(list(family.values())):
                    try: p.kill()
                    except psutil.NoSuchProcess: pass
            child.wait(timeout=10)
    _, alive = psutil.wait_procs(list(family.values()), timeout=3)
    receipt = dict(args=args, cwd=str(cwd), elapsedSeconds=time.monotonic()-start, exitCode=child.returncode,
                   peakSampledRSSBytes=peak, stopReason=reason, livePids=[p.pid for p in alive])
    save(out/(name+"_process.json"), receipt)
    if child.returncode or reason or alive:
        raise RuntimeError("观测进程异常，见 "+str(out/(name+".log")))


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    archive=TOOL/"su2.zip"
    assert sha(archive)=="c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc"
    entries=json.loads((BUILD/"compile_commands.json").read_bytes())
    locks={}; objects=[]
    with zipfile.ZipFile(archive) as z:
        for rel in FILES:
            path=SRC/rel
            member,=[n for n in z.namelist() if n.endswith("/"+rel)]
            assert z.read(member)==path.read_bytes()
            locks[str(path)]=sha(path)
            changed=out/path.name
            changed.write_text(instrument(rel,path.read_text("utf8")),encoding="utf8")
            old,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
            cmd=arguments(old["command"]);obj=out/(path.stem+".obj")
            # 保留原编译目录和参数，依赖文件、替换源码及输出对象只写本次目录。
            cmd[cmd.index("-MQ")+1]=str(obj)
            cmd[cmd.index("-MF")+1]=str(obj)+".d"
            cmd[cmd.index("-o")+1]=str(obj)
            cmd[-1:]=["-iquote",str(path.parent),str(changed)]
            invoke(path.stem,cmd,BUILD,out,120,3)
            objects.append(obj)
    ninja=(BUILD/"build.ninja").read_text("utf8")
    link=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",ninja,re.M)
    assert link
    # 薄归档及其引用的原对象全程只读，命令行上的两个对象优先满足原符号。
    for path in [*BUILD.rglob("*.obj"), BUILD/"SU2_CFD/src/libSU2core.a", BUILD/"Common/src/libSU2Common.a"]:
        locks[str(path)]=sha(path)
    exe=out/"CM011_observe.exe"
    invoke("link",[BIN/"g++.exe","-o",exe,link[1],*objects,*arguments(link[2])],BUILD,out,60,3)
    assert all(sha(Path(k))==h for k,h in locks.items())
    save(out/"build.json",dict(executableSHA256=sha(exe), sourceArchiveSHA256=sha(archive),
         observedSources={rel:sha(out/Path(rel).name) for rel in FILES}, cachedBuildInputs=locks,
         compilerSHA256=sha(BIN/"g++.exe"), sharedFilesChanged=False))


def sample(out, name, exe, trace):
    parent=ROOT/"outputs/v271_validation"/name
    prior=json.loads((parent/"engine_process.json").read_bytes())
    assert prior["status"]=="terminal" and prior["exitCode"]==0 and not prior["timedOut"] and all(prior["inputsStable"].values())
    assert sha(parent/"restart.csv")==prior["outputSHA256"]["restart.csv"]==RESTART[name]
    for n in ("mesh.su2","case.cfg"):
        assert sha(parent/n)==prior["inputSHA256"][str(parent/n)]
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(parent/"restart.csv",out/"initial_flow.csv")
    shutil.copyfile(parent/"mesh.su2",out/"mesh.su2")
    cfg,n=re.subn(r"^ITER= \d+$","ITER= 1",(parent/"case.cfg").read_text(),flags=re.M);assert n==1
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    locks={str(p):sha(p) for p in [exe,out/"case.cfg",out/"mesh.su2",out/"initial_flow.csv"]}
    invoke("step",[exe,"case.cfg","--threads","1"],out,out,30,1,trace)
    assert all(sha(Path(p))==h for p,h in locks.items())
    assert (out/"restart.csv").is_file() and (out/"far_terms.csv").is_file()==trace
    save(out/"sample.json",dict(source=name,trace=trace,inputSHA256=locks,outputSHA256={
        p.name:sha(p) for p in out.iterdir() if p.is_file() and str(p) not in locks}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=["build","observe","off","official"])
    parser.add_argument("new_directory",type=Path)
    parser.add_argument("--source",choices=RESTART)
    parser.add_argument("--build",type=Path)
    args=parser.parse_args();out=args.new_directory.resolve();out.relative_to(BASE)
    if args.action=="build":build(out)
    else:
        assert args.source
        exe=OFFICIAL if args.action=="official" else args.build.resolve()/"CM011_observe.exe"
        if args.action=="official":assert sha(exe)=="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
        else:assert sha(exe)==json.loads((exe.parent/"build.json").read_bytes())["executableSHA256"]
        sample(out,args.source,exe,args.action=="observe")
