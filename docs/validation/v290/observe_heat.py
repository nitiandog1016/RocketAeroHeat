"""只读记录 6 m/s 热敏珠同一物理步的空间、储存和边界强约束装配。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import sys
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v272")]
from observe_flux import BIN, BUILD, FILES, OFFICIAL, SRC, TOOL, arguments, invoke, replace, save
from runSU2NativeDebugger import observe
from runEmbeddedSU2 import sha

RUN=ROOT/"outputs/v290_validation"
SOURCE=ROOT/"outputs/v289_validation/window01"
FRAME="c8a15e3c275f95711557462562237d6534908b5a7d40062693100978a4c8d4ba"
SEED="ea94eb8a3b668a47e1d3e66e59065e704ab05c57ae752fcd6ab75fa112f42445"
MESH="673151f2681de2a150f7de689bea3096e3b6fa6d9b91ddb80dc635e563fe951c"
CFG="ad34fa368a33b1ed5ba3a424497645a261c1b7e332023e82823534b2ff71d2f1"
ENGINE_SHA="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"


def instrument(original):
    head,tail=original.split("void CIntegration::Space_Integration(",1)
    body,end=tail.split("void CIntegration::Time_Integration(",1)
    marker="  /*--- Compute inviscid residuals ---*/"
    code=r'''
  const auto iteration=config->GetInnerIter();
  const bool snapshot=MainSolver==SPECIES_SOL && iMesh==MESH_0 && dual_time &&
      geometry->GetnPointDomain()==148274 && omp_get_num_threads()==1 &&
      config->GetTimeIter()==1 && (iteration==199 || iteration==200);
  auto observeBead = [&](const char* name, int dump) {
    if (!snapshot) return;
    auto scalar=solver_container[MainSolver];
    auto flow=solver_container[FLOW_SOL]->GetNodes();
    auto nodes=scalar->GetNodes();
    long double sum=0;
    auto f=(iteration==200 && dump) ? std::fopen(dump==1?"before_wall.csv":"after_wall.csv","wb") : nullptr;
    if (iteration==200 && dump && !f) std::abort();
    if (f) std::fprintf(f,dump==1?"id,rhs,theta,theta_old,volume,p,u,v,w,rho,dt\n":"id,rhs\n");
    for (unsigned long i=0;i<geometry->GetnPointDomain();++i) {
      double rhs=double(scalar->LinSysRes(i,0)); sum+=rhs;
      if (!f) continue;
      std::fprintf(f,"%lu,%.17g",geometry->nodes->GetGlobalIndex(i),rhs);
      if (dump==1) std::fprintf(f,",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",
        double(nodes->GetSolution(i,0)),double(nodes->GetSolution_time_n(i,0)),double(geometry->nodes->GetVolume(i)),
        double(flow->GetPressure(i)),double(flow->GetVelocity(i,0)),double(flow->GetVelocity(i,1)),
        double(flow->GetVelocity(i,2)),double(flow->GetDensity(i)),double(config->GetDelta_UnstTimeND()));
      std::fputc('\n',f);
    }
    if (f) std::fclose(f);
    auto l=std::fopen("stages.csv","ab"); if (!l) std::abort();
    std::fprintf(l,"%lu,%s,%.17g\n",(unsigned long)iteration,name,double(sum));std::fclose(l);
  };
  observeBead("start",0);
'''
    body=replace(body,marker,code+marker)
    for anchor,tag,dump in [
        ("/*--- Add viscous and convective residuals","beforeTime",0),
        ("/*--- Pick convective and viscous numerics","afterTime",0),
        ("/*--- Weak boundary conditions","beforeFar",0),
        ("/*--- Modification of the system","afterFar",0),
        ("/*--- Strong boundary conditions","beforeWall",1),
        ("/*--- Complete residuals for periodic","afterWall",2)]:
        body=replace(body,anchor,f'observeBead("{tag}",{dump});\n  '+anchor)
    return head+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+end


def build(out):
    cache=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))
    assert sha(TOOL/"su2.zip")==cache["sourceArchiveSHA256"]
    assert all(sha(Path(k))==v for k,v in cache["cachedBuildInputs"].items())
    source=SRC/FILES[0]
    with zipfile.ZipFile(TOOL/"su2.zip") as z:
        member,=[s for s in z.namelist() if s.endswith("/"+FILES[0])]
        assert source.read_bytes()==z.read(member)
    local=out/source.name
    local.write_text(instrument(source.read_text("utf8")),encoding="utf8")
    entries=json.loads((BUILD/"compile_commands.json").read_bytes())
    entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+FILES[0])]
    cmd=arguments(entry["command"]);obj=out/"CIntegration.obj"
    for key,value in [("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)]:
        cmd[cmd.index(key)+1]=str(value)
    cmd[-1:]=["-iquote",str(source.parent),str(local)]
    invoke("compile",cmd,BUILD,out,120,3)
    link=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",
                   (BUILD/"build.ninja").read_text("utf8"),re.M)
    exe=out/"observe_bead.exe"
    invoke("link",[BIN/"g++.exe","-o",exe,link[1],obj,*arguments(link[2])],BUILD,out,60,3)
    assert all(sha(Path(k))==v for k,v in cache["cachedBuildInputs"].items())
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<300*1024**2
    save(out/"build.json",dict(observerSHA256=sha(exe),observedSourceSHA256=sha(local),
        originalSourceSHA256=sha(source),sourceArchiveSHA256=cache["sourceArchiveSHA256"],
        compilerSHA256=sha(BIN/"g++.exe"),cachedInputsVerified=len(cache["cachedBuildInputs"])))


def sample(stage,out,compiled):
    accepted=HERE.parent/"v289/evidence/first_step_200.json"
    assert json.loads(accepted.read_bytes())["acceptedFirstStep"]
    assert sha(OFFICIAL)==ENGINE_SHA and sha(SOURCE/"case.cfg")==CFG and sha(SOURCE/"mesh.su2")==MESH
    assert sha(SOURCE/"seed_00000.csv")==SEED and sha(SOURCE/"restart_00001.csv")==FRAME
    cfg=(SOURCE/"case.cfg").read_text("ascii")
    paths=[HERE/"PLAN.md",Path(__file__).resolve(),accepted,SOURCE/"case.cfg",ROOT/"scripts/runSU2NativeDebugger.py"]
    if stage=="observe":
        binary=compiled/"observe_bead.exe";receipt=compiled/"build.json"
        assert sha(binary)==json.loads(receipt.read_bytes())["observerSHA256"]
        cfg=replace(cfg,"INNER_ITER= 200\n","INNER_ITER= 201\n")
        paths.extend([receipt,compiled/"CIntegration.cpp"])
        initial=SOURCE/"seed_00000.csv";seconds,limit=360,350
    else:
        binary=OFFICIAL;initial=SOURCE/"restart_00001.csv";seconds,limit=60,180
        cfg=replace(cfg,"INNER_ITER= 200\n","INNER_ITER= 1\n")
        cfg=replace(cfg,"MARKER_WALL_SPECIES= (wall, VALUE, 1.0)","MARKER_WALL_SPECIES= (wall, FLUX, 0.0)")
    os.link(initial,out/"seed_00000.csv");os.link(SOURCE/"mesh.su2",out/"mesh.su2")
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    paths.extend([out/"case.cfg",out/"mesh.su2",out/"seed_00000.csv",binary])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    rec=observe([str(binary),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in paths},seconds,limit*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("build","release","observe"))
    cli.add_argument("output",type=Path)
    cli.add_argument("--build",type=Path)
    args=cli.parse_args();folder=args.output.resolve();folder.relative_to(RUN)
    folder.mkdir(parents=True,exist_ok=False)
    if args.stage=="build":build(folder)
    else:sample(args.stage,folder,args.build.resolve() if args.build else None)
