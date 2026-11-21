"""用固定 SU2 8.6.0 标量内核只读观测两份热敏珠终场，共用一个独立可执行文件。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import sys
import zipfile

import numpy as np

from common_layer import BASE, GRIDS, HERE, ROOT, read, save, sha
sys.path.insert(0,str(HERE.parent/"v272"))
from observe_flux import OFFICIAL, TOOL, arguments, invoke, replace

INTEGRATION="SU2_CFD/src/integration/CIntegration.cpp"
SPECIES="SU2_CFD/src/solvers/CSpeciesSolver.cpp"
SCALAR="SU2_CFD/include/solvers/CScalarSolver.inl"
HEADERS=("SU2_CFD/include/numerics/species/species_edge_flux.hpp",
         "SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp")


def stage_source(original):
    head,tail=original.split("void CIntegration::Space_Integration(",1)
    body,end=tail.split("void CIntegration::Time_Integration(",1)
    anchor="  /*--- Compute inviscid residuals ---*/"
    fragment=r'''
  const bool beadLedger = std::getenv("RAH_CM011_LOW_SPEED_LEDGER") && MainSolver==SPECIES_SOL &&
    iMesh==MESH_0 && geometry->GetnDim()==3 && config->GetInnerIter()==0 && omp_get_num_threads()==1;
  auto dumpLedger = [&](const char* filename) {
    if (!beadLedger) return;
    auto stream=std::fopen(filename,"wb");
    if (!stream) std::abort();
    auto solver=solver_container[MainSolver];
    std::fprintf(stream,"id,theta,rhs\n");
    for (unsigned long i=0; i<geometry->GetnPointDomain(); ++i)
      std::fprintf(stream,"%lu,%.17g,%.17g\n",geometry->nodes->GetGlobalIndex(i),
        double(solver->GetNodes()->GetSolution(i,0)),double(solver->LinSysRes(i,0)));
    std::fclose(stream);
  };
  dumpLedger("start.csv");
'''
    body=replace(body,anchor,fragment+anchor)
    for marker,tag in (("/*--- Compute viscous residuals","internal"),("/*--- Weak boundary conditions","beforeFar"),
                       ("/*--- Modification of the system","afterFar"),("/*--- Strong boundary conditions","beforeWall"),
                       ("/*--- Complete residuals for periodic","afterWall"),("//AD::ResumePreaccumulation","end")):
        body=replace(body,marker,'dumpLedger("'+tag+'.csv");\n  '+marker)
    return head+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+end


def edge_source(original):
    anchor="  for (auto color : EdgeColoring) {"
    start,tail=original.split("::EdgeFluxResidual(",1)
    body,tail=tail.split("::BoundaryFluxResidual(",1)
    before=r'''
  std::FILE* beadEdges=nullptr;
  if (std::getenv("RAH_CM011_LOW_SPEED_LEDGER") && nDim==3 && config->GetInnerIter()==0 &&
      omp_get_num_threads()==1) {
    if (nVar!=1 || opt.boundedScalar || !opt.convective || !opt.viscous || opt.muscl || opt.oneSided)
      std::abort();
    beadEdges=std::fopen("edges.csv","wb");
    if(!beadEdges) std::abort();
    std::fprintf(beadEdges,"i,j,nx,ny,nz,conv_i,conv_j,diff_i,diff_j,total_i,total_j\n");
  }
'''
    body=replace(body,anchor,before+anchor)
    marker="      flux.ComputeFlux(opt, iEdge, iPoint, side, jPoint, side, normal, massFlux, updateType,"
    before=r'''
      if constexpr (nDim==3) {
        if (beadEdges) {
          auto adv=opt, dif=opt;
          adv.viscous=false; dif.convective=false;
          const auto a=flux.ComputeFlux(adv,iPoint,side,jPoint,side,normal,massFlux);
          const auto d=flux.ComputeFlux(dif,iPoint,side,jPoint,side,normal,massFlux);
          const auto t=flux.ComputeFlux(opt,iPoint,side,jPoint,side,normal,massFlux);
          std::fprintf(beadEdges,"%lu,%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            geometry->nodes->GetGlobalIndex(iPoint),geometry->nodes->GetGlobalIndex(jPoint),
            double(normal[0]),double(normal[1]),double(normal[2]),double(a.flux_i(0)),double(a.flux_j(0)),
            double(d.flux_i(0)),double(d.flux_j(0)),double(t.flux_i(0)),double(t.flux_j(0)));
        }
      }
'''
    body=replace(body,marker,before+marker)
    marker="AD::ResumePreaccumulation(pausePreacc);"
    body=replace(body,marker,"if (beadEdges) std::fclose(beadEdges);\n  "+marker)
    return start+"::EdgeFluxResidual("+body+"::BoundaryFluxResidual("+tail


def build(out):
    out.mkdir(parents=True,exist_ok=False)
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",TOOL)).resolve()
    builddir=kit/"build-v117";bin_dir=kit/"toolchain/w64devkit/bin"
    cache=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))
    old,=[Path(p) for p in cache["cachedBuildInputs"] if Path(p).as_posix().endswith("/source/"+INTEGRATION)]
    old_kit=old.parents[len(Path(INTEGRATION).parts)]
    locked={kit/Path(k).relative_to(old_kit):h for k,h in cache["cachedBuildInputs"].items()}
    assert all(sha(p)==h for p,h in locked.items())
    assert sha(kit/"su2.zip")==cache["sourceArchiveSHA256"] and sha(bin_dir/"g++.exe")==cache["compilerSHA256"]
    source_root=kit/"source"; paths=(INTEGRATION,SPECIES,SCALAR,*HEADERS)
    with zipfile.ZipFile(kit/"su2.zip") as z:
        for rel in paths:
            member,=[n for n in z.namelist() if n.endswith("/"+rel)]
            assert z.read(member)==(source_root/rel).read_bytes()
    changed={INTEGRATION:stage_source((source_root/INTEGRATION).read_text("utf8")),
             SCALAR:edge_source((source_root/SCALAR).read_text("utf8")),
             SPECIES:replace((source_root/SPECIES).read_text("utf8"),'"../../include/solvers/CScalarSolver.inl"','"CScalarSolver.inl"')}
    for rel,content in changed.items():
        (out/Path(rel).name).write_text(content,encoding="utf8")
    entries=read(builddir/"compile_commands.json");objects=[]
    for rel in (INTEGRATION,SPECIES):
        source=source_root/rel;entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
        args=arguments(entry["command"]);obj=out/(source.stem+".obj")
        assert Path(args[0]).name=="g++.exe";args[0]=bin_dir/"g++.exe"
        for flag,path in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=path
        args[-1:]=["-iquote",source.parent,"-iquote",(source_root/SCALAR).parent,out/source.name]
        invoke(source.stem,args,builddir,out,180,3,bin_directory=bin_dir)
        objects.append(obj)
    link=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",
                   (builddir/"build.ninja").read_text("utf8"),re.M)
    exe=out/"bead_layers.exe"
    invoke("link",[bin_dir/"g++.exe","-o",exe,link[1],*objects,*arguments(link[2])],builddir,out,60,3,bin_directory=bin_dir)
    assert all(sha(p)==h for p,h in locked.items())
    save(out/"build.json",dict(executableSHA256=sha(exe),compilerSHA256=cache["compilerSHA256"],
         sourceArchiveSHA256=cache["sourceArchiveSHA256"],sourceSHA256={k:sha(source_root/k) for k in paths},
         modifiedSHA256={k:sha(out/Path(k).name) for k in changed},cachedInputs=len(locked),license="LGPL-2.1-or-later"))


def sample(level, geometry, builddir, out):
    parent=GRIDS[level-2];g=read(geometry/"geometry.json")["grids"][level-2]
    assert sha(parent/"mesh.su2")==g["meshSHA256"] and sha(parent/"restart.csv")==g["frameSHA256"]
    assert g["mapSHA256"]==sha(geometry/f"grid{level}.npz")
    out.mkdir(exist_ok=False);exe=builddir/"bead_layers.exe"
    assert sha(exe)==read(builddir/"build.json")["executableSHA256"]
    original=parent/"ledger/control"
    reference=read(original/"step_process.json")
    assert reference["exitCode"]==0 and not(reference["livePids"] or reference["stopReason"]) and sha(OFFICIAL)==read(parent/"engine_process.json")["inputSHA256"][str(OFFICIAL)]
    assert sha(original/"initial_flow.csv")==g["frameSHA256"] and sha(original/"mesh.su2")==g["meshSHA256"]
    for src,name in ((parent/"mesh.su2","mesh.su2"),(parent/"restart.csv","initial_flow.csv"),(original/"case.cfg","case.cfg")):
        os.link(src,out/name)
    assert (out/"case.cfg").read_text("ascii")==re.sub(r"^ITER= 7000$","ITER= 1",(parent/"case.cfg").read_text("ascii"),flags=re.M)
    inputs={str(p):sha(p) for p in (exe,out/"mesh.su2",out/"initial_flow.csv",out/"case.cfg",geometry/"geometry.json",
                                     original/"restart.csv",original/"step_process.json",HERE/"PLAN.md",Path(__file__).resolve())}
    invoke("step",[exe,"case.cfg","--threads","1"],out,out,60,1,trace=True,bin_directory=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",TOOL))/"toolchain/w64devkit/bin")
    assert all(sha(Path(p))==v for p,v in inputs.items())
    difference=float(np.max(abs(np.loadtxt(original/"restart.csv",delimiter=",",skiprows=1)-np.loadtxt(out/"restart.csv",delimiter=",",skiprows=1))))
    assert difference < 2e-7 and sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 150*1024**2
    files=("edges.csv","start.csv","internal.csv","beforeFar.csv","afterFar.csv","beforeWall.csv","afterWall.csv","end.csv","restart.csv")
    save(out/"sample.json",dict(grid=level,maxOfficialStepDifference=difference,inputs=inputs,outputs={n:sha(out/n) for n in files}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("action",choices=("build","grid2","grid3"))
    p.add_argument("output",type=Path);p.add_argument("--build",type=Path);p.add_argument("--geometry",type=Path)
    a=p.parse_args();out=a.output.resolve();out.relative_to(BASE)
    if a.action=="build":build(out)
    else:sample(int(a.action[-1]),a.geometry.resolve(),a.build.resolve(),out)
