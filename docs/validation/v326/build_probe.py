"""在私有编译目录插入只读测点，链接原有逐节点温界主入口。"""
import argparse
from pathlib import Path
import re
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from buildThermistorWireDriver import arguments,invoke
from shortSphere65 import read,sha,write_json,DRIVER

SPECIES="SU2_CFD/src/solvers/CSpeciesSolver.cpp"
FLOW="SU2_CFD/src/solvers/CIncEulerSolver.cpp"
TEMPLATE="SU2_CFD/include/solvers/CScalarSolver.inl"


def insert(text,anchor,addition,after=False):
    assert text.count(anchor)==1,anchor
    return text.replace(anchor,anchor+addition if after else addition+anchor)


def flow_source(src):
    head,body=src.split("void CIncEulerSolver::Upwind_Residual(",1)
    anchor="if (bounded_scalar) EdgeMassFluxes[iEdge] = conv_residual[0];"
    body=insert(body,anchor,r'''
    if (static_cast<long>(config->GetInnerIter()) == ShortProbe::request().step) {
      const long i = geometry->nodes->GetGlobalIndex(iPoint), j = geometry->nodes->GetGlobalIndex(jPoint);
      if (i == ShortProbe::request().point || j == ShortProbe::request().point)
        ShortProbe::log("flow_edges.csv","%lu,%ld,%ld,%.17g\n",iEdge,i,j,double(conv_residual[0]));
    }
''',True)
    src=head+"void CIncEulerSolver::Upwind_Residual("+body
    head,rest=src.split("void CIncEulerSolver::BC_Far_Field(",1)
    body,tail=rest.split("void CIncEulerSolver::",1)
    body=insert(body,"/*--- Update residual value ---*/",r'''
      if (static_cast<long>(config->GetInnerIter()) == ShortProbe::request().step)
        ShortProbe::log("flow_far.csv","%lu,%.17g\n",geometry->nodes->GetGlobalIndex(iPoint),double(residual[0]));
''')
    return '#include "probe.hpp"\n'+head+"void CIncEulerSolver::BC_Far_Field("+body+"void CIncEulerSolver::"+tail


def scalar_source(src):
    head,rest=src.split("::EdgeFluxResidual(",1)
    edge,rest=rest.split("::BoundaryFluxResidual(",1)
    edge=insert(edge,"      flux.ComputeFlux(opt, iEdge, iPoint, side, jPoint, side, normal, massFlux, updateType,",r'''
      if constexpr (nDim == 3) {
        if (static_cast<long>(config->GetInnerIter()) == ShortProbe::request().step) {
          const long i = geometry->nodes->GetGlobalIndex(iPoint), j = geometry->nodes->GetGlobalIndex(jPoint);
          if (i == ShortProbe::request().point || j == ShortProbe::request().point) {
            auto adv=opt, diff=opt; adv.viscous=false; diff.convective=false;
            const auto a=flux.ComputeFlux(adv,iPoint,side,jPoint,side,normal,massFlux);
            const auto d=flux.ComputeFlux(diff,iPoint,side,jPoint,side,normal,massFlux);
            const auto t=flux.ComputeFlux(opt,iPoint,side,jPoint,side,normal,massFlux);
            ShortProbe::log("scalar_edges.csv","%lu,%ld,%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
              iEdge,i,j,double(normal[0]),double(normal[1]),double(normal[2]),double(massFlux),
              double(nodes->GetSolution(iPoint,0)),double(nodes->GetSolution(jPoint,0)),
              double(a.flux_i(0)),double(a.flux_j(0)),double(d.flux_i(0)),double(d.flux_j(0)),
              double(t.flux_i(0)),double(t.flux_j(0)));
          }
        }
      }
''')
    boundary,rest=rest.split("::FluidInterfaceFluxResidual(",1)
    anchor="    const auto res = flux.ComputeFlux(opt, iPoint, side_i, iVertex, side_j, normal, massFlux);"
    boundary=insert(boundary,anchor,r'''
    if constexpr (nDim == 3) {
      if (static_cast<long>(config->GetInnerIter()) == ShortProbe::request().step &&
          config->GetMarker_All_KindBC(val_marker)==FAR_FIELD) {
        auto* vel=&ghostFlowNodes->GetPrimitive(iVertex)[prim_idx.Velocity()];
        ShortProbe::log("scalar_far.csv","%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
          geometry->nodes->GetGlobalIndex(iPoint),double(massFlux),double(normal[0]),double(normal[1]),double(normal[2]),
          double(vel[0]),double(vel[1]),double(vel[2]),double(flowNodes->GetDensity(iPoint)),
          double(nodes->GetSolution(iPoint,0)),double(ghostNodes->GetSolution(iVertex,0)),double(res.flux_i(0)));
      }
    }
''',True)
    src=head+"::EdgeFluxResidual("+edge+"::BoundaryFluxResidual("+boundary+"::FluidInterfaceFluxResidual("+rest
    head,rest=src.split("::ImplicitEuler_Iteration(",1)
    body,tail=rest.split("::ExplicitEuler_Iteration(",1)
    body=insert(body,"  PrepareImplicitIteration(geometry, solver_container, config);",r'''
  static std::vector<double> old, spatial;
  static std::vector<std::pair<unsigned long,double>> row;
  static bool captured=false;
  const auto step=static_cast<long>(config->GetInnerIter());
  const bool inspect=step==ShortProbe::request().step;
#pragma omp barrier
#pragma omp master
  {
    if(nVar!=1 || nPointDomain!=101122 || omp_get_num_threads()!=4 || !config->GetBounded_Species()) std::abort();
    old.resize(nPointDomain);spatial.resize(nPointDomain);
    for(unsigned long i=0;i<nPointDomain;++i) {
      old[i]=nodes->GetSolution(i,0);spatial[i]=LinSysRes(i,0);
    }
  }
#pragma omp barrier
''')
    body=insert(body,"  auto iter = System.Solve(",r'''
#pragma omp barrier
#pragma omp master
  {
    if(inspect) {
      unsigned long i=0;while(i<nPointDomain && long(geometry->nodes->GetGlobalIndex(i))!=ShortProbe::request().point) ++i;
      if(i==nPointDomain) std::abort();
      row.clear();
      for(unsigned long j=0;j<nPoint;++j) {
        const auto* block=Jacobian.GetBlock(i,j);
        if(block)row.emplace_back(j,double(block[0]));
      }
      ShortProbe::log("row.csv","%ld,%.17g,%.17g,%.17g,%.17g\n",step,spatial[i],double(LinSysRes(i,0)),
                      double(geometry->nodes->GetVolume(i)),double(nodes->GetDelta_Time(i)));
    }
  }
#pragma omp barrier
''')
    body=insert(body,"  CompleteImplicitIteration(geometry, solver_container, config);",r'''
#pragma omp barrier
#pragma omp master
  {
    unsigned long worst=0,highest=0,below=0,above=0;
    for(unsigned long i=0;i<nPointDomain;++i) {
      const double val=nodes->GetSolution(i,0);
      if(val<nodes->GetSolution(worst,0))worst=i;
      if(val>nodes->GetSolution(highest,0))highest=i;
      below+=(val < -1e-12);above+=(val>ShortProbe::request().ceiling+1e-12);
    }
    if(!below && above)worst=highest;
    ShortProbe::log("extrema.csv","%ld,%lu,%.17g,%.17g,%lu,%lu,%.17g,%.17g,%.17g,%.17g\n",
      step,geometry->nodes->GetGlobalIndex(worst),double(nodes->GetSolution(worst,0)),double(nodes->GetSolution(highest,0)),below,above,
      old[worst],spatial[worst],double(LinSysSol(worst,0)),double(nodes->GetUnderRelaxation(worst)));
    if(!captured && (below||above)) {
      captured=true;
      auto* f=std::fopen("first_update.csv","wb");if(!f) std::abort();
      std::fprintf(f,"id,theta_before,theta_after,spatial,delta,relax\n");
      for(unsigned long i=0;i<nPointDomain;++i)
        std::fprintf(f,"%lu,%.17g,%.17g,%.17g,%.17g,%.17g\n",geometry->nodes->GetGlobalIndex(i),
           old[i],double(nodes->GetSolution(i,0)),spatial[i],double(LinSysSol(i,0)),double(nodes->GetUnderRelaxation(i)));
      std::fclose(f);
    }
    if(inspect) for(const auto& entry:row)
      ShortProbe::log("matrix.csv","%lu,%.17g,%.17g\n",geometry->nodes->GetGlobalIndex(entry.first),entry.second,double(LinSysSol(entry.first,0)));
  }
#pragma omp barrier
''',True)
    return head+"::ImplicitEuler_Iteration("+body+"::ExplicitEuler_Iteration("+tail


def build(out,driver):
    out.mkdir(parents=True,exist_ok=False)
    kit=ROOT/"tools/runtime/su2-thermal-build";src=kit/"source";work=kit/"build-v117"
    record=read(driver/"build.json")
    assert sha(kit/"su2.zip")==record["releaseSHA256"] and sha(driver/"wire.exe")==DRIVER["wire"]
    main=driver/"SU2_CFD.obj";sources=(FLOW,SPECIES,TEMPLATE)
    with ZipFile(kit/"su2.zip") as z:
        for rel in sources:
            name,=[n for n in z.namelist() if n.endswith("/"+rel)]
            assert z.read(name)==(src/rel).read_bytes()
    contents={FLOW:flow_source((src/FLOW).read_text("utf8")),TEMPLATE:scalar_source((src/TEMPLATE).read_text("utf8"))}
    contents[SPECIES]='#include "probe.hpp"\n'+(src/SPECIES).read_text("utf8").replace('"../../include/solvers/CScalarSolver.inl"','"CScalarSolver.inl"')
    for rel,code in contents.items():(out/Path(rel).name).write_text(code,encoding="utf8")
    (out/"probe.hpp").write_bytes((HERE/"probe.hpp").read_bytes())
    compiler=kit/"toolchain/w64devkit/bin/g++.exe";assert sha(compiler)==record["compilerSHA256"]
    entries=read(work/"compile_commands.json");objects=[main]
    for rel in (SPECIES,FLOW):
        entry,=[v for v in entries if v["file"].replace("\\","/").endswith("/"+rel)]
        args=arguments(entry["command"]);args[0]=compiler;obj=out/(Path(rel).stem+".obj")
        for flag,path in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=path
        args[-1:]=["-iquote",(src/rel).parent,"-iquote",(src/TEMPLATE).parent,out/Path(rel).name]
        invoke(Path(rel).stem,args,work,out,180,3,bin_directory=compiler.parent);objects.append(obj)
    link=read(driver/"wire_process.json")["args"];exe=out/"probe.exe"
    assert Path(link[2])==driver/"wire.exe" and Path(link[3])==main
    invoke("probe",[compiler,"-o",exe,*objects,*link[4:]],work,out,180,3,bin_directory=compiler.parent)
    write_json(out/"build.json",dict(exe=sha(exe),mainObject=sha(main),sourceArchive=record["releaseSHA256"],
               compiler=sha(compiler),recipe=sha(HERE/"build_probe.py"),sources={p:sha(src/p) for p in sources},
               generated={p.name:sha(p) for p in out.iterdir() if p.suffix in (".hpp",".inl",".cpp",".obj")}))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("out",type=Path);a=p.parse_args()
    build(a.out.resolve(),ROOT/"outputs/v317_validation/build01")
