"""为固定的低速试件构建双求解标量更新，原始 SU2 源码与其他驱动不改写。"""
import argparse
import gzip
import json
from pathlib import Path
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from buildThermistorWireDriver import arguments, invoke
from shortSphere65 import DRIVER, read, sha, write_json

SPECIES="SU2_CFD/src/solvers/CSpeciesSolver.cpp"
SCALAR="SU2_CFD/include/solvers/CScalarSolver.inl"
HEADER="SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp"


def once(text, anchor, replacement):
    if text.count(anchor)!=1:
        raise ValueError(anchor)
    return text.replace(anchor,replacement)


def adapted(text):
    # 原生边回调直接读取正在使用的质量和扩散通量；保留原矩阵构造及边界。
    head, tail=text.split("::EdgeFluxResidual(",1)
    edge, tail=tail.split("::BoundaryFluxResidual(",1)
    start="  bool pausePreacc = false;"
    edge=once(edge,start,"""
#pragma omp single
  {
    if(nPointDomain!=101122 || nPoint!=nPointDomain || !opt.boundedScalar || !opt.viscous ||
       !opt.convective || !opt.implicit || omp_get_num_threads()!=4 || nVar!=1) std::abort();
    Short65::begin(nPoint,geometry->GetnEdge());
  }
"""+start)
    anchor="      flux.ComputeFlux(opt, iEdge, iPoint, side, jPoint, side, normal, massFlux, updateType,"
    edge=once(edge,anchor,"""
      if constexpr (nDim == 3) {
        auto diffOpt=opt;diffOpt.convective=false;
        const auto d=flux.ComputeFlux(diffOpt,iPoint,side,jPoint,side,normal,massFlux);
        const auto t=flux.ComputeFlux(opt,iPoint,side,jPoint,side,normal,massFlux);
        auto& e=Short65::edges[iEdge];
        e.i=iPoint;e.j=jPoint;e.mass=massFlux;e.high=t.flux_i(0);e.alpha=1;
        const double ki=d.jac_ii(0,0);
        if(!std::isfinite(ki) || std::abs(d.flux_i(0)+d.flux_j(0))>1e-15) std::abort();
        e.extraDiagonal=std::max(0.,-ki);
        const double jump=nodes->GetSolution(jPoint,0)-nodes->GetSolution(iPoint,0);
        e.correction=d.flux_i(0)+std::max(0.,ki)*jump;
      }
"""+anchor)
    anchor="    const auto res = flux.ComputeFlux(opt, iPoint, side_i, iVertex, side_j, normal, massFlux);"
    boundary, tail=tail.split("::FluidInterfaceFluxResidual(",1)
    boundary=once(boundary,anchor,anchor+"""
    if(config->GetMarker_All_KindBC(val_marker)==FAR_FIELD) {
      Short65::nodesData[iPoint].farFlux=res.flux_i(0);
      Short65::nodesData[iPoint].farMass=massFlux;
    }
""")
    text=head+"::EdgeFluxResidual("+edge+"::BoundaryFluxResidual("+boundary+"::FluidInterfaceFluxResidual("+tail
    head, tail=text.split("::ImplicitEuler_Iteration(",1)
    body,tail=tail.split("::ExplicitEuler_Iteration(",1)
    prepare="  PrepareImplicitIteration(geometry, solver_container, config);"
    body=once(body,prepare,"""
  using namespace Short65;
#pragma omp single
  {
    for(unsigned long i=0;i<nPoint;++i) {
      auto& n=nodesData[i];n.global=geometry->nodes->GetGlobalIndex(i);n.old=nodes->GetSolution(i,0);
      n.low=n.final=std::numeric_limits<double>::quiet_NaN();n.native=LinSysRes(i,0);
      n.volume=geometry->nodes->GetVolume(i);n.dt=nodes->GetDelta_Time(i);n.weight=n.volume/n.dt;
      if(solver_container[FLOW_SOL]->GetNodes()->GetDensity(i)!=1 || !std::isfinite(n.weight)) std::abort();
    }
    for(unsigned short m=0;m<config->GetnMarker_All();++m) if(config->GetMarker_All_TagBound(m)=="wall") {
      for(unsigned long v=0;v<geometry->nVertex[m];++v) {
        auto& n=nodesData[geometry->vertex[m][v]->GetNode()];n.wall=true;
        upper=std::max(upper,n.old);++wallCount;
      }
    }
    if(!wallCount || std::abs(upper-.4683950617283952)>1e-13) std::abort();
    reduceEdges();
    for(const auto& e:edges) if(e.extraDiagonal>0) {
      if(!nodesData[e.i].wall) {Jacobian.AddVal2Diag(e.i,e.extraDiagonal);Jacobian.GetBlock(e.i,e.j)[0]-=e.extraDiagonal;}
      if(!nodesData[e.j].wall) {Jacobian.AddVal2Diag(e.j,e.extraDiagonal);Jacobian.GetBlock(e.j,e.i)[0]-=e.extraDiagonal;}
    }
    for(unsigned long i=0;i<nPoint;++i) if(!nodesData[i].wall) LinSysRes(i,0)-=nodesData[i].correction;
  }
"""+prepare)
    solve="  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);"
    body=once(body,solve,r"""
  auto predictorIterations=System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
#pragma omp barrier
#pragma omp single
  {
    const auto step=config->GetInnerIter();
    for(unsigned long i=0;i<nPoint;++i) nodesData[i].low=nodesData[i].old+LinSysSol(i,0);
    auto lo=std::min_element(nodesData.begin(),nodesData.end(),[](const Node& a,const Node& b){return a.low<b.low;});
    auto hi=std::max_element(nodesData.begin(),nodesData.end(),[](const Node& a,const Node& b){return a.low<b.low;});
    if(lo->low < -1e-12 || hi->low > upper+1e-12) {
      unsigned long i=(hi->low-upper > -lo->low?hi:lo)-nodesData.begin();
      auto* f=file("predictor_row.csv","wb");
      std::fprintf(f,"id,coefficient,delta\n");
      for(unsigned long j=0;j<nPoint;++j) {
        const auto* a=Jacobian.GetBlock(i,j);
        if(a) std::fprintf(f,"%lu,%.17g,%.17g\n",nodesData[j].global,double(a[0]),double(LinSysSol(j,0)));
      }
      std::fclose(f);
    }
    audit(step,true);
    limit();
    for(unsigned long i=0;i<nPoint;++i) {
      if(!nodesData[i].wall) LinSysRes(i,0)-=nodesData[i].limited;
      LinSysSol(i,0)=0;
    }
  }
  // 原生 RMS/最大残差记录实际送入第二次求解的右端，预测器残差单列。
  SetResToZero();
  su2double rr[MAXNVAR]={}, rm[MAXNVAR]={};unsigned long im[MAXNVAR]={};
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for(unsigned long i=0;i<nPointDomain;++i) ResidualReductions_PerThread(i,0,LinSysRes(i,0),rr,rm,im);
  END_SU2_OMP_FOR
  ResidualReductions_FromAllThreads(geometry,config,rr,rm,im);
""" + solve)
    end="  CompleteImplicitIteration(geometry, solver_container, config);"
    body=once(body,end,end+"""
#pragma omp barrier
#pragma omp single
  {
    for(unsigned long i=0;i<nPoint;++i) nodesData[i].final=nodes->GetSolution(i,0);
    audit(config->GetInnerIter(),false);
  }
""")
    return head+"::ImplicitEuler_Iteration("+body+"::ExplicitEuler_Iteration("+tail


def build(out):
    out.mkdir(parents=True,exist_ok=False)
    kit=ROOT/"tools/runtime/su2-thermal-build";src=kit/"source";work=kit/"build-v117"
    driver=ROOT/"outputs/v317_validation/build01"
    compiler=kit/"toolchain/w64devkit/bin/g++.exe"
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    original,=[Path(x) for x in cache["cachedBuildInputs"] if Path(x).as_posix().endswith("/source/SU2_CFD/src/integration/CIntegration.cpp")]
    assert all(sha(kit/Path(p).relative_to(original.parents[4]))==h for p,h in cache["cachedBuildInputs"].items())
    assert sha(driver/"wire.exe")==DRIVER["wire"] and sha(compiler)==cache["compilerSHA256"]
    assert sha(kit/"su2.zip")==cache["sourceArchiveSHA256"]
    with ZipFile(kit/"su2.zip") as z:
        for name in (SPECIES,SCALAR,HEADER):
            member,=[k for k in z.namelist() if k.endswith("/"+name)]
            assert z.read(member)==(src/name).read_bytes()
    for name in ("silver.hpp","capacity.hpp","affine.cpp"):
        (out/name).write_bytes((HERE/name).read_bytes())
    invoke("affine_build",[compiler,"-std=c++17","-O2",out/"affine.cpp","-o",out/"affine.exe"],out,out,60,1,bin_directory=compiler.parent)
    invoke("affine",[out/"affine.exe"],out,out,60,1,bin_directory=compiler.parent)
    template=out/"CScalarSolver.inl";template.write_text(adapted((src/SCALAR).read_text("utf8")),encoding="utf8")
    species=out/"CSpeciesSolver.cpp"
    species.write_text('#include "silver.hpp"\n'+once((src/SPECIES).read_text("utf8"),'"../../include/solvers/CScalarSolver.inl"','"CScalarSolver.inl"'),encoding="utf8")
    entry,=[e for e in read(work/"compile_commands.json") if e["file"].replace("\\","/").endswith("/"+SPECIES)]
    args=arguments(entry["command"]);args[0]=compiler;obj=out/"CSpeciesSolver.obj"
    for flag,path in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=path
    args[-1:]=["-iquote",(src/SPECIES).parent,"-iquote",(src/SCALAR).parent,species]
    invoke("species",args,work,out,180,3,bin_directory=compiler.parent)
    link=read(driver/"wire_process.json")["args"];main=driver/"SU2_CFD.obj";exe=out/"silver.exe"
    assert sha(main)==read(ROOT/"outputs/v326_validation/build01/build.json")["mainObject"] and Path(link[3])==main
    invoke("silver",[compiler,"-o",exe,main,obj,*link[4:]],work,out,180,3,bin_directory=compiler.parent)
    write_json(out/"build.json",dict(exe=sha(exe),object=sha(obj),mainObject=sha(main),sourceArchive=sha(kit/"su2.zip"),
               recipe=sha(Path(__file__).resolve()),sources={n:sha(out/n) for n in ("silver.hpp","capacity.hpp","affine.cpp","CScalarSolver.inl","CSpeciesSolver.cpp")}))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("out",type=Path);build(p.parse_args().out.resolve())
