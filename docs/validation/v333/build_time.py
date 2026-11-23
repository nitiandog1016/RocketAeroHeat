"""在 v332 银球的固定对象上增加 BDF1 储热核对和最后一次不更新的端点装配。"""
import gzip
import json
from pathlib import Path
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/"v329"),str(ROOT/"scripts")]
from build_shared import KIT, DRIVER, once, invoke, read, sha, write_json

BASE=ROOT/"outputs/v332_validation"
BUILD=BASE/"build01";PREP=BASE/"prepared01"
BASE_FRAME="c23c9ca1beb024105c600642f3fb016c9b7077599c8cecff227eb6d3ab0dbe63"
PREP_KEY="22ef3bf27abf1a093fae6bfb03145fac1b9f401229f11a3d12042be19cd9a04c"


def integration(source):
    # 保留原生 BDF1 的装配次序。最后额外一行只读取端点，不计作任何物理或内迭代推进。
    at="  /*--- Add viscous and convective residuals, and compute the Dual Time Source term ---*/"
    source=once(source,at,r"""
  const bool snapshot=MainSolver==SPECIES_SOL && config->GetInnerIter()+1==config->GetnInner_Iter();
  auto captureRow = [&](const char* path) {
    if(!snapshot) return;
#pragma omp barrier
#pragma omp master
    {
      auto* fp=std::fopen(path,"wb"); if(!fp) std::abort();
      for(unsigned long i=0;i<geometry->GetnPointDomain();++i) {
        const double value=solver_container[MainSolver]->LinSysRes(i,0);
        if(std::fwrite(&value,sizeof(double),1,fp)!=1) std::abort();
      }
      std::fclose(fp);
    }
#pragma omp barrier
  };
  captureRow("before_time.bin");
"""+at)
    for anchor,name in (("  /*--- Pick convective and viscous numerics","after_time"),("  /*--- Strong boundary conditions","before_wall"),
                        ("  /*--- Complete residuals for periodic boundary","after_wall")):
        source=once(source,anchor,'  captureRow("'+name+'.bin");\n'+anchor)
    anchor="  switch (config->GetKind_TimeIntScheme())"
    return once(source,anchor,r"""
  if(MainSolver==FLOW_SOL && config->GetInnerIter()+1==config->GetnInner_Iter()) {
    if(geometry->GetnPointDomain()!=261708 || omp_get_num_threads()!=4 ||
       config->GetTime_Marching()!=TIME_MARCHING::DT_STEPPING_1ST) std::abort();
#pragma omp barrier
#pragma omp master
    {
      auto* fp=std::fopen("read_flow.bin","wb");if(!fp) std::abort();
      auto* v=solver_container[MainSolver]->GetNodes();
      for(unsigned long i=0;i<geometry->GetnPointDomain();++i) {
        const double row[]={double(geometry->nodes->GetGlobalIndex(i)),double(v->GetPressure(i)),
          double(v->GetVelocity(i,0)),double(v->GetVelocity(i,1)),double(v->GetVelocity(i,2))};
        if(std::fwrite(row,sizeof(double),5,fp)!=5) std::abort();
      }
      std::fclose(fp);
    }
#pragma omp barrier
    return;
  }
"""+anchor)


def scalar(source):
    anchor="    for(unsigned long i=0;i<nPoint;++i) {\n      auto& n=nodesData[i];"
    source=once(source,anchor,"""    if(dynamic_grid || config->GetTime_Marching()!=TIME_MARCHING::DT_STEPPING_1ST) std::abort();
    endpoint=config->GetInnerIter()+1==config->GetnInner_Iter();
"""+anchor)
    anchor="n.weight=n.volume/n.dt;"
    source=once(source,anchor,anchor+"""
      n.previous=nodes->GetSolution_time_n(i,0);
      const auto* flow=su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
      if(flow->GetDensity_time_n(i)!=1 || flow->GetDensity(i)!=1) std::abort();
      n.storage=(n.old-n.previous)*n.volume/config->GetDelta_UnstTimeND();
""")
    start=source.index("auto predictorIterations=");pos=source.index("  CompleteImplicitIteration(",start)
    return source[:pos]+"  if (!Short65::endpoint)\n  "+source[pos:]


def header(source):
    source=once(source,"  bool wall=false;","  double previous=0,storage=0;\n  bool wall=false;")
    source=once(source,"static double upper","static bool endpoint=false;\nstatic double upper")
    source=once(source,"a.raw+=a.divergence+a.farFlux;","a.raw+=a.divergence+a.farFlux+a.storage;")
    source=once(source,"if(!predicted && step==0 && SilverFar::capture())","if(!predicted && endpoint)")
    source=once(source,"n.budget.give,n.budget.take};","n.budget.give,n.budget.take,n.previous,n.storage};")
    source=once(source,"sizeof(double),16,f)!=16","sizeof(double),18,f)!=18")
    anchor='    std::fclose(f);f=file("edges.bin","wb");'
    source=once(source,anchor,r'''    std::fclose(f);f=file("flow_far.csv","wb");
    for(const auto& v:SilverFar::data().vertex) {
      if(v.status!=2 || v.used!=v.mass) std::abort();
      std::fprintf(f,"%lu,%.17g\n",v.global,v.mass);
    }
'''+anchor)
    source=once(source,"double lo=1e30","double storage=0;\n  double lo=1e30")
    source=once(source,"    far+=n.farFlux;","    storage+=n.storage;far+=n.farFlux;")
    source=once(source,"freeSum+wallSum-far-source","freeSum+wallSum-far-source-storage")
    source=once(source,"freeSpatial,wallSpatial","freeTotal,wallTotal")
    source=once(source,"schemeRMS\\n","schemeRMS,storage\\n")
    source=once(source,'%.17g\\n",\n       step','%.17g,%.17g\\n",\n       step')
    return once(source,"std::sqrt(schemeRms/nodesData.size()));","std::sqrt(schemeRms/nodesData.size()),storage);")


def build(out):
    out.mkdir(parents=True,exist_ok=False)
    rec=read(BUILD/"build.json");wc=KIT/"build-v117"
    assert sha(BUILD/"build.json")=="0bbcde57b7bc0aa2889bf1426613a002c498da980f61e06b67b8ae2db7e94c19"
    assert all(sha(BUILD/name)==h for name,h in rec["files"].items())
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    assert sha(KIT/"su2.zip")==rec["archive"]==cache["sourceArchiveSHA256"]
    key=next(Path(k) for k in cache["cachedBuildInputs"] if Path(k).name=="CIntegration.cpp")
    assert all(sha(KIT/Path(k).relative_to(key.parents[4]))==h for k,h in cache["cachedBuildInputs"].items())
    with ZipFile(KIT/"su2.zip") as z:
        name,=[k for k in z.namelist() if k.endswith("/SU2_CFD/src/integration/CIntegration.cpp")]
        original=z.read(name).decode("utf8").replace("\r\n","\n")
    (out/"CIntegration.cpp").write_text(integration(original),encoding="utf8")
    (out/"CScalarSolver.inl").write_text(scalar((BUILD/"CScalarSolver.inl").read_text("utf8")),encoding="utf8")
    (out/"silver.hpp").write_text(header((BUILD/"silver.hpp").read_text("utf8")),encoding="utf8")
    for name in ("far_mass.hpp","capacity.hpp","CSpeciesSolver.cpp"):(out/name).write_bytes((BUILD/name).read_bytes())
    compiler=KIT/"toolchain/w64devkit/bin/g++.exe";assert sha(compiler)==cache["compilerSHA256"]
    for stem,origin in (("CIntegration",DRIVER),("CSpeciesSolver",BUILD)):
        receipt=read(origin/(stem+"_process.json"));args=receipt["args"]
        assert receipt["exitCode"]==0 and not receipt["livePids"] and Path(args[0])==compiler
        args=[str(out/Path(a).name) if Path(a).parent==origin else a for a in args]
        invoke(stem,args,wc,out,180,3,bin_directory=compiler.parent)
    for stem,base,recipe in (("physical",BUILD,"shared"),("release",DRIVER,"wire")):
        args=read(base/(recipe+"_process.json"))["args"]
        assert sha(Path(args[3]))==rec["main"]
        if stem=="physical":args[4]=str(out/"CSpeciesSolver.obj")
        args[2]=str(out/(stem+".exe"));args.insert(4,str(out/"CIntegration.obj"))
        invoke(stem,args,wc,out,180,3,bin_directory=compiler.parent)
    write_json(out/"build.json",dict(baseline=sha(BUILD/"build.json"),recipe=sha(Path(__file__)),
        files={p.name:sha(p) for p in out.iterdir() if p.suffix in (".hpp",".cpp",".inl",".exe",".obj")}))


if __name__=="__main__":build(Path(sys.argv[1]).resolve())
