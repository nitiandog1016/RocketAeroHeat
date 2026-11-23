"""构建固定银球工况的共享远场驱动。复用 v328 原方程，仅替换远场质量来源。"""
import gzip
import json
from pathlib import Path
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/"v328"
sys.path.insert(0,str(PREV))
from build import ROOT, SPECIES, SCALAR, adapted, once, arguments, invoke, read, sha, write_json

KIT=ROOT/"tools/runtime/su2-thermal-build"
FLOW="SU2_CFD/src/solvers/CIncEulerSolver.cpp"
DRIVER=ROOT/"outputs/v317_validation/build01"


def flow_source(source):
    head,rest=source.split("void CIncEulerSolver::BC_Far_Field(",1)
    body,tail=rest.split("void CIncEulerSolver::",1)
    anchor="  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)"
    body=once(body,anchor,"""
#pragma omp single
  { SilverFar::begin(geometry,config->GetInnerIter(),val_marker,geometry->nVertex[val_marker]); }
"""+anchor)
    anchor="      /*--- Update residual value ---*/"
    body=once(body,anchor,"      SilverFar::publish(iVertex,geometry->nodes->GetGlobalIndex(iPoint),Normal,residual[0]);\n"+anchor)
    return '#include "far_mass.hpp"\n'+head+"void CIncEulerSolver::BC_Far_Field("+body+"void CIncEulerSolver::"+tail


def scalar_source(source,node_count=101122):
    text=once(adapted(source),"nPointDomain!=101122",f"nPointDomain!={int(node_count)}")
    head,rest=text.split("::BoundaryFluxResidual(",1)
    body,tail=rest.split("::FluidInterfaceFluxResidual(",1)
    start=body.index("      massFlux = BoundedScalarBCFlux(")
    end=body.index("normal.data());",start)+len("normal.data());")
    original=body[start:end]
    body=once(body,original,"""
      if(config->GetMarker_All_KindBC(val_marker)==FAR_FIELD) {
        massFlux=SilverFar::consume(geometry,config->GetInnerIter(),val_marker,iVertex,
                   geometry->nodes->GetGlobalIndex(iPoint),normal.data());
        LinSysRes.AddBlock(iPoint,nodes->GetSolution(iPoint),-massFlux);
        if(opt.implicit) Jacobian.AddVal2Diag(iPoint,-massFlux);
      } else {
"""+original+"\n      }")
    anchor="      Short65::nodesData[iPoint].farMass=massFlux;"
    body=once(body,anchor,anchor+"""
      auto& used=SilverFar::data().vertex[iVertex];
      used.used=massFlux;used.adv=res.flux_i(0);used.inside=nodes->GetSolution(iPoint,0);
      used.outside=ghostNodes->GetSolution(iVertex,0);
""")
    anchor="  END_SU2_OMP_FOR"
    body=once(body,anchor,anchor+"""
#pragma omp barrier
#pragma omp single
  { if(config->GetMarker_All_KindBC(val_marker)==FAR_FIELD) SilverFar::finish(); }
""")
    return head+"::BoundaryFluxResidual("+body+"::FluidInterfaceFluxResidual("+tail


def build(out,prepared=None):
    out.mkdir(parents=True,exist_ok=False)
    node_count=101122;binding={}
    if prepared is not None:
        spec=read(prepared/"prepared.json")
        assert sha(prepared/"mesh.su2")==spec["files"]["mesh.su2"]
        node_count=spec["geometry"]["nodes"]
        binding=dict(preparedSHA256=sha(prepared/"prepared.json"),meshSHA256=spec["files"]["mesh.su2"],nodes=node_count)
    src=KIT/"source";work=KIT/"build-v117";compiler=KIT/"toolchain/w64devkit/bin/g++.exe"
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    old,=[Path(k) for k in cache["cachedBuildInputs"] if Path(k).as_posix().endswith("/source/SU2_CFD/src/integration/CIntegration.cpp")]
    assert all(sha(KIT/Path(k).relative_to(old.parents[4]))==v for k,v in cache["cachedBuildInputs"].items())
    assert sha(compiler)==cache["compilerSHA256"] and sha(KIT/"su2.zip")==cache["sourceArchiveSHA256"]
    with ZipFile(KIT/"su2.zip") as z:
        for rel in (FLOW,SPECIES,SCALAR):
            key,=[k for k in z.namelist() if k.endswith("/"+rel)]
            assert z.read(key)==(src/rel).read_bytes()
    # 端点记录与 v328 相同的原生列顺序；受限、低阶和原空间行各自保留。
    anchor="  double lo=1e30,"
    header=once((PREV/"silver.hpp").read_text("utf8"),anchor,r"""
  if(!predicted && step==0 && SilverFar::capture()) {
    auto* f=file("nodes.bin","wb");
    for(const auto& n:nodesData) {
      const double a[]={double(n.global),n.old,n.low,n.final,double(n.wall),n.volume,n.dt,n.native,n.raw,
                        n.correction,n.limited,n.farMass,n.farFlux,n.divergence,n.budget.give,n.budget.take};
      if(std::fwrite(a,sizeof(double),16,f)!=16) std::abort();
    }
    std::fclose(f);f=file("edges.bin","wb");
    for(const auto& e:edges) {
      const double a[]={double(nodesData[e.i].global),double(nodesData[e.j].global),e.mass,e.high,e.correction,e.alpha};
      if(std::fwrite(a,sizeof(double),6,f)!=6) std::abort();
    }
    std::fclose(f);
  }
"""+anchor)
    (out/"silver.hpp").write_text('#include "far_mass.hpp"\n'+header,encoding="utf8")
    for path in (HERE/"far_mass.hpp",PREV/"capacity.hpp"):(out/path.name).write_bytes(path.read_bytes())
    (out/"CIncEulerSolver.cpp").write_text(flow_source((src/FLOW).read_text("utf8")),encoding="utf8")
    (out/"CScalarSolver.inl").write_text(scalar_source((src/SCALAR).read_text("utf8"),node_count),encoding="utf8")
    species=once((src/SPECIES).read_text("utf8"),'"../../include/solvers/CScalarSolver.inl"','"CScalarSolver.inl"')
    (out/"CSpeciesSolver.cpp").write_text('#include "silver.hpp"\n'+species,encoding="utf8")
    objs=[];entries=read(work/"compile_commands.json")
    for rel in (SPECIES,FLOW):
        entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
        args=arguments(entry["command"]);args[0]=compiler;obj=out/(Path(rel).stem+".obj")
        for flag,path in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=path
        args[-1:]=["-iquote",(src/rel).parent,"-iquote",(src/SCALAR).parent,out/Path(rel).name]
        invoke(obj.stem,args,work,out,180,3,bin_directory=compiler.parent);objs.append(obj)
    link=read(DRIVER/"wire_process.json")["args"];exe=out/"shared.exe";main=DRIVER/"SU2_CFD.obj"
    assert Path(link[3])==main and sha(main)=="6f0079a7d064d790699f24cdf85cbd80ad4f4c5278039f1020b4a722d122cd74"
    invoke("shared",[compiler,"-o",exe,main,*objs,*link[4:]],work,out,180,3,bin_directory=compiler.parent)
    write_json(out/"build.json",dict(exe=sha(exe),main=sha(main),archive=sha(KIT/"su2.zip"),recipe=sha(Path(__file__)),**binding,
                    files={p.name:sha(p) for p in out.iterdir() if p.suffix in (".hpp",".inl",".cpp",".obj")}))


if __name__=="__main__":build(Path(sys.argv[1]).resolve())
