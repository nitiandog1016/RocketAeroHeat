"""复用固定版本 SU2 对象缓存，只在独立构建副本中加入质量与能量取样。"""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile

DOC = Path(__file__).resolve().parent
ROOT = DOC.parents[2]
BASE = ROOT/"outputs/v370_validation"
sys.path.insert(0,str(DOC.parent/"v367"))
from build_observer import BIN, BUILD, SRC, TOOL, INT, SIMD, UP, arguments, invoke, replace, save, sha

ROE="SU2_CFD/src/numerics/flow/convection/roe.cpp"
DIFFUSION="SU2_CFD/src/numerics/flow/flow_diffusion.cpp"


def instrument(rel,source):
    if rel==SIMD:
        return replace(source,'"flow/convection/upwind.hpp"','"observed_upwind.hpp"')
    if rel==INT:
        marker="/*--- Compute inviscid residuals ---*/"
        guard="if (MainSolver==FLOW_SOL && iMesh==MESH_0) "
        source=replace(source,marker,guard+"HeatLedger::begin(geometry,solver_container[MainSolver],config);\n  "+marker)
        for marker,name in ((marker,"start"),("/*--- Compute viscous residuals","transport"),
                            ("/*--- Compute source term residuals","viscous"),("/*--- Add viscous and convective residuals","source"),
                            ("/*--- Pick convective and viscous numerics","time"),("/*--- Modification of the system","weak"),
                            ("/*--- Strong boundary conditions","fixed"),("/*--- Complete residuals for periodic","wall"),
                            ("//AD::ResumePreaccumulation","end")):
            source=replace(source,marker,guard+f'HeatLedger::sample("{name}",geometry,solver_container[MainSolver]);\n  '+marker)
        call="solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);"
        source=replace(source,call,"HeatLedger::far_window(true);\n        "+call+"\n        HeatLedger::far_window(false);")
        call="solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);"
        source=replace(source,call,call+'\n      HeatLedger::sample("rhs",geometry,solver_container[MainSolver]);')
    elif rel in (ROE,DIFFUSION):
        if rel==ROE:
            start="CNumerics::ResidualType<> CUpwRoeBase_Flow::ComputeResidual("
            end="CUpwRoe_Flow::CUpwRoe_Flow("
            line="return ResidualType<>(Flux, Jacobian_i, Jacobian_j);"
            probe='''if (HeatLedger::at_far) {
      auto& v=HeatLedger::far[HeatLedger::primitive_ids.at(V_i)];
      if(v[0]!=0) std::abort();
      v[0]=1; v[1]=Normal[0]; v[2]=Normal[1]; v[3]=Flux[0]; v[4]=Flux[3];
      v[5]=V_j[4]; v[6]=V_j[5]; v[7]=V_j[1]*UnitNormal[0]+V_j[2]*UnitNormal[1];
    }'''
        else:
            start="CNumerics::ResidualType<> CAvgGrad_Flow::ComputeResidual("
            end="void CAvgGrad_Flow::SetHeatFluxJacobian("
            line="return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);"
            probe='''if (HeatLedger::at_far) {
      auto id=HeatLedger::primitive_ids.at(V_i);
      if(HeatLedger::visc_seen[id]++) std::abort();
      HeatLedger::far[id][8]=-Proj_Flux_Tensor[0]; HeatLedger::far[id][9]=-Proj_Flux_Tensor[3];
    }'''
        before,tail=source.split(start,1)
        body,after=tail.split(end,1)
        assert body.count(line)==(2 if rel==ROE else 1)
        source=before+start+body.replace(line,probe+"\n    "+line)+end+after
    elif rel==UP:
        marker="Base::viscousTerms(iEdge, iPoint, jPoint, V1st, solution_, vector_ij, geometry,"
        source=replace(source,marker,"const auto massBefore=flux(0),energyBefore=flux(nDim+1);\n    "+marker)
        marker="/*--- Stop preaccumulation. ---*/"
        source=replace(source,marker,'''if constexpr(nDim==2) if(HeatLedger::capture) {
      for(size_t k=0;k<Double::Size;++k) {
        if(!updateMask[k]) continue;
        auto& e=HeatLedger::edges[iEdge[k]];
        if(e[0]>=0) std::abort();
        e={double(geometry.nodes->GetGlobalIndex(iPoint[k])),double(geometry.nodes->GetGlobalIndex(jPoint[k])),
           massBefore[k],energyBefore[k],flux(0)[k],flux(3)[k]};
      }
    }
    '''+marker)
    return '#include "ledger.hpp"\n'+source


def build(out):
    out=out.resolve()
    out.relative_to(BASE)
    old=DOC.parent/"v272/evidence/build.json.gz"
    cache=json.loads(gzip.decompress(old.read_bytes()))
    locks={Path(k):h for k,h in cache["cachedBuildInputs"].items()}
    archive=TOOL/"su2.zip"
    compiler=BIN/"g++.exe"
    assert sha(archive)==cache["sourceArchiveSHA256"] and sha(compiler)==cache["compilerSHA256"]
    assert all(sha(p)==h for p,h in locks.items())
    entries=json.loads((BUILD/"compile_commands.json").read_bytes())
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(DOC/"ledger.hpp",out/"ledger.hpp")
    objects=[]
    with ZipFile(archive) as z:
        for rel in (UP,SIMD,INT,ROE,DIFFUSION):
            member,=[k for k in z.namelist() if k.endswith("/"+rel)]
            src=SRC/rel
            raw=z.read(member)
            assert raw==src.read_bytes()
            unit=out/("observed_upwind.hpp" if rel==UP else src.name)
            unit.write_text(instrument(rel,raw.decode("utf-8-sig")),encoding="utf8")
            if rel==UP: continue
            entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
            cmd=arguments(entry["command"])
            obj=out/(src.stem+".obj")
            for flag,val in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):
                cmd[cmd.index(flag)+1]=str(val)
            cmd[-1:]=["-iquote",str(src.parent),"-iquote",str((SRC/UP).parent),str(unit)]
            invoke(src.stem,cmd,BUILD,out,240,3)
            objects.append(obj)
    prev=json.loads((ROOT/"outputs/v272_validation/build01/link_process.json").read_bytes())
    link=prev["args"]
    assert prev["exitCode"]==0 and Path(link[0])==compiler and link[6]=="-Wl,--allow-shlib-undefined"
    exe=out/"ledger.exe"
    invoke("link",[compiler,"-o",exe,link[3],*objects,*link[6:]],BUILD,out,240,3)
    assert sum(p.stat().st_size for p in out.iterdir())<350*1024**2 and all(sha(p)==h for p,h in locks.items())
    save(out/"build.json",dict(executableSHA256=sha(exe),headerSHA256=sha(DOC/"ledger.hpp"),recipeSHA256=sha(Path(__file__)),
         archiveSHA256=sha(archive),cacheManifestSHA256=sha(old),compilerSHA256=sha(compiler)))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    build(p.parse_args().out)
