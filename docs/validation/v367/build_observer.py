"""用固定源包和已核对的只读缓存构建局部能量观察对象，不改变共享求解器。"""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(DOC.parent/"v272"))
from observe_flux import TOOL, SRC, BUILD, BIN, arguments, invoke, replace, save, sha

INT="SU2_CFD/src/integration/CIntegration.cpp"
NK="SU2_CFD/src/integration/CNewtonIntegration.cpp"
SIMD="SU2_CFD/include/numerics_simd/CNumericsSIMD.cpp"
UP="SU2_CFD/include/numerics_simd/flow/convection/upwind.hpp"
BASE=ROOT/"outputs/v367_validation"


def instrument(rel,source):
    if rel==INT:
        begin,rest=source.split("void CIntegration::Space_Integration(",1)
        body,end=rest.split("void CIntegration::Time_Integration(",1)
        for anchor,stage in (("/*--- Compute inviscid residuals","start"),("/*--- Compute viscous residuals","transport"),
                ("/*--- Compute source term residuals","afterViscousCall"),("/*--- Add viscous and convective residuals","source"),
                ("/*--- Pick convective and viscous numerics","time"),("//AD::ResumePreaccumulation","boundary")):
            body=replace(body,anchor,f'if (MainSolver==FLOW_SOL && iMesh==MESH_0) TailBudget::sample("{stage}",geometry,solver_container[MainSolver],config);\n  '+anchor)
        source=begin+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+end
    elif rel==NK:
        line="solvers[FLOW_SOL]->PrepareImplicitIteration(geometry, solvers, config);"
        source=replace(source,line,line+'\n  TailBudget::sample("rhs",geometry,solvers[FLOW_SOL],config);')
    elif rel==SIMD:
        return replace(source,'"flow/convection/upwind.hpp"','"observed_upwind.hpp"')
    else:
        # 中心通量从原生 Roe 耗散扣减前取值；总通量从同一 SIMD 调用的黏性计算后取值。
        comment="/*--- Dissipation terms. ---*/"
        source=replace(source,comment,'''if constexpr (nDim==2) if (TailBudget::capturing) {
      for (size_t k=0;k<Double::Size;++k) TailBudget::central[k]=flux(3)[k];
    }
    '''+comment)
        line="Base::viscousTerms(iEdge, iPoint, jPoint, V1st, solution_, vector_ij, geometry,"
        source=replace(source,line,"const auto energyBeforeViscosity=flux(nDim+1);\n    "+line)
        line="/*--- Stop preaccumulation. ---*/"
        source=replace(source,line,'''if constexpr (nDim==2) if (TailBudget::capturing) {
      for (size_t k=0;k<Double::Size;++k) {
        if (!updateMask[k]) continue;
        const auto a=geometry.nodes->GetGlobalIndex(iPoint[k]), b=geometry.nodes->GetGlobalIndex(jPoint[k]);
        if (!TailBudget::selected(a) && !TailBudget::selected(b)) continue;
        auto& e=TailBudget::edges[iEdge[k]];
        if (e.stored) std::abort();
        e={a,b,TailBudget::central[k],energyBeforeViscosity[k],flux(3)[k],true};
      }
    }
    '''+line)
    # 在上游 include 之后引入观察器，复用原生枚举和 OpenMP 配置。
    lines=source.splitlines(keepends=True)
    end=max(j for j,line in enumerate(lines) if line.startswith("#include "))
    lines.insert(end+1,'#include "assembly.hpp"\n')
    return "".join(lines)


def build(out):
    out=out.resolve();out.relative_to(BASE);out.mkdir(parents=True,exist_ok=False)
    baseline=DOC.parent/"v272/evidence/build.json.gz"
    cache=json.loads(gzip.decompress(baseline.read_bytes()))
    locks={Path(k):v for k,v in cache["cachedBuildInputs"].items()}
    assert all(sha(p)==v for p,v in locks.items())
    compiler=BIN/"g++.exe";archive=TOOL/"su2.zip"
    assert sha(compiler)==cache["compilerSHA256"] and sha(archive)==cache["sourceArchiveSHA256"]
    entries=json.loads((BUILD/"compile_commands.json").read_bytes())
    source_sha={};obj=[]
    shutil.copyfile(DOC/"assembly.hpp",out/"assembly.hpp")
    with ZipFile(archive) as z:
        for rel in (UP,SIMD,INT,NK):
            member,=[k for k in z.namelist() if k.endswith("/"+rel)]
            raw=z.read(member);path=SRC/rel
            assert path.read_bytes()==raw
            source_sha[rel]=sha(path)
            dest=out/("observed_upwind.hpp" if rel==UP else path.name)
            dest.write_text(instrument(rel,raw.decode("utf8")),encoding="utf8")
            if rel==UP:continue
            entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
            cmd=arguments(entry["command"]);target=out/(path.stem+".obj");obj.append(target)
            for key,val in (("-MQ",target),("-MF",str(target)+".d"),("-o",target)):
                cmd[cmd.index(key)+1]=str(val)
            cmd[-1:]=["-iquote",str(path.parent),"-iquote",str((SRC/UP).parent),str(dest)]
            invoke(path.stem,cmd,BUILD,out,180,3)
    link=json.loads((ROOT/"outputs/v272_validation/build01/link_process.json").read_bytes())["args"]
    assert Path(link[0])==compiler and link[6]=="-Wl,--allow-shlib-undefined"
    exe=out/"tail_observe.exe"
    invoke("link",[compiler,"-o",exe,link[3],*obj,*link[6:]],BUILD,out,180,3)
    assert all(sha(p)==v for p,v in locks.items())
    size=sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    assert size<300*1024**2
    save(out/"build.json",dict(executableSHA256=sha(exe),cacheManifestSHA256=sha(baseline),cacheInputs=len(locks),
         compilerSHA256=sha(compiler),upstreamArchiveSHA256=sha(archive),upstreamFiles=source_sha,
         observerSHA256=sha(DOC/"assembly.hpp"),recipeSHA256=sha(Path(__file__)),buildBytes=size))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("out",type=Path);build(p.parse_args().out)
