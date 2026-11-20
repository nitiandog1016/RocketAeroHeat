"""对同一原生标量强边界装配进行只读观测，验证公开二进制的反力读法。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import sys
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"docs/validation/v272")]
from observe_flux import BIN, BUILD, FILES, OFFICIAL, SRC, TOOL, arguments, invoke, replace
from thermistorBead import CP,D,DT,RHO,U,read,sha,write_json
from thermistorReaction import terminal


def instrument(original,n):
    before,tail=original.split("void CIntegration::Space_Integration(",1)
    body,after=tail.split("void CIntegration::Time_Integration(",1)
    anchor="  /*--- Compute inviscid residuals ---*/"
    code=r'''
  const bool bead = std::getenv("RAH_CM011_LOW_SPEED_LEDGER") && MainSolver==SPECIES_SOL &&
       iMesh==MESH_0 && geometry->GetnDim()==3 && geometry->GetnPointDomain()==NODECOUNT &&
       omp_get_num_threads()==1 && config->GetInnerIter()==0;
  auto tallyBead = [&](const char* stage, bool dump) {
    if (!bead) return;
    auto scalar=solver_container[MainSolver];
    auto f = dump ? std::fopen("original_scalar.csv","wb") : nullptr;
    if(dump && !f) std::abort();
    if(f) std::fprintf(f,"id,theta,residual\n");
    long double sum=0;
    for(unsigned long i=0;i<geometry->GetnPointDomain();++i) {
      double rhs=double(scalar->LinSysRes(i,0));
      sum+=rhs;
      if(f) std::fprintf(f,"%lu,%.17g,%.17g\n",geometry->nodes->GetGlobalIndex(i),
                         double(scalar->GetNodes()->GetSolution(i,0)),rhs);
    }
    if(f) std::fclose(f);
    auto log=std::fopen("scalar_stages.csv","ab");
    if(!log) std::abort();
    std::fprintf(log,"%s,%.17g\n",stage,double(sum));std::fclose(log);
  };
  tallyBead("start", false);
'''.replace('NODECOUNT',str(n))
    body=replace(body,anchor,code+anchor)
    for marker,name,dump in (("/*--- Weak boundary conditions","beforeFar","false"),
                             ("/*--- Modification of the system","afterFar","false"),
                             ("/*--- Strong boundary conditions","beforeWall","true"),
                             ("/*--- Complete residuals for periodic","afterWall","false")):
        body=replace(body,marker,f'tallyBead("{name}", {dump});\n  '+marker)
    return before+"void CIntegration::Space_Integration("+body+"void CIntegration::Time_Integration("+after


def audit(parent, probes, out):
    rec=terminal(parent)
    out.mkdir(parents=True,exist_ok=False)
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    assert all(sha(Path(k))==v for k,v in cache["cachedBuildInputs"].items())
    assert sha(TOOL/"su2.zip")==cache["sourceArchiveSHA256"]
    source=SRC/FILES[0]
    with zipfile.ZipFile(TOOL/"su2.zip") as z:
        member,=[k for k in z.namelist() if k.endswith("/"+FILES[0])]
        assert source.read_bytes()==z.read(member)
    initial=np.loadtxt(parent/"restart.csv",delimiter=",",skiprows=1)
    local=out/source.name
    local.write_text(instrument(source.read_text("utf8"),len(initial)),encoding="utf8")
    entries=read(BUILD/"compile_commands.json")
    entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+FILES[0])]
    cmd=arguments(entry["command"]);obj=out/(source.stem+".obj")
    for k,v in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):
        cmd[cmd.index(k)+1]=str(v)
    cmd[-1:]=["-iquote",str(source.parent),str(local)]
    invoke(source.stem,cmd,BUILD,out,120,3)
    link=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",
                   (BUILD/"build.ninja").read_text("utf8"),re.M)
    exe=out/"thermistor_observe.exe"
    invoke("link",[BIN/"g++.exe","-o",exe,link[1],obj,*arguments(link[2])],BUILD,out,60,3)
    assert all(sha(Path(k))==v for k,v in cache["cachedBuildInputs"].items())
    for tag,binary,trace in (("control",OFFICIAL,False),("observed",exe,True)):
        folder=out/tag;folder.mkdir()
        os.link(parent/"mesh.su2",folder/"mesh.su2");os.link(parent/"restart.csv",folder/"initial_flow.csv")
        cfg,k=re.subn(r"^ITER= \d+$","ITER= 1",(parent/"case.cfg").read_text("ascii"),flags=re.M)
        assert k==1
        (folder/"case.cfg").write_text(cfg,encoding="ascii")
        invoke("step",[binary,"case.cfg","--threads","1"],folder,folder,60,1,trace=trace)
    check(parent,probes,out)


def check(parent,probes,out):
    rec=terminal(parent)
    initial=np.loadtxt(parent/"restart.csv",delimiter=",",skiprows=1)
    outputs=[]
    for tag in ("control","observed"):
        folder=out/tag
        run=read(folder/"step_process.json")
        assert run["exitCode"]==0 and run["stopReason"] is None and not run["livePids"]
        assert sha(folder/"initial_flow.csv")==sha(parent/"restart.csv")
        assert sha(folder/"mesh.su2")==sha(parent/"mesh.su2")
        outputs.append(np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1))
    assert sha(out/"control/case.cfg")==sha(out/"observed/case.cfg")
    assert sha(OFFICIAL)==rec["inputSHA256"][str(OFFICIAL)]
    delta=float(np.max(abs(outputs[0]-outputs[1])));assert delta<2e-7
    sample=np.loadtxt(out/"observed/original_scalar.csv",delimiter=",",skiprows=1)
    sample=sample[np.argsort(sample[:,0])]
    assert np.array_equal(sample[:,0],initial[:,0]) and np.array_equal(sample[:,1],initial[:,-1])
    stages=dict((name,float(v)) for name,v in (line.split(",") for line in (out/"observed/scalar_stages.csv").read_text().splitlines()))
    ids=np.load(parent/"nodes.npz")["ids"]
    wall_ids=np.unique(ids[:,:read(parent/"case.json")["meridianIntervals"]+1])
    probe=probes/"reaction0";terminal(probe)
    assert sha(parent/"restart.csv")==sha(probe/"initial_flow.csv")
    row=np.genfromtxt(probe/"restart.csv",delimiter=",",names=True)["Residual_Species_0"]
    scale=12*D**2*RHO*U*CP*DT
    diff_wall=float(max(abs(row[wall_ids]+sample[wall_ids,2]))*scale)
    wall=(stages["beforeWall"]-stages["afterWall"])*scale
    far=(stages["afterFar"]-stages["beforeFar"])*scale
    assert wall>0 and far>0 and abs(stages["beforeFar"]-stages["start"])*scale<1e-11
    assert diff_wall<1e-11 and abs(wall+row[wall_ids].sum()*scale)<1e-10
    result=dict(inputFrameSHA256=sha(parent/"restart.csv"),sampledThetaExact=True,maxOneUpdateDifference=delta,
                maxWallNodeReactionDifferenceW=diff_wall,wallReleaseW=wall,farOutW=far,
                relativeDifference=(wall-far)/wall,stageSums=stages,observerSHA256=sha(out/"thermistor_observe.exe"))
    assert abs(result["relativeDifference"])<.005
    write_json(out/"energy.json",result)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("parent","probes","new_output"):parser.add_argument(arg,type=Path)
    parser.add_argument("--check-saved",action="store_true")
    args=parser.parse_args()
    (check if args.check_saved else audit)(args.parent.resolve(),args.probes.resolve(),args.new_output.resolve())
